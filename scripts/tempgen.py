import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from tqdm import tqdm

from pycbc.types.timeseries import TimeSeries



# GR units
pc=3.08*10**(16)
G=6.67*10**(-11)
c=3*10**8
Ms=1.98*10**(30)


#################################################################################################

###################################################
### Needed for function generate() calculations ###
###################################################
def get_wave(m1, m2, ri, Qi, vi = np.sqrt(G*(5*10**5)*Ms/(3* 10*pc)), Dl = 0, duration = 20, dt = 1/10, acc=20, load=False, incor=False, plot_phi = False, times=[False]):
    
    """
    Docstring:
    Implementation of Sajal's MATHEMATICA code for waveform generation on Python
    
    m1 (type: float):                                          Primary Mass in solar mass units
    m2 (type: float):                                          Secondary Mass in solar mass units
    ri (type: float):                                          initial distance of approach
    Qi (type: float):                                          initial angle of approach
    Dl (type float, Default = 6):                              Luminoscity deiscae (in log10 pc values)
    duration (type: float, default = 20):                      time duration
    delta_t (type: float, default = 1/256):                    Sampling rate
    
    
    acc (type: list, default = [6, 100, 9):                    Accaracy settings for root finding in function get_wave()
                                                               increasing any of these value will increase waveform accuracy, but 
                                                               increase the loading time
    """
    
    m1, m2, ri, Qi, vi, M, Mu, GM, L, phi0, rmin, p, ecc = get_orbit(m1, m2, ri, Qi, vi)
    ep2 = (ecc**2-1)
    
    if times[0] == False:
        times = np.arange(-duration/2, duration/2+dt/2, dt)
        times = np.round(times, 10)
    Times = times[:int(len(times)//2)+1]
    
    nc = len(Times)
    phit = np.zeros(nc)
    
    
    C1 = (ecc-1)/np.sqrt(ecc**2-1)
    C2 = ecc*np.sqrt(ecc**2-1)/2
    factor = 2*p**2/(ecc**2-1)**(3/2)/L

    tqd = range(nc)
    if load:
        tqd = tqdm(range(nc))
    

    for j, z in zip(range(nc), tqd):
        
        t = Times[j]
        phi1, phi2 = 0, np.pi
        
        Ltf = t/factor
        
        ## Finding roots with bisection method    
        for itt in range(acc):
                
            PHI = (phi1 + phi2)/2
            
            eq1, eq2, eq = Equation(phi1, phi0, Ltf, C1, C2, ecc), Equation(phi2, phi0, Ltf, C1, C2, ecc), Equation(PHI, phi0, Ltf, C1, C2, ecc)
            
            
            if eq1 == 0 or eq2 == 0 or eq == 0:
                break
            s1, s2, s = eq1/abs(eq1), eq2/abs(eq2), eq/abs(eq)
                
            if s1 == s2:
                break
            if s1 == s and s1 != s2:
                phi1 = PHI
            else:
                phi2 = PHI
            
        eq1, eq2, eq = abs(Equation(phi1, phi0, Ltf, C1, C2, ecc)), abs(Equation(phi1, phi0, Ltf, C1, C2, ecc)), abs(Equation(phi1, phi0, Ltf, C1, C2, ecc))
        if eq1 < eq:
            PHI = phi1
        if eq2 < eq:
            PHI = phi2
            
        if t == 0:
            PHI = phi0
                    
        phit[j] = PHI
        
            
    if plot_phi:
        plt.plot(Times, phit)
        
    phit2 = 2*phi0 - phit[:-1][::-1]
    phit = np.array(list(phit) + list(phit2))
    
    if plot_phi:
        plt.plot(times, phit)
        plt.show()
    
    D = 10**Dl * pc
    hp, hc = H_pols(phit, phi0, Mu, GM, D, L, ecc, rmin)

    return TimeSeries(hp, delta_t=dt), TimeSeries(hc, delta_t=dt)

def get_orbit(m1, m2, ri, Qi, vi = np.sqrt(G*(5*10**5)*Ms/(3* 10*pc))):
    
    m1 *= Ms
    m2 *= Ms
    ri *= pc
    Qi = 10**Qi
    
    M = m1+m2
    Mu = m1*m2/(m1+m2)
    GM = G*M
    L = ri*vi*np.sin(Qi)

    phi0 = np.arctan( (L*vi*np.cos(Qi))/(L**2/ri - GM) )
    keep = True
    while keep:
        if phi0 < 0:
            phi0 += np.pi
        elif phi0 > np.pi:
            phi0 -= np.pi
        if phi0 >= 0 and phi0 <= np.pi:
            keep = False
            
    rmin = (L**2)*np.cos(phi0)/(L**2/ri + GM*(np.cos(phi0) - 1))
    p = (L**2)/GM
    ecc = (L**2)/(GM*rmin)-1
    
    return m1, m2, ri, Qi, vi, M, Mu, GM, L, phi0, rmin, p, ecc

    
def Equation(phi, phi0, Ltf, C1, C2, ecc):
    
    PHI = phi0 - phi
    return np.arctanh( C1* np.tan(PHI/2)) - C2*np.sin(PHI)/(1+ecc*np.cos(PHI)) - Ltf



def H_pols(phit, phi0, Mu, GM, D, L, pm, rmin, ret='pols'):
    
    factor = 2*G*((GM**2)*Mu)/(D*(L**2)*(c**4))
    PHI = phi0-phit

    
    if ret == 'pols' or ret == 'hc':
        trigs = 3*np.sin(phi0 - 3*phit) + np.sin(phi0+phit)
        part1 = -6*np.sin(2*phit)*(1+pm*np.cos(PHI))**2
        part2 = (pm**2)*(3*np.sin(2*phit))*(np.sin(PHI))**2
        part3 = -pm*(1+pm*np.cos(PHI))*trigs/2
    
        H12 = factor*(part1+part2+part3)
    
    if ret == 'pols' or ret == 'hp':
    
        trigs = 9*np.cos(-phi0 + 3*phit) - np.cos(phi0)*np.cos(phit) + 5*np.sin(phi0)*np.sin(phit)
        part1 = -6*np.cos(2*phit)*(1+pm*np.cos(PHI))**2
        part2 = (pm**2)*(1+3*np.cos(2*phit))*(np.sin(PHI))**2
        part3 = pm*(1+pm*np.cos(PHI))*trigs/2
    
        H11 = factor*(part1+part2+part3)
    
        trigs = 9*np.cos(-phi0 + 3*phit) - 5*np.cos(phi0)*np.cos(phit) + np.sin(phi0)*np.sin(phit)
        part1 = 6*np.cos(2*phit)*(1+pm*np.cos(PHI))**2
        part2 = -(pm**2)*(-1+3*np.cos(2*phit))*(np.sin(PHI))**2
        part3 = -pm*(1+pm*np.cos(PHI))*trigs/2
    
        H22 = factor*(part1+part2+part3)
        
    if ret == 'pols':
    
        return np.array(H11) - np.array(H22), 2*np.array(H12)
    if ret == 'hc':
        return 2*np.array(H12)
    else:
        return np.array(H11) - np.array(H22)
    


def get_amp(m1, m2, ri, Qi, vi = np.sqrt(G*(5*10**5)*Ms/(3* 10*pc)), Dl = 0):
    
    
    m1, m2, ri, Qi, vi, M, Mu, GM, L, phi0, rmin, p, ecc = get_orbit(m1, m2, ri, Qi, vi)
    ep2 = (ecc**2-1)
    
    D = 10**Dl * pc
    hp, _ = H_pols(phi0, phi0, Mu, GM, D, L, ecc, rmin)
    
    return abs(hp)






def get_ht(t, m1, m2, ri, Qi, vi = np.sqrt(G*(5*10**5)*Ms/(3* 10*pc)), Dl = 0, acc=100, method='min'):
    
    """
    Docstring:
    Implementation of Sajal's MATHEMATICA code for waveform generation on Python
    
    m1 (type: float):                                          Primary Mass in solar mass units
    m2 (type: float):                                          Secondary Mass in solar mass units
    ri (type: float):                                          initial distance of approach
    Qi (type: float):                                          initial angle of approach
    Dl (type float, Default = 6):                              Luminoscity deiscae (in log10 pc values)
    duration (type: float, default = 20):                      time duration
    delta_t (type: float, default = 1/256):                    Sampling rate
    
    
    acc (type: list, default = [6, 100, 9):                    Accaracy settings for root finding in function get_wave()
                                                               increasing any of these value will increase waveform accuracy, but 
                                                               increase the loading time
    """
    
    m1, m2, ri, Qi, vi, M, Mu, GM, L, phi0, rmin, p, ecc = get_orbit(m1, m2, ri, Qi, vi)
    ep2 = (ecc**2-1)
    
    C1 = (ecc-1)/np.sqrt(ecc**2-1)
    C2 = ecc*np.sqrt(ecc**2-1)/2
    factor = 2*p**2/(ecc**2-1)**(3/2)/L


    phi1, phi2 = 0, np.pi
        
    Ltf = t/factor
        
    if 1==1:
        ## Finding roots with bisection method    
        for itt in range(acc):
                
            PHI = (phi1 + phi2)/2
            
            eq1, eq2, eq = Equation(phi1, phi0, Ltf, C1, C2, ecc), Equation(phi2, phi0, Ltf, C1, C2, ecc), Equation(PHI, phi0, Ltf, C1, C2, ecc)
            
            if eq1 == 0 or eq2 == 0 or eq == 0:
                break
            s1, s2, s = eq1/abs(eq1), eq2/abs(eq2), eq/abs(eq)
                
            if s1 == s2:
                break
            if s1 == s and s1 != s2:
                phi1 = PHI
            else:
                phi2 = PHI
            
        eq1, eq2, eq = abs(Equation(phi1, phi0, Ltf, C1, C2, ecc)), abs(Equation(phi1, phi0, Ltf, C1, C2, ecc)), abs(Equation(phi1, phi0, Ltf, C1, C2, ecc))
        if eq1 < eq:
            PHI = phi1
        if eq2 < eq:
            PHI = phi2
            
        if t == 0:
            PHI = phi0
                    
        phit = PHI
    
    D = 10**Dl * pc
    hp, hc = H_pols(phit, phi0, Mu, GM, D, L, ecc, rmin)

    return hp, hc