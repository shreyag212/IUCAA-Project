import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from tqdm import tqdm



# GR units
pc=3.08*10**(16)
G=6.67*10**(-11)
c=3*10**8
Ms=1.98*10**(30)


def generate(M, Pars, vi= np.sqrt(G*(5*10**5)*Ms/(3* 10*pc))):
    
    """
    Docstring:
    Geneate H11, H12, and H22 pols of waveforms for a given sample Space
    
    Calls get_wave()
    
    M (type: float):                                           Total Mass (primary + secondary) in solar mass units
    Pars (type: pandas Dataframe):                             ri (in pc) and Qi (in log10 value) paramters as columns
    """
    df = pd.DataFrame(columns = ['M', 'ri', 'Qi', 'H11', 'H12', 'H22'])
        
    tot = Pars.shape[0]
    count = 0
    
    for k, z in zip(range(Pars.shape[0]), tqdm (range (Pars.shape[0]-1), desc="Loading...")):
        ri = Pars['ri'].iloc[k]
        Qi = Pars['Qi'].iloc[k]
            
            
            
            
        h11, h12, h22 = get_wave(M/2, M/2, ri, Qi, vi=vi)
            
        key = pd.DataFrame([[M, ri, Pars['Qi'].iloc[k], '-1000', '-1000', '-1000']], columns = df.columns)
    
        df = pd.concat([df, key], ignore_index=True)
            
        lng = df.shape[0]
            
        df.at[lng-1, 'H11'] = h11
        df.at[lng-1, 'H12'] = h12
        df.at[lng-1, 'H22'] = h22
            
    return df

#################################################################################################

###################################################
### Needed for function generate() calculations ###
###################################################
def get_wave(m1, m2, ri, Qi, vi = 'def', Dl = 6, duration = 20, dt = 1/256, acc=[6, 100, 100]):
    
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
    
    phis = np.linspace(0,2*np.pi + 2*np.pi/acc[1], acc[1])
    times = np.arange(-duration/2, duration/2+dt, dt)
    
    Times = times[:int(len(times)//2)+1]
    nc = len(Times)
    
    m1 *= Ms
    m2 *= Ms
    ri *= pc
    
    M=m1+m2;
    Mu=m1*m2/(m1+m2)
    
    Qi = 10**Qi
    
    if vi == 'def':
        vi= np.sqrt(G*(5*10**5)*Ms/(3* 10*pc))
    
    rgm=2*G*m1/c**2
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
    epsilon = (L**2)/(GM*rmin)-1
    ep2 = (epsilon**2-1)
    
    phit = np.zeros(nc)
    
    phismin = np.linspace(0,phi0, acc[1])
    phisplus = np.linspace(phi0*0.95, 2.1*np.pi, acc[1])
    
    for j in range(nc):
        
        t = Times[j]
        
        if t <= 0:
            phis = phismin.copy()
        else:
            phis = phisplus.copy()
        
        # first grating
        eqn = np.zeros(acc[1])
        for k in range(acc[1]):
            eqn[k] = Equation(phis[k], t, L, p, epsilon, ep2, phi0)
            
        if np.isnan(min(eqn)):
            phit[j] = -9999.78
            Flag = False
        
        else:
            indx = np.argmin(eqn)
            Flag = True
        

        if Flag:
            for g in range(acc[0]):
                # second grating
                if indx == 0:
                    indx += 1
                if indx == len(eqn)-1:
                    indx -= 1
                
                phis = np.linspace(phis[int(indx-1)], phis[int(indx+1)], acc[2])
                eqn = np.zeros(acc[2])
                for k in range(acc[2]):
                    eqn[k] = Equation(phis[k], t, L, p, epsilon, ep2, phi0)
                indx = np.argmin(eqn)
            
            phit[j] = phis[indx]
            
    phit2 = 2*phi0 - phit[:-1][::-1]
    
    phit = list(phit) + list(phit2)
    
        
    nc = len(times)
    
    h11 = np.zeros(nc)
    h22 = np.zeros(nc)
    h12 = np.zeros(nc)
    pm = (L**2)/(GM*rmin) - 1
    D = 10**Dl * pc
    
    
    for j in range(nc):
        h11[j], h12[j], h22[j] = H_pols(phit[j], phi0, Mu, GM, D, L, pm, rmin)

    return h11, h12, h22
        


def Equation(phi, t, L, p, epsilon, ep2, phi0):
    
    PHI = phi0 - phi
    
    part1 = 2*np.arctanh( (epsilon-1)*np.tan(0.5*PHI)/np.sqrt(ep2) )/(ep2)**1.5
    part2 = epsilon*np.sin(PHI)/((-1+epsilon)*(1+epsilon)*(1+epsilon*np.cos(PHI)))
    
    EQN = (p**2)*(part1-part2) - L*(t)
    
    return abs(EQN) 

def H_pols(phit, phi0, Mu, GM, D, L, pm, rmin, ret='pols'):
    
    
    factor = 2*G*((GM**2)*Mu)/(D*(L**2)*(c**4))
    PHI = phi0-phit
    
    trigs = 3*np.sin(phi0 - 3*phit) + np.sin(phi0+phit)
    part1 = -6*np.sin(2*phit)*(1+pm*np.cos(PHI))**2
    part2 = (pm**2)*(3*np.sin(2*phit))*(np.sin(PHI))**2
    part3 = -pm*(1+pm*np.cos(PHI))*trigs/2
    
    H12 = factor*(part1+part2+part3)
    
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
    
    
    if type(phit) == list:
        if phit == -9999.78:
            return 0, 0, 0

    return H11, H12, H22



def get_wave_amp(m1, m2, ri, Qi, vi = 'def', Dl = 6):
    
    """
    Docstring:
    Implementation of Sajal's MATHEMATICA code for waveform generation on Python
    
    m1 (type: float):                                          Primary Mass in solar mass units
    m2 (type: float):                                          Secondary Mass in solar mass units
    ri (type: float):                                          initial distance of approach
    Qi (type: float):                                          initial angle of approach
    Dl (type float, Default = 6):                              Luminoscity deiscae (in log10 pc values)
    """
    
    
    m1 *= Ms
    m2 *= Ms
    ri *= pc
    
    M=m1+m2;
    Mu=m1*m2/(m1+m2)
    
    Qi = 10**Qi
    
    if vi == 'def':
        vi= np.sqrt(G*(5*10**5)*Ms/(3* 10*pc))
    
    rgm=2*G*m1/c**2
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
    epsilon = (L**2)/(GM*rmin)-1
    ep2 = (epsilon**2-1)
    pm = (L**2)/(GM*rmin) - 1
    D = 10**Dl * pc
    
    phit = np.linspace(0, 2, 1000000)*np.pi
    
    h11, h12, h22 = H_pols(phit, phi0, Mu, GM, D, L, pm, rmin, ret='amp')
        
    
    amp1 = max(abs ( np.array(h11) - np.array(h22) ))
    
    amp2 = 2*max(abs(h12))

    return amp1, amp2






def get_ht(m1, m2, ri, Qi, T, vi = 'def', Dl = 0, acc=[6, 100, 100]):
    
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
    
    phis = np.linspace(0,2*np.pi + 2*np.pi/acc[1], acc[1])
    
    m1 *= Ms
    m2 *= Ms
    ri *= pc
    
    M=m1+m2;
    Mu=m1*m2/(m1+m2)
    
    Qi = 10**Qi
    
    if vi == 'def':
        vi= np.sqrt(G*(5*10**5)*Ms/(3* 10*pc))
    
    rgm=2*G*m1/c**2
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
    epsilon = (L**2)/(GM*rmin)-1
    ep2 = (epsilon**2-1)
    
    phit = np.zeros(1)
    
    phismin = np.linspace(0,phi0, acc[1])
    phisplus = np.linspace(phi0*0.95, 2.1*np.pi, acc[1])
    
    for j in range(1):
        
        t = T
        
        if t <= 0:
            phis = phismin.copy()
        else:
            phis = phisplus.copy()
        
        # first grating
        eqn = np.zeros(acc[1])
        for k in range(acc[1]):
            eqn[k] = Equation(phis[k], t, L, p, epsilon, ep2, phi0)
            
        if np.isnan(min(eqn)):
            phit[j] = -9999.78
            Flag = False
        
        else:
            indx = np.argmin(eqn)
            Flag = True
        

        if Flag:
            for g in range(acc[0]):
                # second grating
                if indx == 0:
                    indx += 1
                if indx == len(eqn)-1:
                    indx -= 1
                
                phis = np.linspace(phis[int(indx-1)], phis[int(indx+1)], acc[2])
                eqn = np.zeros(acc[2])
                for k in range(acc[2]):
                    eqn[k] = Equation(phis[k], t, L, p, epsilon, ep2, phi0)
                indx = np.argmin(eqn)
            
            phit[j] = phis[indx]
            
    nc = 1
    
    h11 = np.zeros(nc)
    h22 = np.zeros(nc)
    h12 = np.zeros(nc)
    pm = (L**2)/(GM*rmin) - 1
    D = 10**Dl * pc
    
    for j in range(nc):
        h11[j], h12[j], h22[j] = H_pols(phit[j], phi0, Mu, GM, D, L, pm, rmin)

    return h11, h12, h22
