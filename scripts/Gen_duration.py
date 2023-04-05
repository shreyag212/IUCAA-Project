import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from pycbc.types.timeseries import TimeSeries
from tqdm import tqdm
import os

from tempgen import get_wave, get_wave_amp, get_ht
import warnings
warnings.filterwarnings("ignore")

import argparse

from scipy.signal import tukey


parser = argparse.ArgumentParser(description ='Inputs')
 
parser.add_argument('--min', required=True, type=int, metavar ='', help ='minimum duration')
parser.add_argument('--max', required=True, type=int, metavar ='', help ='maximum duration')

args = parser.parse_args()


min_dur = args.min
max_dur = args.max



df = pd.read_csv('Durations.csv')
df.drop('Unnamed: 0', axis=1, inplace=True)
df_ = df[(df['dur'] >= min_dur) & (df['dur'] <= max_dur)]

fname = 'Temps-' + str(min_dur)+ '-' + str(max_dur) + 's-512Hz.hdf5'
hd = h5py.File(fname, 'w')
hp = hd.create_group('Hp')
hc = hd.create_group('Hc')

ids = list(df_['id'])
for i, z in zip(ids, tqdm(range(len(ids)))):
    
    j = df_[df_['id'] == i].iloc[0]
    
    M, ri, Qi, dur = j['M'], j['ri'], j['Qi'], j['dur']
    
    H11, H12, H22 = get_wave(M/2, M/2, ri, Qi,  duration=dur, dt=1/512, Dl=0)
    
    Hp = np.array(H11) - np.array(H22)
    Hc = np.array(H12)*2


    dset = hp.create_dataset(str(i), data=Hp)
    dset = hc.create_dataset(str(i), data=Hc)
hd.close()



hd = h5py.File('Temps2.hdf5', 'w')
Hp = hd.create_group('Hp')
Hc = hd.create_group('Hc')

files = [i for i in os.listdir() if i[-4:] == 'hdf5' and i[:5] == 'Temps' and i != 'Temps2.hdf5']


oldkeys = []

for file in files:
    print(file)
    hd2 = h5py.File(file, 'r')
    Hp2 = hd2['Hp']
    Hc2 = hd2['Hc']
    
    newkeys = list(Hp2.keys())
    
    for key in newkeys:
        if key not in oldkeys:
            
            hp = np.array(Hp2[key])
            hc = np.array(Hc2[key])
            
            Hp.create_dataset(key, data=hp)
            Hc.create_dataset(key, data=hc)
            
            oldkeys.append(key)
            
    hd2.close()


    
hd.close()

os.remove('Temps.hdf5')
os.rename('Temps2.hdf5', 'Temps.hdf5')