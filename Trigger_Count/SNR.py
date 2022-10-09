import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
from tqdm import tqdm

from pycbc.types.timeseries import TimeSeries
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.events import ThresholdCluster

import sys
home_dir = r'/home/shrey/Desktop/Hyper/'
sys.path.insert(0, home_dir)

## Path for Dataset
path = ''
sys.path.insert(0, path)
from snrs import SNR

dt = 0.00390625
times = np.arange(5, 12282.99609375+dt, dt)

def get_prop_trig(snr, threshold=4, window=1):
    thresh = ThresholdCluster(snr)
    snrv, idx = thresh.threshold_and_cluster(window=window, threshold=threshold)
    
    trigs = times[idx]
    return trigs, abs(snrv)


dt = 1/256
h = open(path + 'H.txt', 'r')
H = h.read().split('\n')
H = [i for i in H if i != '']
h.close()
H = TimeSeries(H, delta_t=dt)

l = open(path + 'L.txt', 'r')
L = l.read().split('\n')
L = [i for i in L if i != '']
l.close()
L = TimeSeries(L, delta_t=dt)

read = h5py.File(path + 'Samples.hdf5', 'r')
H11 = read['H11']
H12 = read['H12']
H22 = read['H22']
keys = list(H11.keys())
keys = np.array(keys).astype(str)

write = h5py.File(path + 'Triggers.hdf5', 'w')
Hw = write.create_group('H')
Htr = Hw.create_group('trig')
Hsn = Hw.create_group('snr')

Lw = write.create_group('L')
Ltr = Lw.create_group('trig')
Lsn = Lw.create_group('snr')

Hpsd = H.psd(4)
Hpsd = interpolate(Hpsd, H.delta_f)
Hpsd = inverse_spectrum_truncation(Hpsd, int(4 * H.sample_rate), trunc_method='hann')

Lpsd = L.psd(4)
Lpsd = interpolate(Lpsd, L.delta_f)
Lpsd = inverse_spectrum_truncation(Lpsd, int(4 * L.sample_rate), trunc_method='hann')

thrsh = 3.5

for i, z in zip(range(len(keys)), tqdm(range(len(keys)), desc="Loading" )):
    index = keys[i]
    hp = np.array(H11[index]) - np.array(H22[index])
    hp = TimeSeries(hp, delta_t=1/256)
    hc = 2*np.array(H12[index])
    hc = TimeSeries(hc, delta_t=1/256)
    
    snrp = SNR(hp, H, psd = Hpsd, ret_ts=True, method='new')
    snrc = SNR(hc, H, psd = Hpsd, ret_ts=True, method='new')
    snrH = np.sqrt(snrp**2 + snrc**2)
    Htrigs = get_prop_trig(snrH, threshold=thrsh)
    
    snrp = SNR(hp, L, psd = Lpsd, ret_ts=True, method='new')
    snrc = SNR(hc, L, psd = Lpsd, ret_ts=True, method='new')
    snrL = np.sqrt(snrp**2 + snrc**2)
    Ltrigs = get_prop_trig(snrL, threshold=thrsh)
    
    dset = Htr.create_dataset(index, data=Htrigs[0])
    dset.attrs['id'] = int(index)
    dset = Hsn.create_dataset(index, data=Htrigs[1])
    dset.attrs['id'] = int(index)
    
    dset = Ltr.create_dataset(index, data=Ltrigs[0])
    dset.attrs['id'] = int(index)
    dset = Lsn.create_dataset(index, data=Ltrigs[1])
    dset.attrs['id'] = int(index)
    
write.close()
