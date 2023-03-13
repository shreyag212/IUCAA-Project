import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from pycbc.types.timeseries import TimeSeries
from tqdm import tqdm

from tempgen import get_wave
import warnings
warnings.filterwarnings("ignore")


df_ = pd.read_csv('/Dataset/lins.csv') ## read params

hd = h5py.File('linS.hdf5', 'w')
hp = hd.create_group('Hp')
hc = hd.create_group('Hc')



M = [1, 10, 25, 50]

Mt = []
vt = []
rt = []
Qt = []
ct = []


count = 0
for m in M:
    
    df = df_[df_['M'] == m].copy()
    for i, z in zip(range(df.shape[0]), tqdm(range(df.shape[0]), desc = 'Loading M: '+str(m))):
            
            ri, Qi = float(df.iloc[i]['ri']), float(df.iloc[i]['Qi'])
            H11, H12, H22 = get_wave(m/2, m/2, ri, Qi)
            Hp = H11-H22
            Hc = 2*H12
            
            dset = hp.create_dataset(str(count), data=Hp)
            dset = hc.create_dataset(str(count), data=Hc)
            
            Mt.append(m)
            rt.append(ri)
            Qt.append(Qi)
            ct.append(count)
            count += 1
        
    
hd.attrs['M'] = Mt
hd.attrs['Qi'] = Qt
hd.attrs['ri'] = rt
hd.attrs['id'] = ct

hd.close()
