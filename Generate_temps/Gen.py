import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from pycbc.types.timeseries import TimeSeries

from tempgen import generate
import warnings
warnings.filterwarnings("ignore")

## Path for Datasets
path = ''



M = [1,10,25,50]

pars = pd.read_csv(path + 'Params.csv') ## read params

hd = h5py.File('Samples.hdf5', 'w')
h11 = hd.create_group('H11')
h12 = hd.create_group('H12')
h22 = hd.create_group('H22')

Matr = []
ratr = []
Qatr = []
index = []
count = 0

for m in M:
    df = generate(m, pars)


    for i in range(df.shape[0]):
    
    
        j = df.iloc[i]
    
        H11, H12, H22 = np.array(j['H11']), np.array(j['H12']), np.array(j['H22'])
    
        dset = h11.create_dataset(str(count), data=H11)
        dset = h12.create_dataset(str(count), data=H12)
        dset = h22.create_dataset(str(count), data=H22)
        
        Matr.append(m)
        ratr.append(float(pars.iloc[i]['ri']))
        Qatr.append(float(pars.iloc[i]['Qi']))
                    
        index.append(count)
        count += 1
        
    
hd.attrs['M'] = Matr
hd.attrs['ri'] = ratr
hd.attrs['Qi'] = Qatr
hd.attrs['id'] = index

hd.close()