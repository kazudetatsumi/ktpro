#!/usr/bin/env python
import h5py 
import numpy as np
import sys

argvs = sys.argv

i=1
sfile = argvs[1]+str(i)+".hdf5"
f = h5py.File(sfile)
avestmk = np.squeeze(f["avestmk"])

for i in range(2,5):
    sfile = argvs[1]+str(i)+".hdf5"
    f = h5py.File(sfile)
    avestmk += np.squeeze(f["avestmk"])


with h5py.File(argvs[1]+'1-4.hdf5', 'w') as hf: 
     hf.create_dataset('avestmk', data=avestmk)

