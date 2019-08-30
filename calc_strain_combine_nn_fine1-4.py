#!/usr/bin/env python
import h5py 
import numpy as np
import sys

argvs = sys.argv


for i in range(1,5):
    sfile = argvs[1]+i+".hdf5"
    f = h5py.File(sfile)
    sdata = f["avestmk"]
    if i == 0:
        avestmk=np.zeros_like(sdata)
    else:
        avestmk += sdata

with h5py.File(argvs[1]+'1-4.hdf5', 'w') as hf: 
     hf.create_dataset('avestmk', data=avestmk)

