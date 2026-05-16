#!/usr/bin/env python
import h5py
import numpy as np


f=h5py.File("Output4D_00_120.hdf5",'r')
data4 = f["data4"][:,:,:,:]
f=h5py.File("Output4D_00_120_brief.hdf5", 'r')
data4_brief = f["data4"][:,:,:,:]
diff = data4 - data4_brief
print(np.max(np.abs(diff)))
print(np.sum(diff))
