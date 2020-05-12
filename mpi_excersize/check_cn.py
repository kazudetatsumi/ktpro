#!/usr/bin/env python
import h5py
import numpy as np


f=h5py.File("Cn_serial.hdf5")
Cn_serial = f["Cn"][:,:,:,:]
f=h5py.File("Cn.hdf5")
Cn = f["Cn"][:,:,:,:]
diff = Cn - Cn_serial
print(np.max(np.abs(diff)))
print(np.sum(diff))
