#!/usr/bin/env python
import numpy as np
import h5py
from ctypes import *
lib = CDLL("./memory_check.so")


datafile = "/home/kazu/desktop/200312/for_cu/all_filled/16h/tst/eliminated_data.hdf5"
f = h5py.File(datafile)
data = f["data4"][:, :, :, :]
condition = np.array(f["condition"], dtype=np.int32)
size = data.shape

#size = 300
#A = np.zeros((size, size, size, size), dtype=np.float64)
lib.memory_check.restype = c_void_p
lib.memory_check.argtypes = [
                             POINTER(c_int),
                             POINTER(c_int),
                             POINTER(c_int),
                             POINTER(c_int),
                             np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                             np.ctypeslib.ndpointer(dtype=np.int32, ndim=4)
                            ]
lib.memory_check(
                 c_int(size[0]),
                 c_int(size[1]),
                 c_int(size[2]),
                 c_int(size[3]),
                 data,
                 condition
                )

