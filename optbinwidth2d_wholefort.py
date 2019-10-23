#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./costfort2d.so")



def calc_cost2d_f90(A, maxw, condition):
    class result(ctypes.Structure):
        _fields_ =[("len0", ctypes.c_int), ("len1", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double)), ("kavearr", ctypes.POINTER(ctypes.c_double)), ("darr", ctypes.POINTER(ctypes.c_double))]

    lib.cost2d.restype = result
    lib.cost2d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)]

    Nmax0 = A.shape[0]
    Nmax1 = A.shape[1]

    result = lib.cost2d(ctypes.byref(ctypes.c_int(maxw[0])), ctypes.byref(ctypes.c_int(maxw[1])), ctypes.byref(ctypes.c_int(Nmax0)), ctypes.byref(ctypes.c_int(Nmax1)), A)
    result_len0 = result.len0
    result_len1 = result.len1
    Cn = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, ))
    kaves = np.ctypeslib.as_array(result.kavearr, shape=(result_len0, result_len1, ))
    deltas = np.ctypeslib.as_array(result.darr, shape=(result_len0, result_len1, ))

    return Cn, kaves, deltas




def run2d_f90():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[:, 0, 0, :]*1.0
    condition = np.ones(data.shape, dtype=bool)
    n = np.sum(data)*1.0
    print "n=", n
    maxxwidth = np.min(np.sum(condition, axis=0)) / 2
    maxywidth = np.min(np.sum(condition, axis=1)) / 2
    maxw = np.array([maxxwidth, maxywidth])
    print maxw
    A = np.cumsum(np.cumsum(data, axis=0), axis=1)

    Cn, kaves, deltas = calc_cost2d_f90(A, maxw, condition)
    print "kaves",kaves[0,0:], kaves.shape
    #print Cn

    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)
    opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1 )
    print "opt_indx", opt_indx

    m = 0.1*n

    ex = (1/m - 1/n) * kaves / (deltas**2*n) 

    Cm = ex + Cn

    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1 )
    print "opt_indx", opt_indx



run2d_f90()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


