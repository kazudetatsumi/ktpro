#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./costfort3d.so")

def calc_cost3d_f90(A, maxw, condition):
    class result(ctypes.Structure):
        _fields_ =[("len0", ctypes.c_int), ("len1", ctypes.c_int), ("len2", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double)), ("kavearr", ctypes.POINTER(ctypes.c_double)), ("darr", ctypes.POINTER(ctypes.c_double))]

    lib.cost3d.restype = result
    lib.cost3d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), 
                           ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=3)]

    Nmax0 = A.shape[0]
    Nmax1 = A.shape[1]
    Nmax2 = A.shape[2]

    result = lib.cost3d(ctypes.byref(ctypes.c_int(maxw[0])), ctypes.byref(ctypes.c_int(maxw[1])), ctypes.byref(ctypes.c_int(maxw[1])), 
                        ctypes.byref(ctypes.c_int(Nmax0)), ctypes.byref(ctypes.c_int(Nmax1)), ctypes.byref(ctypes.c_int(Nmax2)), A)
    result_len0 = result.len0
    result_len1 = result.len1
    result_len2 = result.len2
    Cn = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, result_len2, ))
    kaves = np.ctypeslib.as_array(result.kavearr, shape=(result_len0, result_len1, result_len2, ))
    deltas = np.ctypeslib.as_array(result.darr, shape=(result_len0, result_len1, result_len2, ))

    return Cn, kaves, deltas





def run_simu3d():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    #datafile = "data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:]*1.0 # nqx, nqy, nqz, nomega
    data = np.sum(data[:, :, 0:5, :],axis=2)*1.0
    #data = data[:, :, 0, :]*1.0
    fflag = 1
    condition = np.ones(data.shape, dtype=bool)
    
    n = np.sum(data)*1.0
    print "n=", n


    maxxwidth = np.min(np.sum(condition, axis=0)) / 12
    maxywidth = np.min(np.sum(condition, axis=1)) / 12
    maxzwidth = np.min(np.sum(condition, axis=2)) / 12
    
    maxw = np.array([maxxwidth, maxywidth, maxzwidth])
    A = np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2)
    Cn, kaves, deltas = calc_cost3d_f90(A, maxw, condition)
    opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1 )
    print "opt_indx for Cn", opt_indx

    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

    m = 0.1*n

    ex = (1/m - 1/n) * kaves / (deltas**2*n) 

    Cm = ex + Cn

    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1 )
    print "opt_indx for Cm with m/n=", m/n, ":", opt_indx





run_simu3d()
#run_tst3d()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


