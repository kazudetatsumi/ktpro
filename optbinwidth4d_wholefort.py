#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./costfort4d.so")

#def calc_cost4d_f90(A, maxw, data, condition):
def calc_cost4d_f90(maxw, data, condition):
    class result(ctypes.Structure):
        _fields_ =[("len0", ctypes.c_int), ("len1", ctypes.c_int), ("len2", ctypes.c_int), ("len3", ctypes.c_int),
                  ("arr", ctypes.POINTER(ctypes.c_double)), ("kavearr", ctypes.POINTER(ctypes.c_double)), ("darr", ctypes.POINTER(ctypes.c_double))]

    lib.cost4d.restype = result
    lib.cost4d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                           ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                           #np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.int32, ndim=4)]

    Nmax0 = data.shape[0]
    Nmax1 = data.shape[1]
    Nmax2 = data.shape[2]
    Nmax3 = data.shape[3]

    #result = lib.cost4d(ctypes.byref(ctypes.c_int(maxw[0])), ctypes.byref(ctypes.c_int(maxw[1])), ctypes.byref(ctypes.c_int(maxw[2])), ctypes.byref(ctypes.c_int(maxw[3])),
    #                    ctypes.byref(ctypes.c_int(Nmax0)), ctypes.byref(ctypes.c_int(Nmax1)), ctypes.byref(ctypes.c_int(Nmax2)), ctypes.byref(ctypes.c_int(Nmax3)), A, data, condition)
    result = lib.cost4d(ctypes.byref(ctypes.c_int(maxw[0])), ctypes.byref(ctypes.c_int(maxw[1])), ctypes.byref(ctypes.c_int(maxw[2])), ctypes.byref(ctypes.c_int(maxw[3])),
                        ctypes.byref(ctypes.c_int(Nmax0)), ctypes.byref(ctypes.c_int(Nmax1)), ctypes.byref(ctypes.c_int(Nmax2)), ctypes.byref(ctypes.c_int(Nmax3)),  data, condition)
    result_len0 = result.len0
    result_len1 = result.len1
    result_len2 = result.len2
    result_len3 = result.len3
    Cn = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, result_len2, result_len3, ))
    Cn = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, result_len2, result_len3, ))
    kaves = np.ctypeslib.as_array(result.kavearr, shape=(result_len0, result_len1, result_len2, result_len3, ))
    deltas = np.ctypeslib.as_array(result.darr, shape=(result_len0, result_len1, result_len2, result_len3, ))

    return Cn, kaves, deltas


def run_simu4d():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_1000000.hdf5"
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    #datafile = "data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:]*1.0 # nqx, nqy, nqz, nomega
    #print "size of data is", data.shape
    #data = np.sum(data[:, :, :, :],axis=2)
    condition = np.ones(data.shape, dtype=np.int32)
    print condition.dtype
    
    n = np.sum(data)*1.0
    print "n=", n


    maxxwidth = np.min(np.sum(condition, axis=0)) / 9
    maxywidth = np.min(np.sum(condition, axis=1)) / 9
    maxzwidth = np.min(np.sum(condition, axis=2)) / 9
    maxowidth = np.min(np.sum(condition, axis=3)) / 9
    maxxwidth = 2
    maxywidth = 2
    maxzwidth = 2
    maxowidth = 1
    
    maxw = np.array([maxxwidth, maxywidth, maxzwidth, maxowidth])
    #!A = np.cumsum(np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2), axis=3)
    #Cn, kaves, deltas = calc_cost4d_f90(A, maxw, data, condition)
    Cn, kaves, deltas = calc_cost4d_f90(maxw, data, condition)
    opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1, opt_indx[3] + 1)
    print "opt_indx for Cn", opt_indx
    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

    m = 1.0*n

    ex = (1/(m*1.0) - 1/(n*1.0)) * kaves / (deltas**2*n) 

    Cm = ex + Cn

    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1, opt_indx[3] + 1)
    print "opt_indx for Cm with m/n=", m/n, ":", opt_indx


    len0 = Cn.shape[0]
    len1 = Cn.shape[1]
    len2 = Cn.shape[2]
    len3 = Cn.shape[3]

    lib.delete_array.restype = None
    lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=4), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]
    lib.delete_array(ctypes.byref(ctypes.c_int(len0)), ctypes.byref(ctypes.c_int(len1)), ctypes.byref(ctypes.c_int(len2)), ctypes.byref(ctypes.c_int(len3)), Cn, kaves, deltas)




#run_tst4d()
run_simu4d()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


