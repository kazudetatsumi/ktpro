#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
import sys
lib = ctypes.CDLL("./costfort4d_withcond.so")

def calc_cost4d_f90(maxw, data, condition, usecond, condparam):
    class result(ctypes.Structure):
        _fields_ =[("len0", ctypes.c_int), ("len1", ctypes.c_int), ("len2", ctypes.c_int), ("len3", ctypes.c_int),
                  ("arr", ctypes.POINTER(ctypes.c_double)), ("kavearr", ctypes.POINTER(ctypes.c_double)), ("darr", ctypes.POINTER(ctypes.c_double))]

    lib.cost4d.restype = result
    lib.cost4d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                           ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_bool),
                           ctypes.POINTER(ctypes.c_double),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.int32, ndim=4)]

    Nmax0 = data.shape[0]
    Nmax1 = data.shape[1]
    Nmax2 = data.shape[2]
    Nmax3 = data.shape[3]

    result = lib.cost4d(ctypes.byref(ctypes.c_int(maxw[0])), ctypes.byref(ctypes.c_int(maxw[1])), ctypes.byref(ctypes.c_int(maxw[2])), ctypes.byref(ctypes.c_int(maxw[3])),
                        ctypes.byref(ctypes.c_int(Nmax0)),
                        ctypes.byref(ctypes.c_int(Nmax1)),
                        ctypes.byref(ctypes.c_int(Nmax2)),
                        ctypes.byref(ctypes.c_int(Nmax3)),
                        ctypes.byref(ctypes.c_bool(usecond)),
                        ctypes.byref(ctypes.c_double(condparam)), data, condition)
    result_len0 = result.len0
    result_len1 = result.len1
    result_len2 = result.len2
    result_len3 = result.len3
    Cn = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, result_len2, result_len3, ))
    kaves = np.ctypeslib.as_array(result.kavearr, shape=(result_len0, result_len1, result_len2, result_len3, ))
    deltas = np.ctypeslib.as_array(result.darr, shape=(result_len0, result_len1, result_len2, result_len3, ))

    return Cn, kaves, deltas


def run_simu4d():
    datafile = "./eliminated_data.hdf5"
    f = h5py.File(datafile)
    data = f["data4"][:,:,:,:] # nqx, nqy, nqz, nomega
    print "size of data is", data.shape
    condition = np.array(f["condition"], dtype=np.int32)
    usecond = True
    condparam = float(sys.argv[1])  #  parameter for selection of elements taken into accounts of the statistics
    print "usecond:", usecond
    print "condparam:", condparam
    
    n = np.sum(data)*1.0
    print "n=", n

    maxxwidth = int(data.shape[0] // 2)
    maxywidth = int(data.shape[1] // 2)
    maxzwidth = int(data.shape[2] // 2)
    maxowidth = int(data.shape[3] // 2)
    print "maxwidth:", maxxwidth, maxywidth, maxzwidth, maxowidth

    
    maxw = np.array([maxxwidth, maxywidth, maxzwidth, maxowidth])
    Cn, kaves, deltas = calc_cost4d_f90(maxw, data, condition, usecond, condparam)
    opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1, opt_indx[3] + 1)
    print "opt_indx for Cn", opt_indx
    Cn2 = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

    m = 0.01*n

    ex = (1/(m*1.0) - 1/(n*1.0)) * kaves / (deltas**2*n) 

    Cm = ex + Cn2

    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1, opt_indx[3] + 1)
    print "opt_indx for Cm with m/n=", m/n, ":", opt_indx

    print "---save results in Cn.hdf5---"
    outfile = "Cn.hdf5"
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('Cn', data=Cn)
        hf.create_dataset('kave', data=kaves)
        hf.create_dataset('delta', data=deltas)

    len0 = Cn.shape[0]
    len1 = Cn.shape[1]
    len2 = Cn.shape[2]
    len3 = Cn.shape[3]

    print "---delete_array_pointer---"
    lib.delete_array_pointer.restype = None
    lib.delete_array_pointer.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                 np.ctypeslib.ndpointer(dtype=np.float64, ndim=4), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]
    lib.delete_array_pointer(ctypes.byref(ctypes.c_int(len0)), ctypes.byref(ctypes.c_int(len1)), ctypes.byref(ctypes.c_int(len2)), ctypes.byref(ctypes.c_int(len3)), Cn, kaves, deltas)


run_simu4d()



