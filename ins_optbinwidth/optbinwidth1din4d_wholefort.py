#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib1d = ctypes.CDLL("./costfort1din4d.so")

def calc_cost1d_f90(B, maxw, CDB):
    class result1d(ctypes.Structure):
        _fields_ =[("len31d", ctypes.c_int),
                  ("arr1d", ctypes.POINTER(ctypes.c_double)), ("kavearr1d", ctypes.POINTER(ctypes.c_double)), ("darr1d", ctypes.POINTER(ctypes.c_double))]

    lib1d.cost1d.restype = result1d
    lib1d.cost1d.argtypes = [ctypes.POINTER(ctypes.c_int), 
                           ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.int32, ndim=4)]

    Nmax0 = B.shape[0]
    Nmax1 = B.shape[1]
    Nmax2 = B.shape[2]
    Nmax3 = B.shape[3]

    result1d = lib1d.cost1d(ctypes.byref(ctypes.c_int(maxw)),
                        ctypes.byref(ctypes.c_int(Nmax0)),
                        ctypes.byref(ctypes.c_int(Nmax1)),
                        ctypes.byref(ctypes.c_int(Nmax2)),
                        ctypes.byref(ctypes.c_int(Nmax3)), B, CDB)
    result1d_len3 = result1d.len31d
    Cn = np.ctypeslib.as_array(result1d.arr1d, shape=(result1d_len3, ))
    kaves = np.ctypeslib.as_array(result1d.kavearr1d, shape=(result1d_len3, ))
    deltas = np.ctypeslib.as_array(result1d.darr1d, shape=(result1d_len3, ))
    print Cn

    return Cn, kaves, deltas, result1d_len3


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
    
    n = np.sum(data)*1.0
    print "n=", n


    maxxwidth = np.min(np.sum(condition, axis=0)) / 9
    maxywidth = np.min(np.sum(condition, axis=1)) / 9
    maxzwidth = np.min(np.sum(condition, axis=2)) / 9
    maxowidth = np.min(np.sum(condition, axis=3)) / 9
    maxxwidth = 4
    maxywidth = 4
    maxzwidth = 4
    maxowidth = 4
    maxw = np.array([maxxwidth, maxywidth, maxzwidth, maxowidth])
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))


    # First, in the case of all of the widths are 1.
    deltas[0, 0, 0, 0] = 1
    knonzero = np.extract(np.max(condition) == condition, data)
    kave = np.average(knonzero)
    v = np.var(knonzero)
    cost = (2 * kave - v) / ((1.0)**2)
    kaves[0, 0, 0, 0] = kave
    Cn[0, 0, 0, 0] = cost
    print "Cn[0,0,0,0]:",cost


    # Next, in the case that one of the widths is not 1.

    B = np.cumsum(data, axis=0)
    CDB = np.cumsum(condition, axis=0, dtype='int32')
    Cntmp, kavestmp, deltastmp, result1d_len3  = calc_cost1d_f90(B, maxw[0], CDB)
    Cn[1:, 0, 0, 0] = Cntmp
    kaves[1:, 0, 0, 0] = kavestmp
    deltas[1:, 0, 0, 0] = deltastmp 
    del B, CDB

    lib1d.delete_array1d.restype = None
    lib1d.delete_array1d.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1) ]
    lib1d.delete_array1d(ctypes.byref(ctypes.c_int(result1d_len3)),  Cntmp, kavestmp, deltastmp)

    B = np.cumsum(data, axis=1)
    CDB = np.cumsum(condition, axis=1, dtype='int32')
    BB = B.transpose(1,0,2,3)
    CDBB = CDB.transpose(1,0,2,3)

    Cn[0, 1:, 0, 0], kaves[0, 1:, 0, 0], deltas[0, 1:, 0, 0], result1d_len3 = calc_cost1d_f90(BB, maxw[2], CDBB)
    B = np.transpose(np.cumsum(data, axis=2), (2, 1, 0, 3))
    CDB = np.transpose(np.cumsum(condition, axis=2, dtype='int32'), (2, 1, 0, 3))
    Cn[0, 0, 1:, 0], kaves[0, 0, 1:, 0], deltas[0, 0, 1:, 0], result1d_len3 = calc_cost1d_f90(B, maxw[2], CDB)
    B = np.transpose(np.cumsum(data, axis=3), (3, 1, 2, 0))
    CDB = np.transpose(np.cumsum(condition, axis=3, dtype='int32'), (3, 1, 2, 0))
    Cn[0, 0, 0, 1:], kaves[0, 0, 0, 1:], deltas[0, 0, 0, 1:], result1d_len3 = calc_cost1d_f90(B, maxw[3], CDB)



#run_tst4d()
run_simu4d()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
#plt.show()


