#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./histfort1d.so")


def calc_hist1d(A, nw, condition):
    Nmax = A.shape[0]
    N = int((Nmax - (Nmax % nw)) / nw)
    k = np.zeros((N))
    kcond = np.zeros((N))
    for i in range(0, N):
        ihead = (i+1)*nw - 1
        if i == 0:
            k[i] = A[ihead]
        else:
            k[i] = A[ihead] - A[ihead - nw]

    for i in range(0, N):
        kcond[i] = np.sum(condition[i*nw:(i+1)*nw])
    return k, kcond


def calc_hist1d_f90(A, nw0, condition):
    class result(ctypes.Structure):
        _fields_ =[("len0", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double))]

    lib.hist1d.restype = result
    lib.hist1d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]

    #lib.delete_array.restype = None
    #lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]

    Nmax = A.shape[0]
    print "Nmax", Nmax
    print "nw0", nw0
    result = lib.hist1d(ctypes.byref(ctypes.c_int(nw0)), ctypes.byref(ctypes.c_int(Nmax)), A)
    result_len0 = result.len0
    result_vec = np.ctypeslib.as_array(result.arr, shape=(result_len0, ))
    print result_vec

    #lib.delete_array(ctypes.byref(ctypes.c_int(result_len0)), result_vec)

    kcond = np.zeros((result_len0))
    for i in range(0, result_len0):
        kcond[i] = np.sum(condition[i*nw0:(i+1)*nw0])
    return result_vec, result_len0, kcond


def calc_cost1d(cumdata, maxw, condition, fflag):
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))
    for i in range(1, maxw):
       if fflag == 1:
          k, klen, kcond = calc_hist1d_f90(cumdata, i, condition) 
       else:
          k,  kcond = calc_hist1d(cumdata, i, condition) 
       kave = np.average(k)
       v = np.var(k)
       cost = (2 * kave - v) / (i**2)
       Cn[i] = cost
       kaves[i] = kave
       deltas[i] = (i*1.0)

       if fflag == 1:
          lib.delete_array.restype = None
          lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
          lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)

    deltas[0] = 1.0
    return Cn, kaves, deltas


def make_mappable(maxvalue):
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=maxvalue)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("jet")
    mappable=ScalarMappable(norm=norm,cmap=cmap)
    mappable._A = []
    return mappable



def run_tst1d():
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[0, 0, 0, :]*1.0
    condition = np.ones(data.shape, dtype=bool)
    
    A = np.cumsum(data)

    k, klen, kcond = calc_hist1d_f90(A, 8, condition) 
    k, klen, kcond = calc_hist1d_f90(A, 4, condition) 
    #lib = ctypes.CDLL("./histfort.so")
    #print "hist finished"
    #for i in range(0,100):
    #    print i,np.average(k[:,i,:])

    #plt.figure(figsize=(8, 16))
    #plt.plot(data)
    #plt.figure(figsize=(8, 16))
    #plt.plot(k)
    return k, klen


def run1d():
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    f = h5py.File(datafile)
    fflag = 1
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[0, 0, 0, :]*1.0
    n = np.sum(data)*1.0
    print "n=", n
    condition = np.ones(data.shape, dtype=bool)
    maxw =  int(data.shape[0] / 2)
    cumdata = np.cumsum(data)

    Cn, kaves, deltas = calc_cost1d(cumdata, maxw, condition, fflag)
    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)
    print "argmin(Cn)",np.argmin(Cn)

    m = 0.10*n

    ex = (1/m - 1/n) * kaves / (deltas**2*n) 
    ex[0] = 0.0

    Cm = ex + Cn

    print "argmin(Cm)",np.argmin(Cm)

    if fflag == 1:
        k, klen, kcond = calc_hist1d_f90(cumdata, 3, condition) 
    else:
        k, cond = calc_hist1d(cumdata, 3, condition)

    plt.figure(figsize=(16, 8)) 
    plt.plot(k)
    plt.figure(figsize=(16, 8)) 
    #plt.plot(data)
    plt.plot(Cm)

    if fflag == 1:
       lib.delete_array.restype = None
       lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
       lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)


run1d()
#k, klen = run_tst1d()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


