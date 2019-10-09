#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./costfort1d.so")



def calc_cost1d_f90(A, maxw, condition):
    class result(ctypes.Structure):
        _fields_ =[("len0", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double)), ("kavearr", ctypes.POINTER(ctypes.c_double)), ("darr", ctypes.POINTER(ctypes.c_double))]

    lib.cost1d.restype = result
    lib.cost1d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]

    Nmax = A.shape[0]

    result = lib.cost1d(ctypes.byref(ctypes.c_int(maxw)), ctypes.byref(ctypes.c_int(Nmax)), A)
    result_len0 = result.len0
    Cn = np.ctypeslib.as_array(result.arr, shape=(result_len0, ))
    kaves = np.ctypeslib.as_array(result.kavearr, shape=(result_len0, ))
    deltas = np.ctypeslib.as_array(result.darr, shape=(result_len0, ))

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
    A = np.cumsum(data)

    Cn, kaves, deltas = calc_cost1d_f90(A, maxw, condition)
    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)
    print "Cn", Cn
    print "argmin(Cn)",np.argmin(Cn)

    m = 1.0*n

    ex = (1/m - 1/n) * kaves / (deltas**2*n) 
    ex[0] = 0.0

    Cm = ex + Cn
    print "Cm", Cm[:]
    print maxw
    print "Cm.size", Cm.shape

    print "argmin(Cm)",np.argmin(Cm)

    #!if fflag == 1:
    #    k, klen, kcond = calc_hist1d_f90(cumdata, 3, condition) 
    #else:
    #    k, cond = calc_hist1d(cumdata, 3, condition)

    #plt.figure(figsize=(16, 8)) 
    #plt.plot(k)
    plt.figure(figsize=(16, 8)) 
    #plt.plot(data)
    plt.plot(Cm)

run1d()
#k, klen = run_tst1d()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


