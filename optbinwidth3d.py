#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./histfort3d.so")



def calc_hist3d_f90(A, nw0, nw1, nw2, condition):

    class result(ctypes.Structure):
        _fields_ = [("len0", ctypes.c_int), ("len1", ctypes.c_int), ("len2", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double))]

    lib.hist3d.restype = result
    lib.hist3d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                           ctypes.POINTER(ctypes.c_int),  ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=3)]

    Nmax0 = A.shape[0]
    Nmax1 = A.shape[1]
    Nmax2 = A.shape[2]

    result = lib.hist3d(ctypes.byref(ctypes.c_int(nw0)), ctypes.byref(ctypes.c_int(nw1)), ctypes.byref(ctypes.c_int(nw2)),
                        ctypes.byref(ctypes.c_int(Nmax0)), ctypes.byref(ctypes.c_int(Nmax1)), ctypes.byref(ctypes.c_int(Nmax2)), A)
    result_len0 = result.len0
    result_len1 = result.len1
    result_len2 = result.len2
    result_vec = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, result_len2,))
    #print result_vec

    kcond = np.zeros((result_len0, result_len1, result_len2))
    for i in range(0, result_len0):
        for j in range(0, result_len1):
            for h in range(0, result_len2):
                kcond[i, j, h] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1, h*nw2:(h+1)*nw2])
    return result_vec, result_len0, result_len1, result_len2, kcond


def calc_hist3d(A, nw0, nw1, nw2, condition):
    Nmax = A.shape
    N0 = int((Nmax[0] - (Nmax[0] % nw0)) / nw0)
    N1 = int((Nmax[1] - (Nmax[1] % nw1)) / nw1)
    N2 = int((Nmax[2] - (Nmax[2] % nw2)) / nw2)
    k = np.zeros((N0, N1, N2))
    kcond = np.zeros((N0, N1, N2))
    for i in range(0, N0):
        ihead = (i+1)*nw0 - 1
        for j in range(0, N1):
            jhead = (j+1)*nw1 - 1
            for h in range(0, N2):
                hhead = (h+1)*nw2 - 1
                if i == 0 and j == 0 and h == 0:
                    k[i, j, h] = A[ihead, jhead, hhead]

                elif i != 0 and j == 0 and h == 0:                                            # procedure like 1D
                    k[i, j, h] = A[ihead, jhead, hhead] - A[ihead - nw0, jhead, hhead]
                elif i == 0 and j != 0 and h == 0:                                            # procedure like 1D
                    k[i, j, h] = A[ihead, jhead, hhead] - A[ihead, jhead - nw1, hhead]
                elif i == 0 and j == 0 and h != 0:                                            # procedure like 1D
                    k[i, j, h] = A[ihead, jhead, hhead] - A[ihead, jhead, hhead - nw2]

                elif i != 0 and j != 0 and h == 0:                                                                                                         # procedure like 2D
                    k[i, j, h] = A[ihead, jhead, hhead] - A[ihead - nw0, jhead, hhead] - A[ihead, jhead - nw1, hhead] + A[ihead - nw0, jhead - nw1, hhead]
                elif i != 0 and j == 0 and h != 0:                                                                                                         # procedure like 2D
                    k[i, j, h] = A[ihead, jhead, hhead] - A[ihead - nw0, jhead, hhead] - A[ihead, jhead, hhead - nw2] + A[ihead - nw0, jhead, hhead - nw2]
                elif i == 0 and j != 0 and h != 0:                                                                                                         # procedure like 2D
                    k[i, j, h] = A[ihead, jhead, hhead] - A[ihead, jhead - nw1, hhead] - A[ihead, jhead, hhead - nw2] + A[ihead, jhead - nw1, hhead - nw2]

                else:
                    k[i, j, h] = A[ihead, jhead, hhead] \
                              - A[ihead - nw0, jhead, hhead] - A[ihead, jhead - nw1, hhead] - A[ihead, jhead, hhead - nw2] \
                              + A[ihead - nw0, jhead - nw1, hhead] + A[ihead, jhead - nw1, hhead - nw2] + A[ihead - nw0, jhead, hhead - nw2] \
                              - A[ihead - nw0, jhead - nw1, hhead - nw2]
                
    for i in range(0, N0):
        for j in range(0, N1):
            for h in range(0, N2):
                kcond[i, j, h] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1, h*nw2:(h+1)*nw2])
    return k, kcond



def calc_cost3d(A, data, maxw, condition, fflag):
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))
    print maxw
    for i in range(1, maxw[0]):
        for j in range(1, maxw[1]):
            for h in range(1, maxw[2]):
                if fflag == 1:
                    k, klen0, klen1, klen2, kcond = calc_hist3d_f90(A, i, j, h, condition)
                else:
                    k, kcond = calc_hist3d(A, i, j, h, condition)

                knonzero = np.extract(np.max(kcond) == kcond, k)
                if i == 1 and j == 1 and h == 1:
                    print "shape of k matrix with zero elements",k.shape
                    print "total number of k is", k.shape[0]*k.shape[1]*k.shape[2]
                    print "number of nonzero k is",knonzero.shape 

                kave = np.average(knonzero)
                v = np.var(knonzero)
                cost = (2 * kave - v) / ((i*j*h)**2*1.0)
                if fflag == 1:
                    lib.delete_array3.restype = None
                    lib.delete_array3.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                                      ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=3)]
                    lib.delete_array3(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)),
                                         ctypes.byref(ctypes.c_int(klen2)), k)
                Cn[i, j, h] = cost
                kaves[i, j, h] = kave
                deltas[i, j, h] = (i*j*h*1.0)
                print "cost with (i, j, h) = ", i, j, h, ":", cost, "kave", kave, "v", v
    return Cn, kaves, deltas


def make_mappable(maxvalue):
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=maxvalue)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("jet")
    mappable=ScalarMappable(norm=norm,cmap=cmap)
    mappable._A = []
    return mappable



def run_simu3d():
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
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
    Cn, kaves, deltas = calc_cost3d(A, data, maxw, condition, fflag)

    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

    m = 1.0*n

    ex = (1/m - 1/n) * kaves / (deltas**2*n) 
    ex[0, :, :] = 0.0
    ex[:, 0, :] = 0.0
    ex[:, :, 0] = 0.0

    Cm = ex + Cn

    print "opt bin index", np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    if fflag == 1:
        k, klen0, klen1, klen2, kcond = calc_hist3d_f90(A, opt_indx[0], opt_indx[1], opt_indx[2], condition) 
    else:
        k, kcond = calc_hist3d(A, opt_indx[0], opt_indx[1], opt_indx[2], condition) 

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k[:,:,0]), vmax=np.max(k[:,:,0]), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data[:, :, 0]), vmax=np.max(data[:, :, 0]), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)

    if fflag == 1:
        lib.delete_array3.restype = None
        lib.delete_array3.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=3)]
        lib.delete_array3(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)),
                          ctypes.byref(ctypes.c_int(klen2)), k)



def run_tst3d():
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[:, :, 0, :]*1.0
    fflag = 1
    condition = np.ones(data.shape, dtype=bool)
    
    A = np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2)

    if fflag == 1:
        k, klen0, klen1, klen2,  kcond = calc_hist3d_f90(A, 4, 4, 4, condition) 
    else:
        k, kcond = calc_hist3d(A, 2, 2, 2, condition) 

    #for i in range(0,100):
    #    print i,np.average(k[:,i,:])
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k[:, :, 1]), vmax=np.max(k[:, :, 1]), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data[:, :, 1]), vmax=np.max(data[:, :, 1]), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)

    if fflag == 1:
        lib.delete_array3.restype = None
        lib.delete_array3.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=3)]
        lib.delete_array3(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), ctypes.byref(ctypes.c_int(klen2)), k)


run_simu3d()
#run_tst3d()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


