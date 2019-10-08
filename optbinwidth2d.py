#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./histfort2d.so")


def calc_hist2d_f90(A, nw0, nw1,  condition):

    class result(ctypes.Structure):
        _fields_ = [("len0", ctypes.c_int), ("len1", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double))]

    lib.hist2d.restype = result
    lib.hist2d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                           ctypes.POINTER(ctypes.c_int),  ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)]

    Nmax0 = A.shape[0]
    Nmax1 = A.shape[1]
    print nw0, nw1, Nmax0, Nmax1

    result = lib.hist2d(ctypes.byref(ctypes.c_int(nw0)), ctypes.byref(ctypes.c_int(nw1)),
                        ctypes.byref(ctypes.c_int(Nmax0)), ctypes.byref(ctypes.c_int(Nmax1)), A)
    result_len0 = result.len0
    result_len1 = result.len1
    result_vec = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, ))
    print result_vec

    kcond = np.zeros((result_len0, result_len1))
    for i in range(0, result_len0):
        for j in range(0, result_len1):
            kcond[i, j] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1])
    return result_vec, result_len0, result_len1, kcond


def calc_hist2d(A, nw0, nw1, condition):
    Nmax = A.shape
    N0 = int((Nmax[0] - (Nmax[0] % nw0)) / nw0)
    N1 = int((Nmax[1] - (Nmax[1] % nw1)) / nw1)
    k = np.zeros((N0, N1))
    kcond = np.zeros((N0, N1))
    for i in range(0, N0):
        ihead = (i+1)*nw0 - 1
        for j in range(0, N1):
            jhead = (j+1)*nw1 - 1
            if i == 0 and j == 0:
                k[i, j] = A[ihead, jhead]
            elif j == 0 and i != 0:
                k[i, j] = A[ihead, jhead] - A[ihead - nw0, jhead]
            elif i == 0 and j != 0:
                k[i, j] = A[ihead, jhead] - A[ihead, jhead - nw1]
            else:
                k[i, j] = A[ihead, jhead] - A[ihead - nw0, jhead] - A[ihead, jhead - nw1] + A[ihead - nw0, jhead - nw1]
    for i in range(0, N0):
        for j in range(0, N1):
            kcond[i, j] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1])
    return k, kcond


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


def calc_cost2d(A, maxw, condition):
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))
    for i in range(1, maxw[0]):
       for j in range(1, maxw[1]):
          k, kcond = calc_hist2d(A, i, j, condition)
          #strict condition for nonzero,  probably better.
          knonzero = np.extract(np.max(kcond) == kcond, k)
          # soft condition for nonzero
          #knonzero = np.extract(kcond, k)
          if i == 1 and j ==1:
             print "shape of k matrix with zero elements",k.shape
             print "total number of k is", k.shape[0]*k.shape[1]
             print "number of nonzero k is",knonzero.shape 
          kave = np.average(knonzero)
          v = np.var(knonzero)
          
          # in order to effectively remove the artificial zero elements in original data, we divide the average with fracdata
          #kf = k.flatten()                            
          #kave = np.average(kf)/fracdata
          #v =  np.sum((kf - kave)**2)/(kf.shape[0]*1.0*fracdata)
          cost = (2 * kave - v) / ((i*j)**2*1.0)
          print "cost with (i, j) = ", i, j, ":", cost, "kave", kave, "v", v
          Cn[i, j] = cost
          kaves[i, j] = kave
          deltas[i, j] = (i*j*1.0)
    return Cn, kaves, deltas


def calc_cost2d_f90(A, maxw, condition, fflag):
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))
    for i in range(1, maxw[0]):
        for j in range(1, maxw[1]):
            if fflag == 1:
                k, klen0, klen1, kcond = calc_hist2d_f90(A, i, j, condition)
            else:
                k, kcond = calc_hist2d(A, i, j, condition)
            knonzero = np.extract(np.max(kcond) == kcond, k)
            if i == 1 and j == 1:
                print "shape of k matrix with zero elements", k.shape
                print "total number of k is", k.shape[0]*k.shape[1]
                print "number of nonzero k is", knonzero.shape
            kave = np.average(knonzero)
            v = np.var(knonzero)
          
            cost = (2 * kave - v) / ((i*j)**2*1.0)
            print "cost with (i, j) = ", i, j, ":", cost, "kave", kave, "v", v
            Cn[i, j] = cost
            kaves[i, j] = kave
            deltas[i, j] = (i*j*1.0)
            if fflag == 1:
                lib.delete_array2.restype = None
                lib.delete_array2.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)]
                lib.delete_array2(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), k)
    return Cn, kaves, deltas




def make_mappable(maxvalue):
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=maxvalue)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("jet")
    mappable=ScalarMappable(norm=norm,cmap=cmap)
    mappable._A = []
    return mappable


def run_tst2d():
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    f = h5py.File(datafile)
    fflag = 1
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[0, :, 0, :]*1.0
    condition = np.ones(data.shape, dtype=bool)
   
    A = np.cumsum(np.cumsum(data, axis=0), axis=1)

    if fflag == 1:
        k, klen0, klen1,  kcond = calc_hist2d_f90(A, 2, 2, condition)
    else:
        k, kcond = calc_hist2d(A, 2, 2, condition)

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k), vmax=np.max(k), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data), vmax=np.max(data), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)

    if fflag == 1:
        lib.delete_array2.restype = None
        lib.delete_array2.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)]
        lib.delete_array2(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), k)


def run2d_f90():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[:, 0, 0, :]*1.0
    fflag = 1
    condition = np.ones(data.shape, dtype=bool)
    n = np.sum(data)*1.0
    print "n=", n
    maxxwidth = np.min(np.sum(condition, axis=0)) / 2
    maxywidth = np.min(np.sum(condition, axis=1)) / 2
    maxw = np.array([maxxwidth, maxywidth])
    print maxw
    A = np.cumsum(np.cumsum(data, axis=0), axis=1)

    Cn, kaves, deltas = calc_cost2d_f90(A, maxw, condition, fflag)

    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)
    print "argmin(Cn)",np.argmin(Cn)

    m = 1.0*n

    ex = (1/m - 1/n) * kaves / (deltas**2*n) 
    ex[0, :] = 0.0
    ex[:, 0] = 0.0

    Cm = ex + Cn

    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    print "opt_indx", opt_indx

    if fflag == 1:
        k, klen0, klen1, kcond = calc_hist2d_f90(A, 9, 2, condition) 
    else:
        k, kcond = calc_hist2d(A, 9, 2, condition)

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k), vmax=np.max(k), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data), vmax=np.max(data), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)

    if fflag == 1:
        lib.delete_array2.restype = None
        lib.delete_array2.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)]
        lib.delete_array2(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), k)


#run_tst2d()
run2d_f90()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


