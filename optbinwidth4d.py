#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./histfort4d.so")

def calc_hist4d_f90(A, nw0, nw1, nw2, nw3,  condition):

    class result(ctypes.Structure):
        _fields_ = [("len0", ctypes.c_int), ("len1", ctypes.c_int), ("len2", ctypes.c_int),  ("len3", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double))]

    lib.hist4d.restype = result
    lib.hist4d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), 
                           ctypes.POINTER(ctypes.c_int),  ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]

    Nmax0 = A.shape[0]
    Nmax1 = A.shape[1]
    Nmax2 = A.shape[2]
    Nmax3 = A.shape[3]

    result = lib.hist4d(ctypes.byref(ctypes.c_int(nw0)), ctypes.byref(ctypes.c_int(nw1)), ctypes.byref(ctypes.c_int(nw2)), ctypes.byref(ctypes.c_int(nw3)),
                        ctypes.byref(ctypes.c_int(Nmax0)), ctypes.byref(ctypes.c_int(Nmax1)), ctypes.byref(ctypes.c_int(Nmax2)), ctypes.byref(ctypes.c_int(Nmax3)), A)
    result_len0 = result.len0
    result_len1 = result.len1
    result_len2 = result.len2
    result_len3 = result.len3
    result_vec = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, result_len2, result_len3, ))
    #print result_vec

    kcond = np.zeros((result_len0, result_len1, result_len2, result_len3))
    for i in range(0, result_len0):
        for j in range(0, result_len1):
            for h in range(0, result_len2):
                for l in range(0, result_len3):
                    kcond[i, j, h, l] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1, h*nw2:(h+1)*nw2, l*nw3:(l+1)*nw3])
    return result_vec, result_len0, result_len1, result_len2, result_len3, kcond


def calc_hist4d(A, data, nw0, nw1, nw2, nw3, condition):
    Nmax = A.shape
    N0 = int((Nmax[0] - (Nmax[0] % nw0)) / nw0)
    N1 = int((Nmax[1] - (Nmax[1] % nw1)) / nw1)
    N2 = int((Nmax[2] - (Nmax[2] % nw2)) / nw2)
    N3 = int((Nmax[3] - (Nmax[3] % nw3)) / nw3)
    k = np.zeros((N0, N1, N2, N3))
    kcond = np.zeros((N0, N1, N2, N3))
    
    ## make process faster and more complicated from here
    if nw0 == 1 and nw1 == 1 and nw2 == 1 and nw3 == 1:
        k = data
        kcond = condition
    elif nw0 != 1 and nw1 == 1 and nw2 == 1 and nw3 == 1:
        B = np.cumsum(data, axis=0)
        for j in range(0, N1):
            for h in range(0, N2):
                for l in range(0, N3):
                    k[:, j, h, l], kcond[:, j, h, l] = calc_hist1d(B[:, j, h, l], nw0, condition[:, j, h, l])
    elif nw0 == 1 and nw1 != 1 and nw2 == 1 and nw3 == 1:
        B = np.cumsum(data, axis=1)
        for i in range(0, N0):
            for h in range(0, N2):
                for l in range(0, N3):
                    k[i, :, h, l], kcond[i, :, h, l] = calc_hist1d(B[i, :, h, l], nw1, condition[i, :, h, l])
    elif nw0 == 1 and nw1 == 1 and nw2 != 1 and nw3 == 1:
        B = np.cumsum(data, axis=2)
        for i in range(0, N0):
            for j in range(0, N1):
                for l in range(0, N3):
                    k[i, j, :, l], kcond[i, j, :, l] = calc_hist1d(B[i, j, :, l], nw2, condition[i, j, :, l])
    elif nw0 == 1 and nw1 == 1 and nw2 == 1 and nw3 != 1:
        B = np.cumsum(data, axis=3)
        for i in range(0, N0):
            for j in range(0, N1):
                for h in range(0, N2):
                    k[i, j, h, :], kcond[i, j, h, :] = calc_hist1d(B[i, j, h, :], nw3, condition[i, j, h, :])
    ## make process faster and more complicated up to here
    else:
        for i in range(0, N0):
            for j in range(0, N1):
                for h in range(0, N2):
                    for l in range(0, N3):
                        kcond[i, j, h, l] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1, h*nw2:(h+1)*nw2, l*nw3:(l+1)*nw3])
        for i in range(0, N0):
            ihead = (i+1)*nw0 - 1
            for j in range(0, N1):
                jhead = (j+1)*nw1 - 1
                for h in range(0, N2):
                    hhead = (h+1)*nw2 - 1
                    for l in range(0, N3):
                        lhead = (l+1)*nw3 - 1
                        if i == 0 and j == 0 and h == 0 and l==0 :
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead]

                        elif i != 0 and j == 0 and h == 0 and l == 0:                                            # procedure like 1D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead - nw0, jhead, hhead, lhead]
                        elif i == 0 and j != 0 and h == 0 and l == 0:                                            # procedure like 1D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead, jhead - nw1, hhead, lhead]
                        elif i == 0 and j == 0 and h != 0 and l == 0:                                            # procedure like 1D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead, jhead, hhead - nw2, lhead]
                        elif i == 0 and j == 0 and h == 0 and l != 0:                                            # procedure like 1D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead, jhead, hhead, lhead - nw3]

                        elif i != 0 and j != 0 and h == 0 and l == 0:                                                                                                         # procedure like 2D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead - nw0, jhead, hhead, lhead] - A[ihead, jhead - nw1, hhead, lhead] + A[ihead - nw0, jhead - nw1, hhead, lhead]
                        elif i != 0 and j == 0 and h != 0 and l == 0:                                                                                                         # procedure like 2D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead - nw0, jhead, hhead, lhead] - A[ihead, jhead, hhead - nw2, lhead] + A[ihead - nw0, jhead, hhead - nw2, lhead]
                        elif i == 0 and j != 0 and h != 0 and l == 0:                                                                                                         # procedure like 2D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead, jhead - nw1, hhead, lhead] - A[ihead, jhead, hhead - nw2, lhead] + A[ihead, jhead - nw1, hhead - nw2, lhead]
                        elif i != 0 and j == 0 and h == 0 and l != 0:                                                                                                         # procedure like 2D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead - nw0, jhead, hhead, lhead] - A[ihead, jhead, hhead, lhead - nw3] + A[ihead - nw0, jhead, hhead, lhead - nw3]
                        elif i == 0 and j != 0 and h == 0 and l != 0:                                                                                                         # procedure like 2D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead, jhead - nw1, hhead, lhead] - A[ihead, jhead, hhead, lhead - nw3] + A[ihead, jhead - nw1, hhead, lhead - nw3]
                        elif i == 0 and j == 0 and h != 0 and l != 0:                                                                                                         # procedure like 2D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] - A[ihead, jhead, hhead - nw2, lhead] - A[ihead, jhead, hhead, lhead - nw3] + A[ihead, jhead, hhead - nw2, lhead - nw3]

                        elif i != 0 and j != 0 and h != 0 and l == 0:                                                                                                         # procedure like 3D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] \
                                      - A[ihead - nw0, jhead, hhead, lhead] - A[ihead, jhead - nw1, hhead, lhead] - A[ihead, jhead, hhead - nw2, lhead] \
                                      + A[ihead - nw0, jhead - nw1, hhead, lhead] + A[ihead, jhead - nw1, hhead - nw2, lhead] + A[ihead - nw0, jhead, hhead - nw2, lhead] \
                                      - A[ihead - nw0, jhead - nw1, hhead - nw2, lhead]
                        elif i != 0 and j != 0 and h == 0 and l != 0:                                                                                                         # procedure like 3D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] \
                                      - A[ihead - nw0, jhead, hhead, lhead] - A[ihead, jhead - nw1, hhead, lhead] - A[ihead, jhead, hhead, lhead - nw3] \
                                      + A[ihead - nw0, jhead - nw1, hhead, lhead] + A[ihead, jhead - nw1, hhead, lhead - nw3] + A[ihead - nw0, jhead, hhead, lhead - nw3] \
                                      - A[ihead - nw0, jhead - nw1, hhead, lhead - nw3]
                        elif i != 0 and j == 0 and h != 0 and l != 0:                                                                                                         # procedure like 3D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] \
                                      - A[ihead - nw0, jhead, hhead, lhead] - A[ihead, jhead, hhead - nw2, lhead] - A[ihead, jhead, hhead, lhead - nw3] \
                                      + A[ihead - nw0, jhead, hhead - nw2, lhead] + A[ihead, jhead, hhead - nw2, lhead - nw3] + A[ihead - nw0, jhead, hhead, lhead - nw3] \
                                      - A[ihead - nw0, jhead, hhead - nw2, lhead - nw3]
                        elif i == 0 and j != 0 and h != 0 and l != 0:                                                                                                         # procedure like 3D
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] \
                                      - A[ihead, jhead - nw1, hhead, lhead] - A[ihead, jhead, hhead - nw2, lhead] - A[ihead, jhead, hhead, lhead - nw3] \
                                      + A[ihead, jhead - nw1, hhead - nw2, lhead] + A[ihead, jhead, hhead - nw2, lhead - nw3] + A[ihead, jhead - nw1, hhead, lhead - nw3] \
                                      - A[ihead, jhead - nw1, hhead - nw2, lhead - nw3]

                        else:
                            k[i, j, h, l] = A[ihead, jhead, hhead, lhead] \
                                      - A[ihead - nw0, jhead, hhead, lhead] - A[ihead, jhead - nw1, hhead, lhead] - A[ihead, jhead, hhead - nw2, lhead] - A[ihead, jhead, hhead, lhead - nw3] \
                                      + A[ihead - nw0, jhead - nw1, hhead, lhead] + A[ihead - nw0, jhead, hhead - nw2, lhead] + A[ihead - nw0, jhead, hhead, lhead - nw3] \
                                      + A[ihead, jhead - nw1, hhead - nw2, lhead] + A[ihead, jhead - nw1, hhead, lhead - nw3] + A[ihead, jhead, hhead - nw2, lhead - nw3] \
                                      - A[ihead, jhead - nw1, hhead - nw2, lhead - nw3] - A[ihead - nw0, jhead, hhead - nw2, lhead - nw3] - A[ihead - nw0, jhead - nw1, hhead, lhead - nw3] - A[ihead - nw0, jhead - nw1, hhead - nw2, lhead] \
                                      + A[ihead - nw0, jhead - nw1, hhead - nw2, lhead - nw3]


    return k, kcond

def calc_cost4d(A, data, maxw, condition, fflag):
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))
    print maxw
    for i in range(1, maxw[0]):
        for j in range(1, maxw[1]):
            for h in range(1, maxw[2]):
                for l in range(1, maxw[3]):
                    if fflag == 1:
                        k, klen0, klen1, klen2, klen3, kcond = calc_hist4d_f90(A, i, j, h, l, condition)
                    else:
                        k, kcond = calc_hist4d(A, data, i, j, h, l, condition)
                    #strict condition for nonzero,  probably better.
                    knonzero = np.extract(np.max(kcond) == kcond, k)
                    # soft condition for nonzero
                    #knonzero = np.extract(kcond, k)
                    if i == 1 and j == 1 and h == 1 and l == 1:
                        print "shape of k matrix with zero elements",k.shape
                        print "total number of k is", k.shape[0]*k.shape[1]*k.shape[2]*k.shape[3]
                        print "number of nonzero k is",knonzero.shape 
                    kave = np.average(knonzero)
                    v = np.var(knonzero)
                    cost = (2 * kave - v) / ((i*j*h*l)**2*1.0)
                    if fflag == 1:
                        lib.delete_array4.restype = None
                        lib.delete_array4.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                                          ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]
                        lib.delete_array4(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), ctypes.byref(ctypes.c_int(klen2)),
                                             ctypes.byref(ctypes.c_int(klen3)), k)
                    print "cost with (i, j, h, l) = ", i, j, h, l, ":", cost, "kave", kave, "v", v
                    Cn[i, j, h, l] = cost
                    kaves[i, j, h, l] = kave
                    deltas[i, j, h, l] = (i*j*h*l*1.0)
    return Cn, kaves, deltas


def make_mappable(maxvalue):
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=maxvalue)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("jet")
    mappable=ScalarMappable(norm=norm,cmap=cmap)
    mappable._A = []
    return mappable




def run_simu4d():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    #datafile = "data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:]*1.0 # nqx, nqy, nqz, nomega
    #print "size of data is", data.shape
    #data = np.sum(data[:, :, :, :],axis=2)
    condition = np.ones(data.shape, dtype=bool)
    fflag = 1
    
    n = np.sum(data)*1.0
    print "n=", n


    maxxwidth = np.min(np.sum(condition, axis=0)) / 4
    maxywidth = np.min(np.sum(condition, axis=1)) / 4
    maxzwidth = np.min(np.sum(condition, axis=2)) / 4
    maxowidth = np.min(np.sum(condition, axis=3)) / 4
    maxxwidth = 2
    maxywidth = 2
    maxzwidth = 2
    maxowidth = 3
    
    maxw = np.array([maxxwidth, maxywidth, maxzwidth, maxowidth])
    A = np.cumsum(np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2), axis=3)
    Cn, kaves, delstas = calc_cost4d(A, data, maxw, condition, fflag)
    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

    m = 1.0*n

    ex = (1/(m*1.0) - 1/(n*1.0)) * kaves / (delstas**2*n) 
    ex[0, :, :, :] = 0.0
    ex[:, 0, :, :] = 0.0
    ex[:, :, 0, :] = 0.0
    ex[:, :, :, 0] = 0.0

    #Cm = ex + Cn
    Cm =  Cn
    Cm[0:2, 0:2, 0:2, 0:2] = 0.0

    print "opt bin index", np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)

    if fflag == 1:
        k, klen0, klen1, klen2, klen3, kcond = calc_hist4d_f90(A, data, opt_indx[0], opt_indx[1], opt_indx[2], opt_indx[3], condition) 
    else:
        k, kcond = calc_hist4d(A, data, opt_indx[0], opt_indx[1], opt_indx[2], opt_indx[3], condition) 

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k[:, 0, 0, :]), vmax=np.max(k[:, 0, 0, :]), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data[:, 0, 0, :]), vmax=np.max(data[:, 0, 0, :]), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)

    if fflag == 1:
        lib.delete_array4.restype = None
        lib.delete_array4.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]
        lib.delete_array4(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), ctypes.byref(ctypes.c_int(klen2)),
                          ctypes.byref(ctypes.c_int(klen3)), k)


def run_tst4d():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:]*1.0 # nqx, nqy, nqz, nomega
    condition = np.ones(data.shape, dtype=bool)
    fflag = 1
    
    A = np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2)
    print "cumsum finished"

    if fflag == 1:
        k, klen0, klen1, klen2, klen3, kcond = calc_hist4d_f90(A, 2, 2, 2, 2, condition)
    else:
        k, kcond = calc_hist4d(A, data, 2, 2, 2, 2, condition)
    print "hist finished"
    #for i in range(0,100):
    #    print i,np.average(k[:,i,:])

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k[:, :, 2, 2]), vmax=np.max(k[:, :, 2, 2]), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data[:, :, 2, 2]), vmax=np.max(data[:, :, 2, 2]), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)

    if fflag == 1:
        lib.delete_array4.restype = None
        lib.delete_array4.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]
        lib.delete_array4(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), ctypes.byref(ctypes.c_int(klen2)),
                          ctypes.byref(ctypes.c_int(klen3)), k)


#run_tst4d()
run_simu4d()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


