#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("./histfort_.so")


def get2ddata(f, xi, xf, yi, yf):
    data = np.genfromtxt(f,  delimiter=',', dtype=None)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    print "max intensity", np.max(z)
    dx = 0.005
    dy = 0.1
    xlin = np.arange(min(x), max(x)+dx, dx)
    nx = xlin.shape[0]
    ylin = np.arange(min(y), max(y)+dy, dy)
    ny = ylin.shape[0]
    print nx, ny
    
    karr = np.zeros((nx, ny))
    karr2 = np.zeros((nx, ny))


    # here we check whether we correctly read the whole lines in the input file
    # and generate a matrix "condition" which describes whether the element location is included in the input file or not.

    for _x, _y, _z in zip(x, y, z):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        karr[xx, yy] = _z + 0.00000001

    condition = karr > 0.0000000001
    print condition.shape
    karrnonzero = np.extract(condition, karr)
    ndata = x.shape[0]
    if karrnonzero.shape[0] != x.shape[0]:
         print "num of nonzero karr is not num of data", karrnonzero.shape
    else:
         print "num of nonzero karr matches  num of data", karrnonzero.shape


    #here we regenerate the data matrix purely.

    for _x, _y, _z in zip(x, y, z):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        karr2[xx, yy] = _z

    return karr2[xi:xf, yi:yf], condition[xi:xf, yi:yf]
    #return data


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


def run2d():
    txtfile = "/home/kazu/data/20min_fine.txt"
    xi = 119
    xf = 218
    yi = 25
    yf = 340
    data, condition = get2ddata(txtfile, xi, xf, yi, yf)
    n = np.sum(data)
    print "n=", n
    maxw = np.array([int(data.shape[0] / 2), int(data.shape[1]) / 2])
    print maxw
    #print data[0:10]
    A = np.cumsum(np.cumsum(data, axis=0), axis=1)

    Cn, kaves, delstas = calc_cost2d(A, maxw, condition)
    print "opt bin index", np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    k = calc_hist2d(A, opt_indx[0], opt_indx[1]) 
    #k = calc_hist2d(A, 4, 4) 

    plt.figure(figsize=(16, 8))
    #plt.plot(np.sum(data, axis=0))
    #plt.pcolor(data, vmax=np.max(data)/10)
    plt.pcolor(np.transpose(k), vmax=np.max(k), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(16, 8))
    plt.pcolor(np.transpose(data), vmax=np.max(data)/500, cmap='jet')

    #from matplotlib.colors import Normalize
    #norm = Normalize(vmin=0, vmax=np.max(data)/100)
    #from matplotlib.cm import ScalarMappable, get_cmap
    #cmap = get_cmap("jet")
    #mappable=ScalarMappable(norm=norm,cmap=cmap)
    #mappable._A = []
    mappable = make_mappable(np.max(data)/500)
    plt.colorbar(mappable)


def runexorg():
    txtfile = "/home/kazu/data/20min_fine.txt"
    #xi = 119
    xi = 119
    #xf = 222
    xf = 218
    yi = 25
    yf = 340
    #yf = 530
    
    data, condition = get2ddata(txtfile, xi, xf, yi, yf)
    n = np.sum(data)*1.0
    print "n=", n


    maxxwidth = np.min(np.sum(condition, axis=0)) / 2
    maxywidth = np.min(np.sum(condition, axis=1)) / 2
    
    #maxw = np.array([int(data.shape[0] / 2), int(data.shape[1]) / 2])
    maxw = np.array([maxxwidth, maxywidth])
    print maxw
    A = np.cumsum(np.cumsum(data, axis=0), axis=1)

    Cn, kaves, delstas = calc_cost2d(A, maxw, condition)
    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

    m = 15.0*n

    ex = (1/m - 1/n) * kaves / (delstas**2*n) 
    ex[0, :] = 0.0
    ex[:, 0] = 0.0

    #print "ex, size",ex.shape
    #print "ex", ex[0:5, 0:5]
    #print Cn[0:5, 0:5]
    Cm = ex + Cn

    print "opt bin index", np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    k, kcond = calc_hist2d(A, opt_indx[0], opt_indx[1], condition) 
    #plt.figure(figsize=(16, 8))
    #plt.plot(Cn[10,:])
    #plt.plot(Cm[10,:])

    plt.figure(figsize=(16, 8))
    plt.pcolor(np.transpose(k), vmax=np.max(k), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(16, 8))
    plt.pcolor(np.transpose(data), vmax=np.max(data)/1000, cmap='jet')

    mappable = make_mappable(np.max(data)/1000)
    plt.colorbar(mappable)


def run_tst():
    txtfile = "/home/kazu/data/20min_fine.txt"
    txtfile = "/home/kazu/desktop/191031/all.txt"
    xi = 0
    xf = 630
    #xf = 727
    yi = 0
    yf = 441
    #yf = 471
    data, condition = get2ddata(txtfile, xi, xf, yi, yf)
    plt.figure(figsize=(16, 8))
    plt.pcolor(np.transpose(data), vmax=np.max(data/1000), cmap='jet')
    mappable = make_mappable(np.max(data)/1000)
    plt.colorbar(mappable)


def runex():
    txtfile = "/home/kazu/data/20min_fine.txt"
    txtfile = "/home/kazu/desktop/191031/all.txt"
    xi = 100
    #xi = 109
    #xf = 222
    #xf = 218
    xf = 212
    #yi = 25
    yi = 100
    #yf = 340
    yf = 255
    #yf = 530
    ### for 20min_fine.txt
    ###xi = 119
    ###xf = 218
    ###yi = 25 
    ###yf = 340

    data, condition = get2ddata(txtfile, xi, xf, yi, yf)
    n = np.sum(data)*1.0
    print "n=", n
    fflag = 0



    maxxwidth = np.min(np.sum(condition, axis=0)) / 2
    maxywidth = np.min(np.sum(condition, axis=1)) / 2
    
    #maxw = np.array([int(data.shape[0] / 2), int(data.shape[1]) / 2])
    maxw = np.array([maxxwidth, maxywidth])
    print maxw
    A = np.cumsum(np.cumsum(data, axis=0), axis=1)

    Cn, kaves, delstas = calc_cost2d_f90(A, maxw, condition, fflag)
    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

    m = 15.0*n

    ex = (1/m - 1/n) * kaves / (delstas**2*n) 
    ex[0, :] = 0.0
    ex[:, 0] = 0.0

    #print "ex, size",ex.shape
    #print "ex", ex[0:5, 0:5] #print Cn[0:5, 0:5]
    Cm = ex + Cn

    print "opt bin index", np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)

    if fflag == 1:
        k, klen0, klen1, kcond = calc_hist2d_f90(A, opt_indx[0], opt_indx[1], condition) 
    else:
        k, kcond = calc_hist2d(A, opt_indx[0], opt_indx[1], condition) 
    #plt.figure(figsize=(16, 8))
    #plt.plot(Cn[10,:])
    #plt.plot(Cm[10,:])

    plt.figure(figsize=(16, 8))
    #plt.pcolor(np.transpose(k), vmax=np.max(k), cmap='jet')
    plt.pcolor(np.transpose(k), vmax=np.max(k), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(16, 8))
    plt.pcolor(np.transpose(data), vmax=np.max(data), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)

    if fflag == 1:
        lib.delete_array2.restype = None
        lib.delete_array2.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)]
        lib.delete_array2(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), k)
#run_tst()
#run_tst2d()
#run2d()
#run2d_f90()
#runexorg()
#run_tst()
#run2d()
runex()
#run_simu()
#k, klen = run_tst1d()
#lib.delete_array.restype = None
#lib.delete_array.argtypes = [ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)]
#lib.delete_array(ctypes.byref(ctypes.c_int(klen)), k)
plt.show()


