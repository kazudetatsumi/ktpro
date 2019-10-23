#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py

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
       xx =  np.where(abs(xlin - _x) < 0.0000001)
       yy =  np.where(abs(ylin - _y) < 0.0000001)
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
       xx =  np.where(abs(xlin - _x) < 0.0000001)
       yy =  np.where(abs(ylin - _y) < 0.0000001)
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


def calc_hist4d(A, nw0, nw1, nw2, nw3, condition):
    Nmax = A.shape
    N0 = int((Nmax[0] - (Nmax[0] % nw0)) / nw0)
    N1 = int((Nmax[1] - (Nmax[1] % nw1)) / nw1)
    N2 = int((Nmax[2] - (Nmax[2] % nw2)) / nw2)
    N3 = int((Nmax[3] - (Nmax[3] % nw3)) / nw3)
    k = np.zeros((N0, N1, N2, N3))
    kcond = np.zeros((N0, N1, N2, N3))
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

    for i in range(0, N0):
        for j in range(0, N1):
            for h in range(0, N2):
                for l in range(0, N3):
                    kcond[i, j, h, l] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1, h*nw2:(h+1)*nw2, l*nw3:(l+1)*nw3])
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


def calc_cost3d(A, maxw, condition):
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))
    print maxw
    for i in range(1, maxw[0]):
        for j in range(1, maxw[1]):
            for h in range(1, maxw[2]):
                k, kcond = calc_hist3d(A, i, j, h, condition)
                #strict condition for nonzero,  probably better.
                knonzero = np.extract(np.max(kcond) == kcond, k)
                # soft condition for nonzero
                #knonzero = np.extract(kcond, k)
                if i == 1 and j == 1 and h == 1:
                    print "shape of k matrix with zero elements",k.shape
                    print "total number of k is", k.shape[0]*k.shape[1]*k.shape[2]
                    print "number of nonzero k is",knonzero.shape 
                kave = np.average(knonzero)
                v = np.var(knonzero)
          
                cost = (2 * kave - v) / ((i*j*h)**2*1.0)
                print "cost with (i, j, h) = ", i, j, h, ":", cost, "kave", kave, "v", v
                Cn[i, j, h] = cost
                kaves[i, j, h] = kave
                deltas[i, j, h] = (i*j*h*1.0)
    return Cn, kaves, deltas


def calc_cost4d(A, maxw, condition):
    Cn = np.zeros((maxw))
    kaves = np.zeros((maxw))
    deltas = np.zeros((maxw))
    print maxw
    for i in range(2, maxw[0]):
        for j in range(2, maxw[1]):
            for h in range(2, maxw[2]):
                for l in range(2, maxw[3]):
                    k, kcond = calc_hist4d(A, i, j, h, l, condition)
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


def run2d():
    txtfile = "/home/kazu/data/20min_fine.txt"
    xi = 119
    xf = 218
    yi = 25
    yf = 340
    data, condition = get2ddata(txtfile, xi, xf, yi, yf)
    n = np.sum(data)
    print "n=", n
    maxw =  np.array([int(data.shape[0] / 2), int(data.shape[1]) / 2])
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

def runex():
    txtfile = "/home/kazu/data/20min_fine.txt"
    #xi = 119
    xi = 119
    xf = 222
    yi = 25
    #yf = 340
    yf = 530
    
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
    data = get2ddata(txtfile)
    plt.pcolor(np.transpose(data), vmax=np.max(data)/1000, cmap='jet')

def run_simu():
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[:, 0, 0, :]
    condition = np.ones(data.shape, dtype=bool)
    
    n = np.sum(data)*1.0
    print "n=", n


    maxxwidth = np.min(np.sum(condition, axis=0)) / 2
    maxywidth = np.min(np.sum(condition, axis=1)) / 2
    
    maxw = np.array([maxxwidth, maxywidth])
    A = np.cumsum(np.cumsum(data, axis=0), axis=1)

    Cn, kaves, delstas = calc_cost2d(A, maxw, condition)
    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)

    m = 1*n

    ex = (1/m - 1/n) * kaves / (delstas**2*n) 
    ex[0, :] = 0.0
    ex[:, 0] = 0.0

    Cm = ex + Cn

    print "opt bin index", np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    k, kcond = calc_hist2d(A, opt_indx[0], opt_indx[1], condition) 

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k), vmax=np.max(k), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data), vmax=np.max(data), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)

def run_simu3d():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    datafile = "data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = np.sum(data[:, :, 0:5, :],axis=2)
    condition = np.ones(data.shape, dtype=bool)
    
    n = np.sum(data)*1.0
    print "n=", n


    maxxwidth = np.min(np.sum(condition, axis=0)) / 12
    maxywidth = np.min(np.sum(condition, axis=1)) / 12
    maxzwidth = np.min(np.sum(condition, axis=2)) / 12
    
    maxw = np.array([maxxwidth, maxywidth, maxzwidth])
    A = np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2)
    Cn, kaves, delstas = calc_cost3d(A, maxw, condition)

    m = 1.0*n

    ex = (1/m - 1/n) * kaves / (delstas**2*n) 
    ex[0, :, :] = 0.0
    ex[:, 0, :] = 0.0
    ex[:, :, 0] = 0.0

    Cm = ex + Cn

    print "opt bin index", np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    opt_indx = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    k, kcond = calc_hist3d(A, opt_indx[0], opt_indx[1], opt_indx[2], condition) 

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k[:,0,:]), vmax=np.max(k[:,0,:]), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data[:, 0, :]), vmax=np.max(data[:, 0, :]), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)


def run_simu4d():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    datafile = "data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    #data = np.sum(data[:, :, :, :],axis=2)
    condition = np.ones(data.shape, dtype=bool)
    
    n = np.sum(data)*1.0
    print "n=", n


    #maxxwidth = np.min(np.sum(condition, axis=0)) / 15
    #maxywidth = np.min(np.sum(condition, axis=1)) / 15
    #maxzwidth = np.min(np.sum(condition, axis=2)) / 15
    #maxowidth = np.min(np.sum(condition, axis=3)) / 26
    maxxwidth = 5
    maxywidth = 5
    maxzwidth = 5
    maxowidth = 3
    
    maxw = np.array([maxxwidth, maxywidth, maxzwidth, maxowidth])
    A = np.cumsum(np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2), axis=3)
    Cn, kaves, delstas = calc_cost4d(A, maxw, condition)

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
    k, kcond = calc_hist4d(A, opt_indx[0], opt_indx[1], opt_indx[2], opt_indx[3], condition) 

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k[:, 0, 0, :]), vmax=np.max(k[:, 0, 0, :]), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data[:, 0, 0, :]), vmax=np.max(data[:, 0, 0, :]), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)


def run_tst3d():
    datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[:, :, 0, :]
    condition = np.ones(data.shape, dtype=bool)
    
    A = np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2)

    k, kcond = calc_hist3d(A, 2, 2, 2, condition) 
    #for i in range(0,100):
    #    print i,np.average(k[:,i,:])

    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(k[:,:,1]), vmax=np.max(k[:,:,1]), cmap='jet')
    mappable = make_mappable(np.max(k))
    plt.colorbar(mappable)
    plt.figure(figsize=(8, 16))
    plt.pcolor(np.transpose(data[:, :, 1]), vmax=np.max(data[:, :, 1]), cmap='jet')

    mappable = make_mappable(np.max(data))
    plt.colorbar(mappable)


def run_tst4d():
    #datafile = "/home/kazu/cscl/phonopy_222/m200200200/data3.hdf5"
    datafile = "data3_100000000.hdf5"
    f = h5py.File(datafile)
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    condition = np.ones(data.shape, dtype=bool)
    
    A = np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2)
    print "cumsum finished"

    k, kcond = calc_hist4d(A, 6, 6, 6, 6, condition) 
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


#run_tst()
#run2d()
#runex()
#run_simu()
#run_simu3d()

#run_tst()
#run2d()
#runex()
#run_simu()
run_simu4d()
#run_tst4d()
plt.show()
