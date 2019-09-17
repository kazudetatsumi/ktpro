#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def get2ddata(f):
    data = np.genfromtxt(f,  delimiter=',', dtype=None)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    dx = 0.005
    dy = 0.1
    xlin = np.arange(min(x), max(x)+dx, dx)
    nx = xlin.shape[0]
    ylin = np.arange(min(y), max(y)+dy, dy)
    ny = ylin.shape[0]
    print nx, ny
    
    karr = np.zeros((nx, ny))

    for _x, _y, _z in zip(x, y, z):
       xx =  np.where(abs(xlin - _x) < 0.0000001)
       yy =  np.where(abs(ylin - _y) < 0.0000001)
       #print "check xx yy", xx, yy
       #karr[xx, yy] = _z
       karr[xx, yy] = _z + 0.00000000001

    condition = karr > 0.000000000001
    karrnonzero = np.extract(condition, karr)
    if karrnonzero.shape[0] != x.shape[0]:
         print "num of nonzero karr is not num of data", karrnonzero.shape
    else:
         print "num of nonzero karr matches  num of data", karrnonzero.shape


    #return karr[100:225,75:]
    #data = karr[100:225, 75:]
    #print "orig_data",data.shape
    #data = np.concatenate((data, np.zeros_like(data)), axis = 0)
    #data = np.concatenate((data, np.zeros_like(data)), axis = 0)
    #print "dummy_zero", data.shape
    #return karr[20:685, 105:]
    #return karr[100:560, 0:300]
    #return karr[:, 105:]
    #return karr[75:450, :]
    return karr[10:690, 115:]
    #return data



def calc_hist2d(A, nw0, nw1):
    Nmax = A.shape
    N0 = int((Nmax[0] - (Nmax[0] % nw0)) / nw0)  
    N1 = int((Nmax[1] - (Nmax[1] % nw1)) / nw1)  
    k = np.zeros((N0, N1))
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
    return k
        

def calc_cost2d(A, maxw):
    Cn = np.zeros((maxw))
    for i in range(1, maxw[0]):
       for j in range(1, maxw[1]):
          k = calc_hist2d(A, i, j)
          condition = k >= 0.0000000001
          knonzero = np.extract( condition , k)   # eliminate zero value element 
          if i == 1 and j ==1:
             print "shape of k matrix with zero elements",k.shape
             print "number of nonzero k",knonzero.shape
          kave = np.average(knonzero)
          v = np.var(knonzero)
          #kave = np.average(k)
          #v = np.var(k)
          cost = (2 * kave - v) / ((i*j)**2)
          Cn[i, j] = cost
    return Cn


def run2d():
    txtfile = "/home/kazu/data/20min_fine.txt"
    data = get2ddata(txtfile)
    maxw =  np.array([int(data.shape[0] / 2), int(data.shape[1]) / 2])
    print maxw
    #print data[0:10]
    A = np.cumsum(np.cumsum(data, axis=0), axis=1)

    Cn = calc_cost2d(A, maxw)
    print "opt bin index", np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    k = calc_hist2d(A, opt_indx[0], opt_indx[1]) 
    #k = calc_hist2d(A, 4, 4) 

    plt.figure(figsize=(16, 8))
    #plt.plot(np.sum(data, axis=0))
    #plt.pcolor(data, vmax=np.max(data)/10)
    plt.pcolor(np.transpose(k), vmax=np.max(k), cmap='jet')
    plt.figure(figsize=(16, 8))
    plt.pcolor(np.transpose(data), vmax=np.max(data)/100, cmap='jet')


def run_tst():
    txtfile = "/home/kazu/data/20min_fine.txt"
    data = get2ddata(txtfile)

run_tst()
#run2d()
#plt.show()
