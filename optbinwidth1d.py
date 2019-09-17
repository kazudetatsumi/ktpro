#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def getdata(f):
    data = np.genfromtxt(f,  delimiter=',', dtype=None)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    dx = 0.005
    dy = 0.1
    xlin = np.arange(min(x), max(x), dx)
    nx = xlin.shape[0]
    ylin = np.arange(min(y), max(y), dy)
    ny = ylin.shape[0]
    
    karr = np.zeros((nx, ny))

    for _x, _y, _z in zip(x, y, z):
       xx =  np.where(abs(xlin - _x) < 0.0000001)
       yy =  np.where(abs(ylin - _y) < 0.0000001)
       karr[xx, yy] = _z

    return np.sum(karr[390:400,:], axis=0)


def calc_hist(cumdata, nw):
    Nmax = cumdata.shape[0]
    N = int((Nmax - (Nmax % nw)) / nw)  
    k = np.zeros((N))
    for i in range(0, N):
        ihead = (i+1)*nw - 1
        if i == 0:
            k[i] = cumdata[ihead]
        else:
            k[i] = cumdata[ihead] - cumdata[ihead - nw]
    return k
        

def calc_cost(cumdata, maxw):
    Cn = np.zeros((maxw))
    for i in range(2, maxw):
       k = calc_hist(cumdata, i)
       kave = np.average(k)
       v = np.var(k)
       cost = (2 * kave - v) / (i**2)
       Cn[i] = cost
    return Cn


def run1d():
    txtfile = "/home/kazu/data/5min_fine.txt"
    data = getdata(txtfile)
    maxw =  int(data.shape[0] / 2)
    print data[0:10]
    cumdata = np.cumsum(data)
    print cumdata[0:10]

    Cn = calc_cost(cumdata, maxw)
    print np.argmin(Cn)
    k = calc_hist(cumdata, 22)

    plt.figure(figsize=(16, 8))
    plt.plot(k)
    plt.figure(figsize=(16, 8))
    plt.plot(data)


run1d()
plt.show()
