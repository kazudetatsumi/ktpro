#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def gauss(x):
    return np.exp(-x**2/2.0)/(2.0*np.pi)**0.5

def hist(dx):
    binedges = np.arange(-3, 3, dx)
    binheights = gauss(binedges)
    return binheights, binedges

def calcise(x, y, binheights, binedges):
    dxx = x[1] - x[0]
    dx = binedges[1] - binedges[0]
    ise = 0.0
    for be, bh in zip(binedges, binheights):
        for _x, _y in zip(x, y):
            if _x >= be-dx*0.5 and _x < be+dx*0.5:
                ise += (_y - bh)**2*dxx
    return(ise)


def run():
    x = np.arange(-3, 3, 0.0001)
    y = gauss(x)
    ise = []
    binwidth = []
    for dx in np.arange(2.0, 0.02, -0.02):
        binheights, binedges = hist(dx)
        ise.append(calcise(x, y, binheights, binedges))
        binwidth.append(dx)
        print(calcise(x, y, binheights, binedges), dx)
        #plt.plot(x,y)
        #plt.bar(binedges, binheights, dx, fc=None, fill=False)
        #plt.show()
    plt.plot(binwidth, ise)
    plt.savefig('binwidth_ise.pdf')
    plt.show()



run()
