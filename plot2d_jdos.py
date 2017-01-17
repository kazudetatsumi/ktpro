#!/usr/bin/env python
from matplotlib.colors import LinearSegmentedColormap

import sys,os,glob,string,array,math,h5py
import numpy as np
from scipy.interpolate import griddata

from matplotlib.colors import LinearSegmentedColormap

import matplotlib.pyplot as plt
import matplotlib.cm as cm


f = '/home/kazu/bsi3n4_m/phono3py_113_fc2_224_sym_monk_shift/log-jdos-g-a'
jdoshead = '/home/kazu/bsi3n4_m/phono3py_113_fc2_224_sym_monk_shift/jdos-m141432-g'

def parse_gp(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if "Grid point:" in line:
                values = line.split() 
                data.append(values[2])
    return np.array(data)

def getjdos(a):
    i = 0
    tmpfile = jdoshead+str(a[0])+".dat"
    tmpdata = np.loadtxt(tmpfile)
    tmpsize = tmpdata.shape
    ngpsize = a.shape
    ngp = ngpsize[0]
    print ngp
    print tmpsize[0]
    print tmpsize[1]
    jdata = np.zeros((ngp, tmpsize[0], tmpsize[1]))
    for x in np.nditer(a): 
       infile = jdoshead+str(x)+".dat"
       print infile
       jdata[i,:,:] = np.loadtxt(infile)
       i += 1
    return(jdata)
gp = parse_gp(f)
jdos = getjdos(gp)


tmpsize = jdos.shape
print tmpsize

#x = np.zeros((tmpsize[0],tmpsize[1]))
#y = np.zeros((tmpsize[0],tmpsize[1]))
k=0
Z = np.zeros((tmpsize[0],tmpsize[1]))
for i in range(0,tmpsize[0]):
   for j in range(0,tmpsize[1]):
      Z[i,j] = jdos[i,j,1] + jdos[i,j,2]
x = np.arange(tmpsize[1])
y = np.arange(tmpsize[0])
X, Y = np.meshgrid(x, y)

plt.figure(figsize=(6,3))
plt.pcolor(X,Y,Z,vmin=np.min(Z),vmax=np.max(Z),cmap=cm.bwr)
plt.colorbar(ticks=[np.min(Z),np.max(Z)])
plt.show()
