#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys,os,glob,string,array,math,h5py
import numpy as np  

f = h5py.File("co_l23_ddscs.hdf5")
ddscs = f['ddscs']
#print ddscs[3,0,:,0]
cddscs = np.sum(ddscs,axis=2)
tmpsize = cddscs.shape
print tmpsize

    
#plt.figure(figsize=(6,3))
F = ( cddscs[3,1,:] - cddscs[3,0,:] ) / ( cddscs[0,0,:] + cddscs[0,1,:] )
#plt.plot(F,label="tst")
#plt.show()
np.savetxt('out.csv',F)
