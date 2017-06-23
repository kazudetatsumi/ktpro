#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py

f = "kappa-m8820.hdf5"
g = "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/kappa-m121212.hdf5"
h = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/kappa-m8810.hdf5"
mytemp = 300

def parse_gamma(filename,temp):
    f = h5py.File(filename,'r')
    temperature = f["temperature"].value
    i=0 
    for t in temperature:
        if t == temp:
            tindex=i 
        i += 1
    gamma=f["gamma"][tindex,:,:]
    omega=f["frequency"][:,:]
    gammashape=gamma.shape
    gamma1=gamma.reshape( (gammashape[0]*gammashape[1]) )
    omega1=omega.reshape( (gammashape[0]*gammashape[1]) )
    dgamma=[]
    domega=[]
    j=0
    for i in omega1:
    	if i <= 12:
    		dgamma.append(gamma1[j])
    		domega.append(i)
    	j=j+1
    print np.average(dgamma),np.std(dgamma)
    return(omega1,gamma1) 

def plotset():
   plt.title("tau vs omega")
   plt.legend(loc='upper right')
   plt.xlim(0,12)
   plt.ylim(0,0.12)
   plt.yticks([0,0.03,0.06,0.09,0.12])
   plt.xlabel("Omega. [THz]")
   plt.ylabel("Tau [ns]")


def run():
   omega1,gamma1=parse_gamma(f,mytemp) 
   gomega1,ggamma1=parse_gamma(g,mytemp) 
   aomega1,agamma1=parse_gamma(h,mytemp) 

   plt.figure(figsize=(4,12.5))

   plt.subplot(3,1,1)
   plt.plot(omega1,1/(2*gamma1*2*np.pi)*0.001,'.',label="beta",markersize=6,c='blue',fillstyle='full')
   plotset()

   plt.subplot(3,1,2)
   plt.plot(gomega1,1/(2*ggamma1*2*np.pi)*0.001,'.',label="gamma",markersize=6,c='green',fillstyle='none')
   plotset()

   plt.subplot(3,1,3)
   plt.plot(aomega1,1/(2*agamma1*2*np.pi)*0.001,'.',label="alpha",markersize=6,c='red',fillstyle='none',linewidth=1)
   plotset()

   #plt.show()
   plt.savefig("tau-omega2_fc2_338.eps")


run()
