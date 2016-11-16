#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py


f = h5py.File("kappa-m8820.hdf5")
#kxx = f['kappa'][:,0]
#kzz = f['kappa'][:,2]
#t = f['temperature'][:]

g = h5py.File("/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/kappa-m121212.hdf5")

h = h5py.File("/home/kazu/asi3n4/phono3py_112_fc2_222_sym_monk_shift/kappa-m8810.hdf5")

gamma=f['gamma'][30,:,:]
omega=f['frequency'][:,:]
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


ggamma=g['gamma'][30,:,:]
gomega=g['frequency'][:,:]
ggammashape=ggamma.shape
ggamma1=ggamma.reshape( (ggammashape[0]*ggammashape[1]) )
gomega1=gomega.reshape( (ggammashape[0]*ggammashape[1]) )

dggamma=[]
dgomega=[]
j=0
for i in gomega1:
	if i <= 12:
		dggamma.append(ggamma1[j])
		dgomega.append(i)
	j=j+1
print np.average(dggamma),np.std(dggamma)

agamma=h['gamma'][30,:,:]
aomega=h['frequency'][:,:]
agammashape=agamma.shape
agamma1=agamma.reshape( (agammashape[0]*agammashape[1]) )
aomega1=aomega.reshape( (agammashape[0]*agammashape[1]) )

dagamma=[]
daomega=[]
j=0
for i in aomega1:
	if i <= 12:
		dagamma.append(agamma1[j])
		daomega.append(i)
	j=j+1
print np.average(dagamma),np.std(dagamma)

plt.figure(figsize=(4,12))
plt.subplot(3,1,1)
plt.title("tau vs omega")
plt.plot(omega1,1/gamma1*0.001,'.',label="beta",markersize=3,c='blue',fillstyle='none')
#plt.plot(gomega1,1/ggamma1*0.001,'.',label="gamma",markersize=3,c='green',fillstyle='none')
#plt.plot(aomega1,1/agamma1*0.001,'.',label="alpha",markersize=3,c='red',fillstyle='none')
plt.legend(loc='upper right')
plt.xlim(0,12)
plt.ylim(0,1.5)
plt.yticks([0,0.5,1.0,1.5])
#plt.yticks([0,0.01,0.02,0.03])
plt.xlabel("Omega. [THz]")
plt.ylabel("Tau [ns]")


plt.subplot(3,1,2)
plt.title("tau vs omega")
plt.plot(gomega1,1/ggamma1*0.001,'.',label="gamma",markersize=3,c='green',fillstyle='none')
plt.legend(loc='upper right')
plt.xlim(0,12)
plt.ylim(0,1.5)
plt.yticks([0,0.5,1.0,1.5])
plt.xlabel("Omega. [THz]")
plt.ylabel("Tau [ns]")


plt.subplot(3,1,3)
plt.title("tau vs omega")
plt.plot(aomega1,1/agamma1*0.001,'.',label="alpha",markersize=3,c='red',fillstyle='none')
plt.legend(loc='upper right')
plt.xlim(0,12)
plt.ylim(0,1.5)
plt.yticks([0,0.5,1.0,1.5])
plt.xlabel("Omega. [THz]")
plt.ylabel("Tau [ns]")
#plt.show()
plt.savefig("tau-omega2.eps")
