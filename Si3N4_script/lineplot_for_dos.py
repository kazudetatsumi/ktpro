#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

cdata=np.loadtxt('/home/kazu/bc3n4_m/phono3py_113_fc2_224_sym/total_dos.dat',comments='#',dtype='float')
sdata=np.loadtxt('/home/kazu/bsi3n4_m/phono3py_113_fc2_224_sym_monk_shift/total_dos.dat',comments='#',dtype='float')
gdata=np.loadtxt('/home/kazu/bge3n4_m/phono3py_113_fc2_224_sym/total_dos.dat',comments='#',dtype='float')

print cdata.shape

ci=0
for x in cdata[:,0]:
	if x <= 2.5:
		ci=ci+1
si=0
for x in sdata[:,0]:
	if x <= 2.5:
		si=si+1
gi=0
for x in gdata[:,0]:
	if x <= 2.5:
		gi=gi+1


def func(x,a):
	return a*x*x

parameter_optimal, covaiance = scipy.optimize.curve_fit(func,cdata[0:ci,0],cdata[0:ci,1])
print "parameter =", parameter_optimal
y = func(cdata[0:ci,0],parameter_optimal[0])
plt.figure(figsize=(6,3))
#plt.subplot(3,2,1)
plt.plot(cdata[:,0],cdata[:,1]/2,label="bC3N4")
plt.plot(sdata[:,0],sdata[:,1]/4,label="bSi3N4")
plt.plot(gdata[:,0],gdata[:,1]/2,label="bGe3N4")
#plt.plot(data[0:i,0],y,label="a*omega**2")
plt.xlim(0,42)
plt.ylim(0,6.0)
#plt.yticks([0,0.4,0.8,1.2])
plt.legend(loc='upper left')

parameter_optimal, covaiance = scipy.optimize.curve_fit(func,sdata[0:si,0],sdata[0:si,1])
print "parameter =", parameter_optimal
#y = func(adata[0:ai,0],parameter_optimal[0])
#plt.subplot(3,2,3)
#plt.plot(adata[:,0],adata[:,1],label="atotal")
#plt.plot(adata[0:ai,0],y,label="a*omega**2")
#plt.xlim(0,35)
#plt.ylim(0,22)
#plt.yticks([0,4,8,12,16,20])
#plt.legend(loc='upper left')

parameter_optimal, covaiance = scipy.optimize.curve_fit(func,gdata[0:gi,0],gdata[0:gi,1])
print "parameter =", parameter_optimal
#y = func(gdata[0:gi,0],parameter_optimal[0])
#plt.subplot(3,2,5)
#plt.plot(gdata[:,0],gdata[:,1],label="gtotal")
#plt.plot(gdata[0:gi,0],y,label="a*omega**2")
#plt.xlim(0,35)
#plt.ylim(0,11)
#plt.yticks([0,0.4,0.8,1.2])
#plt.legend(loc='upper left')


plt.xlabel("Omega [THz]")
plt.ylabel("DOS [1/THz]")
#plt.show()
plt.savefig("tdos-CSiGe.eps")
