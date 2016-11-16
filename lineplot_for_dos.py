#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

data=np.loadtxt('total_dos.dat',comments='#',dtype='float')
adata=np.loadtxt('/home/kazu/asi3n4/phono3py_112_fc2_222_sym_monk_shift/total_dos.dat',comments='#',dtype='float')
gdata=np.loadtxt('/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/total_dos.dat',comments='#',dtype='float')

print data.shape

i=0
for x in data[:,0]:
	if x <= 7:
		i=i+1
ai=0
for x in adata[:,0]:
	if x <= 5:
		ai=ai+1
gi=0
for x in gdata[:,0]:
	if x <= 6:
		gi=gi+1


def func(x,a):
	return a*x*x

parameter_optimal, covaiance = scipy.optimize.curve_fit(func,data[0:i,0],data[0:i,1])
print "parameter =", parameter_optimal
y = func(data[0:i,0],parameter_optimal[0])
plt.figure(figsize=(20,20))
plt.subplot(3,2,1)
plt.plot(data[:,0],data[:,1],label="total")
plt.plot(data[0:i,0],y,label="a*omega**2")
plt.xlim(0,35)
plt.ylim(0,11)
#plt.yticks([0,0.4,0.8,1.2])
plt.legend(loc='upper left')

parameter_optimal, covaiance = scipy.optimize.curve_fit(func,adata[0:ai,0],adata[0:ai,1])
print "parameter =", parameter_optimal
y = func(adata[0:ai,0],parameter_optimal[0])
plt.subplot(3,2,3)
plt.plot(adata[:,0],adata[:,1],label="atotal")
plt.plot(adata[0:ai,0],y,label="a*omega**2")
plt.xlim(0,35)
plt.ylim(0,22)
plt.yticks([0,4,8,12,16,20])
plt.legend(loc='upper left')

parameter_optimal, covaiance = scipy.optimize.curve_fit(func,gdata[0:gi,0],gdata[0:gi,1])
print "parameter =", parameter_optimal
y = func(gdata[0:gi,0],parameter_optimal[0])
plt.subplot(3,2,5)
plt.plot(gdata[:,0],gdata[:,1],label="gtotal")
plt.plot(gdata[0:gi,0],y,label="a*omega**2")
plt.xlim(0,35)
plt.ylim(0,11)
#plt.yticks([0,0.4,0.8,1.2])
plt.legend(loc='upper left')


plt.xlabel("Omega [THz]")
plt.ylabel("DOS [1/THz]")
#plt.show()
plt.savefig("tdos-three-phases.eps")
