#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'

homedir = "/home/kazu/"
gdataL0=np.loadtxt(homedir + '/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/jdos-m222222-g0-t300.dat',comments='#',dtype='float')
gdataL1=np.loadtxt(homedir + '/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/jdos-m222222-g2420-05L-t300.dat',comments='#',dtype='float')
gdataL2=np.loadtxt(homedir + '/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/jdos-m222222-g5324-L-t300.dat',comments='#',dtype='float')
gz=2


plt.figure(figsize=(9,12))
plt.subplot(6,1,1)
plt.plot(gdataL0[:,0],gdataL0[:,1]/gz**2,color="green",label="gSi3N4_1")
plt.plot(gdataL0[:,0],gdataL0[:,2]/gz**2,':',color="green",label="gSi3N4_2")
#plt.plot(data[0:i,0],y,label="a*omega**2")
#plt.xlim(0,42)
plt.ylim(0,28)
plt.yticks([0,10,20])
#plt.yticks([0,0.4,0.8,1.2])
#plt.legend(loc='upper right')



plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")
plt.subplot(6,1,2)
plt.plot(gdataL1[:,0],gdataL1[:,1]/gz**2,color="green",label="gSi3N4_1")
plt.plot(gdataL1[:,0],gdataL1[:,2]/gz**2,':',color="green",label="gSi3N4_2")
plt.ylim(0,28)
plt.yticks([0,10,20])
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.subplot(6,1,3)
plt.plot(gdataL2[:,0],gdataL2[:,1]/gz**2,color="green",label="gSi3N4_1")
plt.plot(gdataL2[:,0],gdataL2[:,2]/gz**2,':',color="green",label="gSi3N4_2")
plt.ylim(0,28)
plt.yticks([0,10,20])
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.subplot(6,1,4)
#plt.plot(cdataK0[:,0],cdataK0[:,1]/cz**2,color="red",label="aSi3N4_1")
#plt.plot(cdataK0[:,0],cdataK0[:,2]/cz**2,':',color="red",label="aSi3N4_2")
#plt.plot(sdataK0[:,0],sdataK0[:,1]/sz**2,color="blue",label="bSi3N4_1")
#plt.plot(sdataK0[:,0],sdataK0[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
plt.ylim(0,28)
plt.yticks([0,10,20])
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.subplot(6,1,5)
#lt.plot(cdataK1[:,0],cdataK1[:,1]/cz**2,color="red",label="aSi3N4_1")
#lt.plot(cdataK1[:,0],cdataK1[:,2]/cz**2,':',color="red",label="aSi3N4_2")
#lt.plot(sdataK1[:,0],sdataK1[:,1]/sz**2,color="blue",label="bSi3N4_1")
#lt.plot(sdataK1[:,0],sdataK1[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
plt.ylim(0,28)
plt.yticks([0,10,20])
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.subplot(6,1,6)
#lt.plot(cdataK2[:,0],cdataK2[:,1]/cz**2,color="red",label="aSi3N4_1")
#lt.plot(cdataK2[:,0],cdataK2[:,2]/cz**2,':',color="red",label="aSi3N4_2")
#lt.plot(sdataK2[:,0],sdataK2[:,1]/sz**2,color="blue",label="bSi3N4_1")
#lt.plot(sdataK2[:,0],sdataK2[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
plt.ylim(0,28)
plt.yticks([0,10,20])
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

#plt.show()
plt.savefig("wjdos-gSi3N4.eps")
