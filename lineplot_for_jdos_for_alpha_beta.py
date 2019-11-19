#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

homedir = "/home/kazu/"
cdataA0=np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g0-A.dat',comments='#',dtype='float')
sdataA0=np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g0-A.dat',comments='#',dtype='float')
cdataA1=np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g784-05A.dat',comments='#',dtype='float')
sdataA1=np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g1568-05A.dat',comments='#',dtype='float')
cdataA2=np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g1568-A.dat',comments='#',dtype='float')
sdataA2=np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g3136-A.dat',comments='#',dtype='float')
cdataK0=np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g0-K.dat',comments='#',dtype='float')
sdataK0=np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g0-K.dat',comments='#',dtype='float')
cdataK1=np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g56-05K.dat',comments='#',dtype='float')
sdataK1=np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g56-05K.dat',comments='#',dtype='float')
cdataK2=np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g98-K.dat',comments='#',dtype='float')
sdataK2=np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g98-K.dat',comments='#',dtype='float')
cz=4
sz=2


plt.figure(figsize=(9,12))
plt.subplot(6,1,1)
plt.plot(cdataA0[:,0],cdataA0[:,1]/cz**2,color="red",label="aSi3N4_1")
plt.plot(cdataA0[:,0],cdataA0[:,2]/cz**2,':',color="red",label="aSi3N4_2")
plt.plot(sdataA0[:,0],sdataA0[:,1]/sz**2,color="blue",label="bSi3N4_1")
plt.plot(sdataA0[:,0],sdataA0[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
#plt.plot(data[0:i,0],y,label="a*omega**2")
#plt.xlim(0,42)
plt.ylim(0,50)
#plt.yticks([0,0.4,0.8,1.2])
#plt.legend(loc='upper right')



plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")
plt.subplot(6,1,2)
plt.plot(cdataA1[:,0],cdataA1[:,1]/cz**2,color="red",label="aSi3N4_1")
plt.plot(cdataA1[:,0],cdataA1[:,2]/cz**2,':',color="red",label="aSi3N4_2")
plt.plot(sdataA1[:,0],sdataA1[:,1]/sz**2,color="blue",label="bSi3N4_1")
plt.plot(sdataA1[:,0],sdataA1[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
plt.ylim(0,50)
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.subplot(6,1,3)
plt.plot(cdataA2[:,0],cdataA2[:,1]/cz**2,color="red",label="aSi3N4_1")
plt.plot(cdataA2[:,0],cdataA2[:,2]/cz**2,':',color="red",label="aSi3N4_2")
plt.plot(sdataA2[:,0],sdataA2[:,1]/sz**2,color="blue",label="bSi3N4_1")
plt.plot(sdataA2[:,0],sdataA2[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
plt.ylim(0,50)
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.subplot(6,1,4)
plt.plot(cdataK0[:,0],cdataK0[:,1]/cz**2,color="red",label="aSi3N4_1")
plt.plot(cdataK0[:,0],cdataK0[:,2]/cz**2,':',color="red",label="aSi3N4_2")
plt.plot(sdataK0[:,0],sdataK0[:,1]/sz**2,color="blue",label="bSi3N4_1")
plt.plot(sdataK0[:,0],sdataK0[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
plt.ylim(0,50)
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.subplot(6,1,5)
plt.plot(cdataK1[:,0],cdataK1[:,1]/cz**2,color="red",label="aSi3N4_1")
plt.plot(cdataK1[:,0],cdataK1[:,2]/cz**2,':',color="red",label="aSi3N4_2")
plt.plot(sdataK1[:,0],sdataK1[:,1]/sz**2,color="blue",label="bSi3N4_1")
plt.plot(sdataK1[:,0],sdataK1[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
plt.ylim(0,50)
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.subplot(6,1,6)
plt.plot(cdataK2[:,0],cdataK2[:,1]/cz**2,color="red",label="aSi3N4_1")
plt.plot(cdataK2[:,0],cdataK2[:,2]/cz**2,':',color="red",label="aSi3N4_2")
plt.plot(sdataK2[:,0],sdataK2[:,1]/sz**2,color="blue",label="bSi3N4_1")
plt.plot(sdataK2[:,0],sdataK2[:,2]/sz**2,':',color="blue",label="bSi3N4_2")
plt.ylim(0,50)
#plt.legend(loc='upper right')
plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")

plt.show()
#plt.savefig("jdos-aSi3N4_bSi3N4.eps")