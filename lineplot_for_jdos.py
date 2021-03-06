#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

homedir = "/home/kazu/stern/"
cdata=np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g0.dat',comments='#',dtype='float')
sdata=np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g0.dat',comments='#',dtype='float')
gdata=np.loadtxt(homedir + '/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/jdos-m222222-g0.dat',comments='#',dtype='float')
cdosdata=np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/total_dos_m141416.dat',comments='#',dtype='float')
sdosdata=np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/total_dos_m141432.dat',comments='#',dtype='float')
gdosdata=np.loadtxt(homedir + '/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/total_dos.dat',comments='#',dtype='float')


plt.figure(figsize=(9,8))
plt.subplot(2,1,1)
plt.plot(cdata[:,0],cdata[:,1]/4.0,color="red",label="aSi3N4_1")
plt.plot(cdata[:,0],cdata[:,2]/4.0,':',color="red",label="aSi3N4_2")
plt.plot(sdata[:,0],sdata[:,1],color="blue",label="bSi3N4_1")
plt.plot(sdata[:,0],sdata[:,2],':',color="blue",label="bSi3N4_2")
plt.plot(gdata[:,0],gdata[:,1],color="green",label="gSi3N4_1")
plt.plot(gdata[:,0],gdata[:,2],':',color="green",label="gSi3N4_2")
#plt.plot(data[0:i,0],y,label="a*omega**2")
#plt.xlim(0,42)
#plt.ylim(0,6.0)
#plt.yticks([0,0.4,0.8,1.2])
plt.legend(loc='upper right')



plt.xlabel("Omega [THz]")
plt.ylabel("JDOS [1/THz]")
plt.subplot(2,1,2)
plt.plot(cdosdata[:,0],cdosdata[:,1]/2.0,color="red",label="aSi3N4_dos")
plt.plot(sdosdata[:,0],sdosdata[:,1],color="blue",label="bSi3N4_dos")
plt.plot(gdosdata[:,0],gdosdata[:,1],color="green",label="gSi3N4_dos")
plt.xlim(0,70)
plt.legend(loc='upper right')
plt.show()
#plt.savefig("jdos-aSi3N4_bSi3N4.eps")
