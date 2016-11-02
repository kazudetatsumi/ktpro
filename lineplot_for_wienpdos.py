#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import csv

data001=np.loadtxt('partial_dos_001.dat',comments='#',dtype='float')
data100=np.loadtxt('partial_dos_100.dat',comments='#',dtype='float')
adata001=np.loadtxt('/home/kazu/asi3n4/phono3py_112_fc2_222_sym_monk_shift/partial_dos_001.dat',comments='#',dtype='float')
adata100=np.loadtxt('/home/kazu/asi3n4/phono3py_112_fc2_222_sym_monk_shift/partial_dos_100.dat',comments='#',dtype='float')

print data001.shape



plt.figure(figsize=(20,20))
plt.subplot(6,2,1)
plt.plot(data001[:,0],data001[:,1],label="Si1_001")
plt.plot(data100[:,0],data100[:,1],label="Si1_100")
plt.xlim(0,12)
plt.ylim(0,0.2)
plt.legend(loc='upper left')
plt.subplot(6,2,3)
plt.plot(data001[:,0],data001[:,7],label="N7_001")
plt.plot(data100[:,0],data100[:,7],label="N7_100")
plt.xlim(0,12)
plt.ylim(0,0.2)
plt.legend()
plt.subplot(6,2,5)
plt.plot(data001[:,0],data001[:,10],label="N10_001")
plt.plot(data100[:,0],data100[:,10],label="N10_100")
plt.legend(loc='upper left')
plt.xlim(0,12)
plt.ylim(0,0.2)

plt.subplot(6,2,2)
plt.plot(adata001[:,0],adata001[:,1],label="aSi1_001")
plt.plot(adata100[:,0],adata100[:,1],label="aSi1_100")
plt.xlim(0,12)
plt.ylim(0,0.2)
plt.legend(loc='upper left')

plt.subplot(6,2,4)
plt.plot(adata001[:,0],adata001[:,4],label="aSi4_001")
plt.plot(adata100[:,0],adata100[:,4],label="aSi4_100")
plt.xlim(0,12)
plt.ylim(0,0.2)
plt.legend(loc='upper left')

plt.subplot(6,2,6)
plt.plot(adata001[:,0],adata001[:,13],label="aN13_001")
plt.plot(adata100[:,0],adata100[:,13],label="aN13_100")
plt.xlim(0,12)
plt.ylim(0,0.2)
plt.legend(loc='upper left')

plt.subplot(6,2,8)
plt.plot(adata001[:,0],adata001[:,14],label="aN14_001")
plt.plot(adata100[:,0],adata100[:,14],label="aN14_100")
plt.xlim(0,12)
plt.ylim(0,0.2)
plt.legend(loc='upper left')

plt.subplot(6,2,10)
plt.plot(adata001[:,0],adata001[:,17],label="aN17_001")
plt.plot(adata100[:,0],adata100[:,17],label="aN17_100")
plt.xlim(0,12)
plt.ylim(0,0.2)
plt.legend(loc='upper left')

plt.subplot(6,2,12)
plt.plot(adata001[:,0],adata001[:,20],label="aN20_001")
plt.plot(adata100[:,0],adata100[:,20],label="aN20_100")
plt.xlim(0,12)
plt.ylim(0,0.2)
plt.legend(loc='upper left')


plt.xlabel("Omega [THz]")
plt.ylabel("PDOS [1/THz]")
#plt.show()
plt.savefig("pdos-twoe-phases-xyz.eps")
