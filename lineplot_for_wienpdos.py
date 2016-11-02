#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import csv
import sys
argvs=sys.argv



data001=np.loadtxt(argvs[1],comments='#',dtype='float')

print data001.shape



#plt.figure(figsize=(20,20))
#plt.subplot(6,2,1)
plt.plot(data001[:,0],data001[:,1],label="total")
plt.plot(data001[:,0],data001[:,2],label="Mg")
plt.plot(data001[:,0],data001[:,3],label="O")
plt.legend(loc='upper right')
#plt.xlim(0,12)
#plt.ylim(0,0.2)


plt.xlabel("Energy [eV]")
plt.ylabel("DOS [1/eV]")
#plt.show()
plt.savefig(argvs[1]+".eps")
