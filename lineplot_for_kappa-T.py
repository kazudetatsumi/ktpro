#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import csv


c = h5py.File("/home/kazu/bc3n4_m/phono3py_113_fc2_224_sym/kappa-m141432.hdf5")
kxx_c = c['kappa'][:,0]
kzz_c = c['kappa'][:,2]
t = c['temperature'][:]

s = h5py.File("/home/kazu/bsi3n4_m/phono3py_113_fc2_224_sym_monk_shift/kappa-m141432.hdf5")
kxx_s = s['kappa'][:,0]
kzz_s = s['kappa'][:,2]

g = h5py.File("/home/kazu/bge3n4_m/phono3py_113_fc2_224_sym/kappa-m141432.hdf5")
kxx_g = g['kappa'][:,0]
kzz_g = g['kappa'][:,2]








plt.figure(figsize=(8,6))
plt.title("kappa vs temp")
plt.plot(t,kxx_c,"c",label="kxx_c")
plt.plot(t,kzz_c,"c",label="kzz_c")
plt.plot(t,kxx_s,"b",label="kxx_s")
plt.plot(t,kzz_s,"b",label="kzz_s")
plt.plot(t,kxx_g,"r",label="kxx_g")
plt.plot(t,kzz_g,"r",label="kzz_g")

plt.legend()
plt.xlim(0,1000)
#plt.ylim(0,300)
#plt.ylim(0,700)
plt.xlabel("Temp. [K]")
plt.ylabel("Kappa [W/K.m]")
plt.yscale("log")
#plt.yticks([5,10,20,40,80,160])
#plt.xscale("log")
#plt.show()
plt.savefig("kappa-temp_csige.eps")
