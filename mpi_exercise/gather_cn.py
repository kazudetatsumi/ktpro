#!/usr/bin/env python
import h5py
import numpy as np

psize = 32

for i in range(0, psize):
    f=h5py.File("Cn_rank"+str(i)+".hdf5")
    Cn_tmp = f["Cn"]
    if i == 0:
        size = Cn_tmp.shape
        Cn_tot = np.zeros((psize, size[0], size[1], size[2], size[3]))
        Cn_merged = np.zeros((size[0], size[1], size[2], size[3]))
    Cn_tot[i, :, :, :, :] = Cn_tmp[:, :, :, :]

CONFLICT_FOUND = 0
for i in range(0, psize):
    for j in range(0, size[0]):
        for k in range(0, size[1]):
            for l in range(0, size[2]):
                for m in range(0, size[3]):
                    if i == 0 and j ==0 and k == 0 and l ==0 and m == 0:
                         Cn_merged[j,k,l,m] = Cn_tot[i,j,k,l,m]
                    elif i != 0 and j == 0 and k == 0 and l == 0 and m == 0:
                         print("i, j, k, l, m:", i, j, k, l, m)
                    elif abs(Cn_tot[i,j,k,l,m]) > 0.0000000001 and abs(Cn_merged[j,k,l,m]) < 0.00000001:
                         Cn_merged[j,k,l,m] = Cn_tot[i,j,k,l,m]
                    elif abs(Cn_tot[i,j,k,l,m]) > 0.0000000001 and abs(Cn_merged[j,k,l,m]) > 0.00000001:
                         print("Conflict: Cn_tot[i,j,k,l,m]=",Cn_tot[i,j,k,l,m],"Cn_merged[j,k,l,m]=",Cn_merged[j,k,l,m], "at", i, j, k, l, m)
                         CONFLICT_FOUND = 1

if CONFLICT_FOUND == 0:
    print("Check finished. Cn_rankX.hdf5 files show no inconsistency")

f=h5py.File("Cn.hdf5")
Cn_serial = f["Cn"]
print(Cn_serial.shape)
print(Cn_merged.shape)

for i in range(0, size[0]):
    for j in range(0, size[1]):
        for k in range(0, size[2]):
            for l in range(0, size[3]):
                if abs(Cn_merged[i, j, k, l] - Cn_serial[i, j, k, l]) > 0.000000001:
                    print("Conflict: Cn_merged[i,j,k,l]=",Cn_merged[i,j,k,l], "Cn_serial[i,j,k,l]=",Cn_serial[i,j,k,l], "at", i, j, k, l)
                    CONFLICT_FOUND = 1

if CONFLICT_FOUND == 0:
    print("Check finished. Cn_merged has no significant differences from Cn_seiral")
