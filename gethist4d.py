#!/usr/bin/env python
# This script generate a coarser 4D matrix than the input matrix.
# The input 4D matrix is considered as a numpy array.
# For example, an input matrix of Nx x Ny x Nz x Ne, with a parameter of nx, ny
# , nz, and ne, this script generates a coarser matrix of Nx//nx x Ny//ny x Nz/
# /nz x Ne//ne.
# I wrote this script for the visualization of the histogram with the
# bin-widths optimized by optbinwidth4d_wholefort.py.
# Kazuyoshi TATSUMI 2020/4/24
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes
lib = ctypes.CDLL("/home/kazu/ktpro/histfort4d.so")

def calc_hist4d_f90(A, data, nw0, nw1, nw2, nw3,  condition):

    class result(ctypes.Structure):
        _fields_ = [("len0", ctypes.c_int), ("len1", ctypes.c_int), ("len2", ctypes.c_int),  ("len3", ctypes.c_int), ("arr", ctypes.POINTER(ctypes.c_double))]

    lib.hist4d.restype = result
    lib.hist4d.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), 
                           ctypes.POINTER(ctypes.c_int),  ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]

    Nmax0 = A.shape[0]
    Nmax1 = A.shape[1]
    Nmax2 = A.shape[2]
    Nmax3 = A.shape[3]

    result = lib.hist4d(ctypes.byref(ctypes.c_int(nw0)), ctypes.byref(ctypes.c_int(nw1)), ctypes.byref(ctypes.c_int(nw2)), ctypes.byref(ctypes.c_int(nw3)),
                        ctypes.byref(ctypes.c_int(Nmax0)), ctypes.byref(ctypes.c_int(Nmax1)), ctypes.byref(ctypes.c_int(Nmax2)), ctypes.byref(ctypes.c_int(Nmax3)), A)
    result_len0 = result.len0
    result_len1 = result.len1
    result_len2 = result.len2
    result_len3 = result.len3
    result_vec = np.ctypeslib.as_array(result.arr, shape=(result_len0, result_len1, result_len2, result_len3, ))
    #print result_vec

    kcond = np.zeros((result_len0, result_len1, result_len2, result_len3))
    for i in range(0, result_len0):
        for j in range(0, result_len1):
            for h in range(0, result_len2):
                for l in range(0, result_len3):
                    kcond[i, j, h, l] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1, h*nw2:(h+1)*nw2, l*nw3:(l+1)*nw3])
    return result_vec, result_len0, result_len1, result_len2, result_len3, kcond


def save_hist(outfile, k, kcond):
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('data4', data=k)
        hf.create_dataset('condition', data=kcond)


def plot_crosssection(xi, xe, yi, ye, zi, ze, data4):
    fig = plt.figure(figsize=(6, 8))
    fig.suptitle("crosssections of 4D INS data #Ei24", fontsize="x-large")
    ax = fig.add_subplot(2, 1, 1)
    #ax.pcolor(np.transpose(data4[:, 42, 7, :]), vmax=np.max(data4[:, 42, 7, :])/198, cmap='jet')
    #tmpdata = np.transpose(np.sum(np.sum(data4[42:45, :, 15:18, :], axis=0), axis=1))
    #tmpdata = np.transpose(np.sum(np.sum(data4[172:176, :, 15:17, :], axis=0), axis=1))
    #tmpdata = np.transpose(data4[58, :, 8, :])
    tmpdata = np.transpose(np.sum(np.sum(data4[:, :, :, :], axis=3), axis=2))
    #ax.pcolor(tmpdata, vmax=np.max(tmpdata)/1.5, cmap='jet')
    ax.pcolor(tmpdata, vmax=np.max(tmpdata)/150, cmap='jet')
    ax.text(2, 48, 'qy=58, qz=8', color='white')
    ax.set_xlabel('qx')
    ax.set_ylabel('E')
    ax.axis('tight')
    ax.axhline(y=20, color='white', lw=0.5)
    ax.axhline(y=52, color='white', lw=0.5)
    ax.axvline(x=34, color='white', lw=0.5)
    ax.axvline(x=61, color='white', lw=0.5)
    ax.xaxis.set_label_coords(0.5, 1.145)
    ax.tick_params(direction="in", color="white", top=True, labeltop=True, labelbottom=False)
    ax = fig.add_subplot(2, 1, 2)
    #ax.pcolor(np.transpose(data4[73, :, 7, :]), vmax=np.max(data4[73, :, 7, :])/400, cmap='jet')
    #tmpdata = np.transpose(np.sum(np.sum(data4[:, :, :, :], axis=0), axis=0))
    #tmpdata = np.transpose(np.sum(np.sum(data4[:, :, :, :], axis=0), axis=1))
    tmpdata = np.transpose(np.sum(np.sum(data4[72:76, :, 5:10, :], axis=0), axis=1))
    ax.pcolor(tmpdata, vmax=np.max(tmpdata)/12, cmap='jet')
    #ax.pcolor(np.transpose(data4[40, :, 16, :]), vmax=np.max(data4[40, :, 16, :])/18, cmap='jet')
    ax.text(2, 38, 'qx=73, qz=6', color='white')
    ax.set_xlabel('qy')
    ax.set_ylabel('E')
    ax.axis('tight')
    ax.xaxis.set_label_coords(0.5, 1.145)
    ax.tick_params(direction="in", color="white", top=True, labeltop=True, labelbottom=False)
    fig.subplots_adjust(top=0.90)


def run():
    nx = 3
    ny = 2
    nz = 2
    ne = 4
    head ="./"
    #datafile = head + "out_hw_all.hdf5"
    datafile = head + "Output4D_00_840.hdf5"
    outfile = head + "hist.hdf5"
    f = h5py.File(datafile)
    data = f["data4"][:]*1.0 # nqx, nqy, nqz, nomega
    condition = np.ones(data.shape, dtype=bool)
    print(data.shape)
    
    A = np.cumsum(np.cumsum(np.cumsum(np.cumsum(data, axis=0), axis=1), axis=2), axis=3)
    print("cumsum finished")

    k, klen0, klen1, klen2, klen3, kcond = calc_hist4d_f90(A, data, nx, ny, nz, ne, condition)
    print("hist finished")

    save_hist(outfile, k, kcond)

    lib.delete_array4.restype = None
    lib.delete_array4.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                  ctypes.POINTER(ctypes.c_int), np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]
    lib.delete_array4(ctypes.byref(ctypes.c_int(klen0)), ctypes.byref(ctypes.c_int(klen1)), ctypes.byref(ctypes.c_int(klen2)),
                      ctypes.byref(ctypes.c_int(klen3)), k)


def run_plot():
    nx = 1
    ny = 1
    nz = 1
    ne = 1
    head ="./"
    #outfile = head + "hist.hdf5"
    outfile = "./Output4D_00_1560.hdf5"
    f = h5py.File(outfile)
    data4 = f["data4"]
    print(data4.shape)
    condition = f["condition"]
    xi = 60/nx
    xe = 84/nx
    yi = 37/ny
    ye = 68/ny
    zi = 8/nz
    ze = 27/nz
    ei = 5/ne
    ee = 35/ne
    plot_crosssection(xi, xe, yi, ye, zi, ze, data4)
    #plot_crosssection(xi, xe, yi, ye, zi, ze, condition)


#run()
run_plot()
plt.show()
#plt.savefig("crosssection_of_histgram.png")


