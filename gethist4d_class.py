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
import h5py
import os
import ctypes
lib = ctypes.CDLL("/home/kazu/ktpro/histfort4d.so")


class gethist4d_class:

    def __init__(self, nw):
        self.nw = nw

    def calc_hist4d_f90(self):
        class result(ctypes.Structure):
            _fields_ = [("len0", ctypes.c_int), ("len1", ctypes.c_int),
                        ("len2", ctypes.c_int),  ("len3", ctypes.c_int),
                        ("arr", ctypes.POINTER(ctypes.c_double))]

        lib.hist4d.restype = result
        lib.hist4d.argtypes = [ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)]

        Nmax0 = self.A.shape[0]
        Nmax1 = self.A.shape[1]
        Nmax2 = self.A.shape[2]
        Nmax3 = self.A.shape[3]

        result = lib.hist4d(ctypes.byref(ctypes.c_int(self.nw[0])),
                            ctypes.byref(ctypes.c_int(self.nw[1])),
                            ctypes.byref(ctypes.c_int(self.nw[2])),
                            ctypes.byref(ctypes.c_int(self.nw[3])),
                            ctypes.byref(ctypes.c_int(Nmax0)),
                            ctypes.byref(ctypes.c_int(Nmax1)),
                            ctypes.byref(ctypes.c_int(Nmax2)),
                            ctypes.byref(ctypes.c_int(Nmax3)), self.A)
        result_len0 = result.len0
        result_len1 = result.len1
        result_len2 = result.len2
        result_len3 = result.len3
        self.k = np.ctypeslib.as_array(result.arr,
                                       shape=(result_len0, result_len1,
                                              result_len2, result_len3, ))

        self.kcond = np.zeros((result_len0, result_len1, result_len2,
                               result_len3))
        for i in range(0, result_len0):
            for j in range(0, result_len1):
                for h in range(0, result_len2):
                    for l in range(0, result_len3):
                        self.kcond[i, j, h, l] = \
                          np.sum(self.condition[i*self.nw[0]:(i+1)*self.nw[0],
                                                j*self.nw[1]:(j+1)*self.nw[1],
                                                h*self.nw[2]:(h+1)*self.nw[2],
                                                l*self.nw[3]:(l+1)*self.nw[3]])

    def save_hist(self):
        with h5py.File(self.outfile, 'w') as hf:
            hf.create_dataset('data4', data=self.k)
            hf.create_dataset('condition', data=self.kcond)

    def histprocess(self, head=",/"):
        self.datafile = head + "eliminated_data.hdf5"
        self.outfile = head + "hist_eliminated.hdf5"
        if not os.path.isfile(self.outfile):
            f = h5py.File(self.datafile, 'r')
            self.data = f["data4"][:]*1.0 # nqx, nqy, nqz, nomega
            self.condition = f["condition"][:]
            self.A = np.cumsum(np.cumsum(np.cumsum(np.cumsum(self.data, axis=0),
                               axis=1), axis=2), axis=3)
            self.calc_hist4d_f90()
            self.save_hist()
        #else:
        #   print(self.outfile + " already exists, skipping histogramming")


def samplerun():
    nws = np.array([2, 4, 2, 5])
    project = gethist4d_class(nws)
    project.histprocess()


#samplerun()
