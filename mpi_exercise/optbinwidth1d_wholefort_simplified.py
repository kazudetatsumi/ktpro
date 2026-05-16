#!/usr/bin/env python
import numpy as np
import h5py
from ctypes import *
lib = CDLL("./costfort1d.so")


def calc_cost1d_f90(A, maxw, condition):
    lib.cost1d.restype = c_void_p
    lib.cost1d.argtypes = [
                           POINTER(c_int),
                           POINTER(c_int),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)
                          ]
    lenA = c_int(A.shape[0])
    Cn = np.zeros(maxw, dtype=np.float64)
    kaves = np.zeros(maxw, dtype=np.float64)
    deltas = np.zeros(maxw, dtype=np.float64)
    lib.cost1d(c_int(maxw), lenA, A, Cn, kaves, deltas)
    return Cn, kaves, deltas


def run1d():
    datafile = "/home/kazu/WORK/vasp-phonopy/cscl/phonopy_222/m200200200/data3_100000000.hdf5"
    f = h5py.File(datafile)
    fflag = 1
    data = f["data3"][:] # nqx, nqy, nqz, nomega
    data = data[0, 0, 0, :]*1.0
    n = np.sum(data)*1.0
    print("n=", n)
    condition = np.ones(data.shape, dtype=bool)
    maxw =  int(data.shape[0] / 2)
    A = np.cumsum(data)

    Cn, kaves, deltas = calc_cost1d_f90(A, maxw, condition)
    Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)
    print("Cn:", Cn)
    print("argmin(Cn)",np.argmin(Cn)+1) # Because fortran array is 1-origin while python array is 0-origin.

    m = 0.1*n

    ex = (1/m - 1/n) * kaves / (deltas**2*n) 

    Cm = ex + Cn

    print("Cm with m = ", m, Cm)
    print("argmin(Cm)",np.argmin(Cm)+1) # Because fortran array is 1-origin while python array is 0-origin.


run1d()

