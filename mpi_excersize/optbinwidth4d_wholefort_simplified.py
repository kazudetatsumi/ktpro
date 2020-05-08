#!/usr/bin/env python
import numpy as np
import h5py
from ctypes import *
lib = CDLL("./costfort4d.so")


def calc_cost4d_f90(maxw, data, condition, usecond):
    lib.cost4d.restype = c_void_p
    lib.cost4d.argtypes = [
                           POINTER(c_int),
                           POINTER(c_int),
                           POINTER(c_int),
                           POINTER(c_int),
                           POINTER(c_int),
                           POINTER(c_int),
                           POINTER(c_int),
                           POINTER(c_int),
                           POINTER(c_bool),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.int32, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           ]

    Nmax0 = data.shape[0]
    Nmax1 = data.shape[1]
    Nmax2 = data.shape[2]
    Nmax3 = data.shape[3]
    Cn = np.zeros((maxw[0], maxw[1], maxw[2], maxw[3]), dtype=np.float64)
    kaves = np.zeros((maxw[0], maxw[1], maxw[2], maxw[3]), dtype=np.float64)
    deltas = np.zeros((maxw[0], maxw[1], maxw[2], maxw[3]), dtype=np.float64)

    lib.cost4d(
               c_int(maxw[0]),
               c_int(maxw[1]),
               c_int(maxw[2]),
               c_int(maxw[3]),
               c_int(Nmax0),
               c_int(Nmax1),
               c_int(Nmax2),
               c_int(Nmax3),
               c_bool(usecond),
               data,
               condition,
               Cn,
               kaves,
               deltas
               )

    print("---check the calculated cost function---")
    opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
    opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1, opt_indx[3] + 1)
    print("opt_indx for Cn", opt_indx)
    print("---save results in Cn.hdf5---")
    outfile = "Cn.hdf5"
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('Cn', data=Cn)
        hf.create_dataset('kave', data=kaves)
        hf.create_dataset('delta', data=deltas)


def run():
    datafile = "/home/kazu/desktop/200312/for_cu/with_cond/orthotope_opt/16h/eliminated_data.hdf5"
    f = h5py.File(datafile)
    data = f["data4"][:, :, :, :]  # nqx, nqy, nqz, nomega
    print("size of data is", data.shape)
    condition = np.array(f["condition"], dtype=np.int32)
    usecond = False
    print("usecond:", usecond)
    
    n = np.sum(data)*1.0
    print("n=", n)

    maxxwidth = int(data.shape[0] // 2)
    maxywidth = int(data.shape[1] // 2)
    maxzwidth = int(data.shape[2] // 2)
    maxowidth = int(data.shape[3] // 2)
    print("maxwidth:", maxxwidth, maxywidth, maxzwidth, maxowidth)
    
    maxw = np.array([maxxwidth, maxywidth, maxzwidth, maxowidth])
    calc_cost4d_f90(maxw, data, condition, usecond)


run()

