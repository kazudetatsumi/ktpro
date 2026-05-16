#!/usr/bin/env python
# This is a form of class generated based on optbinwidth4d_withcond_mpi.py
# I wrote this because I need to use the same INS intensities but change the
# mask i.e., an array of 'condition' many times and separate the input file
# into two files for the INS intensities and the mask respectively, for saving
# the amount of hdd.
# Kazuyoshi TATSUMI 2021/01/27
import numpy as np
import h5py
import os
from ctypes import *
from mpi4py import MPI
lib = CDLL("/home/kazu/ktpro/costfort4d_mpi.so")


class Opt_Bin_Width:

    def __init__(self, datafile, condfile, usecond, savecn=True):
        self.datafile = datafile
        self.condfile = condfile
        self.usecond = usecond
        self.savecn = savecn

    def calc_cost4d_f90(self, data, condition, maxw):
        print(maxw)
        lib.cost4d.restype = c_void_p
        lib.cost4d.argtypes = [
                               POINTER(c_int32),
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

        comm = MPI.COMM_WORLD
        comm = comm.py2f()
        rank = MPI.COMM_WORLD.Get_rank()
        size = MPI.COMM_WORLD.Get_size()
        lib.cost4d(
                   c_int32(comm),
                   c_int(maxw[0]),
                   c_int(maxw[1]),
                   c_int(maxw[2]),
                   c_int(maxw[3]),
                   c_int(Nmax0),
                   c_int(Nmax1),
                   c_int(Nmax2),
                   c_int(Nmax3),
                   c_bool(self.usecond),
                   data,
                   condition,
                   Cn,
                   kaves,
                   deltas
                   )
        print(rank)
        print("---check the calculated cost function---")
        opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
        opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1,
                    opt_indx[3] + 1)
        print("opt_indx for Cn", opt_indx)
        print("---save results in Cn.hdf5---")
        outfile = "Cn_rank"+str(rank)+".hdf5"
        with h5py.File(outfile, 'w') as hf:
            hf.create_dataset('Cn', data=Cn)
            hf.create_dataset('kave', data=kaves)
            hf.create_dataset('delta', data=deltas)

        MPI.COMM_WORLD.barrier()
        print("size", size)
        if rank == 0:
            print("rank", rank)
            for i in range(1, size):
                f = h5py.File("Cn_rank"+str(i)+".hdf5", 'r')
                Cn_tmp = f["Cn"]
                Cn = Cn + Cn_tmp
                kaves_tmp = f["kave"]
                kaves = kaves + kaves_tmp
                deltas_tmp = f["delta"]
                deltas = deltas + deltas_tmp
            Cn[0, 0, 0, 0] = Cn_tmp[0, 0, 0, 0]
            kaves[0, 0, 0, 0] = kaves_tmp[0, 0, 0, 0]
            deltas[0, 0, 0, 0] = deltas_tmp[0, 0, 0, 0]
            opt_indx = np.unravel_index(np.argmin(Cn, axis=None), Cn.shape)
            opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1,
                        opt_indx[3] + 1)
            print("opt_indx for Cn", opt_indx, "with the Cn =",
                  Cn[opt_indx[0] - 1, opt_indx[1] - 1, opt_indx[2] - 1,
                      opt_indx[3] - 1])
            if self.savecn:
                with h5py.File("Cn.hdf5", 'w') as hf:
                    hf.create_dataset('Cn', data=Cn)
                    hf.create_dataset('kave', data=kaves)
                    hf.create_dataset('delta', data=deltas)
        MPI.COMM_WORLD.barrier()
        os.remove(outfile)

    def read_h5py(self, infile, prop):
        f = h5py.File(infile, 'r')
        return(f[prop])

    def preprocess(self):
        data = np.array(self.read_h5py(self.datafile, 'data4'))*1.0
        condition = np.array(self.read_h5py(self.condfile, 'condition'),
                             dtype=np.int32)
        print("size of data is", data.shape)
        print("usecond:", self.usecond)
        n = np.sum(data)*1.0
        print("n=", n)

        maxxwidth = int(data.shape[0] // 2)
        maxywidth = int(data.shape[1] // 2)
        maxzwidth = int(data.shape[2] // 2)
        maxowidth = int(data.shape[3] // 2)
        maxowidth = int(data.shape[3] // 2)
        print("maxwidth:", maxxwidth, maxywidth, maxzwidth, maxowidth)
        maxw = np.array([maxxwidth, maxywidth, maxzwidth, maxowidth])
        return data, condition, maxw
