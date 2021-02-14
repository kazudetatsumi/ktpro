#!/usr/bin/env python
# This script calculates an integrated squared error between pdf and histogram
# for each of binindx sets whose binidx le maxw.
# Kazuyoshi TATSUMI 2020.
import numpy as np
import h5py
import subprocess
import os
from ctypes import *
from mpi4py import MPI
lib = CDLL("/home/kazu/ktpro/isefortsearch_mpi.so")


def calc_ise_search_f90(data, condition, pdfs, maxw, sumpdf):

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
                           POINTER(c_int),
                           POINTER(c_double),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.int32, ndim=4),
                           np.ctypeslib.ndpointer(dtype=np.float64, ndim=4)
                           ]

    Nmax0 = data.shape[0]
    Nmax1 = data.shape[1]
    Nmax2 = data.shape[2]
    Nmax3 = data.shape[3]
    ise = np.zeros((maxw[0], maxw[1], maxw[2], maxw[3]), dtype=np.float64)

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
             c_double(sumpdf),
             data,
             pdfs,
             condition,
             ise
            )

    outfile = "ise_rank"+str(rank)+".hdf5"
    with h5py.File(outfile, 'w') as hf: 
        hf.create_dataset('ise', data=ise)

    MPI.COMM_WORLD.barrier()
    print("size", size)
    if rank == 0:
        print("rank", rank)
        for i in range(1, size):
            f = h5py.File("ise_rank"+str(i)+".hdf5", 'r')
            ise_tmp = f["ise"]
            ise = ise + ise_tmp
        opt_indx = np.unravel_index(np.argmin(ise, axis=None), ise.shape)
        opt_indx = (opt_indx[0] + 1, opt_indx[1] + 1, opt_indx[2] + 1, opt_indx[3] + 1)
        print("opt_indx for Cn", opt_indx, "with the ise =", ise[opt_indx[0] - 1, opt_indx[1] - 1, opt_indx[2] - 1, opt_indx[3] - 1])
        with h5py.File("ise.hdf5", 'w') as hf: 
            hf.create_dataset('ise', data=ise)
    MPI.COMM_WORLD.barrier()
    os.remove(outfile)


def read_h5py(outfile, term):
    f = h5py.File(outfile, 'r')
    return np.array(f[term])


def getbinidx2(logfile):
    line = subprocess.getoutput("grep \"with the\" " + logfile)
    values = line.split("(")
    values2 = values[1].split(")")
    values3 = np.array(values2[0].split(","), dtype='int32')
    return values3


def preconditioning(rootdir):
    data4 = read_h5py(rootdir+"try1/26m/eliminated_data.hdf5", "data4")
    condition = read_h5py(rootdir+"try1/26m/eliminated_data.hdf5", "condition")
    condition = np.array(condition, dtype="int32")
    maxnw_data = [data4.shape[0]/2, data4.shape[1]/2, data4.shape[2]/2, data4.shape[3]/2]
    maxnw = getbinidx2(rootdir+"try1/26m/std-26m.log")*2
    maxnw = [min(maxnw[i], maxnw_data[i]) for i in range(0,4)]
    #maxnw = [2, 2, 2, 2]
    print(maxnw)
    pdf = read_h5py(rootdir + "pdf.hdf5", "pdf")
    sumpdf = np.sum(pdf)
    pdfs = np.array(pdf[:, :, :, 150:])
    calc_ise_search_f90(data4, condition, pdfs, maxnw, sumpdf)
    #ise, avepdfs, fracise = calc_ise_search_f90(data4, condition, pdfs, maxnw, sumpdf)
    #return np.sum(data4), ise, avepdfs, fracise


def run():
    rootdir = "/home/kazu/desktop/200903/morethan10meV/" +\
              "samedataasfilled_again2_ddscs/"
    preconditioning(rootdir)
    #nt, ise, avepdfs, fracise = preconditioning(rootdir)
    #print(nt, ise, avepdfs, fracise)


run()
