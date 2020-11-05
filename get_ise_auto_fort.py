#!/usr/bin/env python
# This script calculates an integrated squared error between pdf and histogram. 
# python + fortran version of get_ise_auto.py.
# This version is much faster than get_ise_auto.py or get_ise_auto_mpi.py
# fortran source file of the library is isefort.f90
# Kazuyoshi TATSUMI 2020.
import numpy as np
import h5py
from ctypes import *
import sys

xi = 120
xe = 172
yi = 61
ye = 145
zi = 16
ze = 53
ei = 20
ee = 70


def calc_ise_f90(data, condition, pdfs, nw, sumpdf):
    lib = cdll.LoadLibrary("/home/kazu/ktpro/isefort4d.so")
    class result(Structure):
        _fields_ = [("isepon", c_double), ("avepdfs", c_double), ("fracise", c_double)]

    lib.ise4d.restype = result
    lib.ise4d.argtypes = [
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
                           np.ctypeslib.ndpointer(dtype=np.int32, ndim=4)
                           ]

    Nmax0 = data.shape[0]
    Nmax1 = data.shape[1]
    Nmax2 = data.shape[2]
    Nmax3 = data.shape[3]

    result = lib.ise4d(
            c_int(nw[0]),
            c_int(nw[1]),
            c_int(nw[2]),
            c_int(nw[3]),
            c_int(Nmax0),
            c_int(Nmax1),
            c_int(Nmax2),
            c_int(Nmax3),
            c_double(sumpdf),
            data,
            pdfs,
            condition
            )

    return  result.isepon, result.avepdfs, result.fracise


def read_h5py(outfile, term):
    f = h5py.File(outfile, 'r')
    return np.array(f[term])


def search_ns_in_resultfile(infile, tcount):
    f = open(infile, 'r')
    extraFound = 0
    line = f.readline()
    nsfound = 0
    while extraFound == 0:
        line = f.readline()
        if line.startswith('extrapolation'):
            extraFound = 1
        else:
            values = line.split()
            if abs(tcount - float(values[0]))/tcount <= 0.01:
                nsfound = 1
                return(np.array(values[1:5], dtype=int))
    if nsfound == 0:
        print("ns is not found!!")
        sys.exit()


def preconditioning(workdir):
    data4 = read_h5py(workdir+"/eliminated_data.hdf5", "data4")
    condition = read_h5py(workdir+"/eliminated_data.hdf5", "condition")
    condition = np.array(condition, dtype="int32")
    tcount = np.sum(data4)
    nw = search_ns_in_resultfile(workdir+"../result.txt_vec", tcount)
    pdf = read_h5py("pdf.hdf5", "pdf")
    sumpdf = np.sum(pdf)
    pdfs = np.array(pdf[xi:xe, yi:ye, zi:ze, ei:ee]) # this np.array treatment is required for correct data transfer to fortran libraly
    ise, avepdfs, fracise = calc_ise_f90(data4, condition, pdfs, nw, sumpdf)
    return np.sum(data4), ise, avepdfs, fracise


def run():
    dirlist=["no0.5","no1","no1.5","no2","no2.5","no3","no3.5", "no4","no5", "no4.5","no10", "no24", "no73", "no120", "no124", "no211"]
    #dirlist=["no1"]
    for dirname in dirlist:
        workdir="/home/kazu/desktop/200312/for_cu_new/orthotope_opt_20-70meV_ddscs_again2/"+dirname+"/"
        nt, ise, avepdfs, fracise = preconditioning(workdir)
        print(nt, ise, avepdfs, fracise)

run()
