#!/usr/bin/env python
# This script calculates an integrated squared error between pdf and histogram. 
# You should list the work dirs of optimal bin-widths calculation  and this script calculates the histogram with the optimal bin-widths.
# You should prepare pdf.hdf5.
# Kazuyoshi TATSUMI 2020.
import numpy as np
import h5py
import ctypes
import sys
lib = ctypes.CDLL("/home/kazu/ktpro/histfort4d.so")

xi = 120
xe = 172
yi = 61
ye = 145
zi = 16
ze = 53
ei = 20
ee = 70


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

    kcond = np.zeros((result_len0, result_len1, result_len2, result_len3))
    for i in range(0, result_len0):
        for j in range(0, result_len1):
            for h in range(0, result_len2):
                for l in range(0, result_len3):
                    kcond[i, j, h, l] = np.sum(condition[i*nw0:(i+1)*nw0, j*nw1:(j+1)*nw1, h*nw2:(h+1)*nw2, l*nw3:(l+1)*nw3])
    return result_vec, result_len0, result_len1, result_len2, result_len3, kcond


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


def get_hist(workdir):
    data4 = read_h5py(workdir+"/eliminated_data.hdf5", "data4")
    condition = read_h5py(workdir+"/eliminated_data.hdf5", "condition")
    tcount = np.sum(data4)
    ns = search_ns_in_resultfile(workdir+"../result.txt_vec", tcount)
    A = np.cumsum(np.cumsum(np.cumsum(np.cumsum(data4, axis=0), axis=1), axis=2), axis=3)
    datahist, klen0, klen1, klen2, klen3, condhist = calc_hist4d_f90(A, data4, ns[0], ns[1], ns[2], ns[3], condition)
    return(datahist, condhist, ns)


def get_ise(workdir):
    # get histogram data
    datahist, condhist, ns = get_hist(workdir)
    #print("histogram is obtained")

    # get pdf data
    pdf = read_h5py("pdf.hdf5", "pdf")
    pdfs = pdf[xi:xe, yi:ye, zi:ze, ei:ee]
    #print("pdf is obtained")

    histsize = datahist.shape
    ise = 0.
    sumhist = 0.
    sumpdfs = 0.
    sumpdf = np.sum(pdf)
    sumhist = np.sum(condhist[condhist == np.max(condhist)]/np.max(condhist)*1.0)
    sumdata = np.sum(datahist[condhist == np.max(condhist)])
    for hx in range(0, histsize[0]):
        for hy in range(0, histsize[1]):
            for hz in range(0, histsize[2]):
                for hw in range(0, histsize[3]):
                    if condhist[hx, hy, hz, hw] == np.max(condhist):
                        sumpdfs += np.sum(pdfs[hx*ns[0]:(hx+1)*ns[0], hy*ns[1]:(hy+1)*ns[1], hz*ns[2]:(hz+1)*ns[2], hw*ns[3]:(hw+1)*ns[3]])
    for hx in range(0, histsize[0]):
        for hy in range(0, histsize[1]):
            for hz in range(0, histsize[2]):
                for hw in range(0, histsize[3]):
                    if condhist[hx, hy, hz, hw] == np.max(condhist):
                        for px in range(hx*ns[0], (hx+1)*ns[0]):
                            for py in range(hy*ns[1], (hy+1)*ns[1]):
                                for pz in range(hz*ns[2], (hz+1)*ns[2]):
                                    for pw in range(hw*ns[3], (hw+1)*ns[3]):
                                        ise += (pdfs[px, py, pz, pw] - datahist[hx, hy, hz, hw]/(1.0*sumdata/sumpdfs*sumpdf*np.prod(ns)))**2
    ise = ise / (np.prod(ns)*sumhist*1.0)
    #print("average_pdfs:",sumpdfs/((np.prod(ns)*sumhist*1.0)))
    #print("ise",ise)
    #print("sq_ise/average_pdf (%)", ise**0.5/(sumpdfs/((np.prod(ns)*sumhist*1.0)))*100)
    print(np.sum(datahist),ise)


def run():
    dirlist=["no0.5","no1","no1.5","no2","no2.5","no3","no3.5", "no4","no5", "no4.5","no10", "no24", "no73", "no120", "no124", "no211"]
    #dirlist=["no1"]
    for dirname in dirlist:
        workdir="/home/kazu/desktop/200312/for_cu_new/orthotope_opt_20-70meV_ddscs_again2/"+dirname+"/"
        get_ise(workdir)

run()
