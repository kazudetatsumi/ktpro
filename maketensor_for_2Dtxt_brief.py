#!/usr/bin/env python
import numpy as np
import pickle
import sys
import os


def get2ddata(txtfile):
    #   genformtxt is too slow and requires too large memory for the purpose of
    #   extracting only one column. So, I comment out and use a for-loop without readlines().
    # data = np.genfromtxt(f, dtype=float, comments='#', delimiter=',')
    # intensity = data[:, 4]  # data[:, 4] are mask if the values are 1e*30
    intensity = []
    error = []
    energy = []
    q = []
    for il, line in enumerate(open(txtfile)):
        if "## XRANGE" in line:
            qbins = line[:-1].split("=")[1].split(",")
            Nq = int((float(qbins[1]) - float(qbins[0]))/float(qbins[2])) + 1
        if "## YRANGE" in line:
            ebins = line[:-1].split("=")[1].split(",")
            Ne = int((float(ebins[1]) - float(ebins[0]))/float(ebins[2])) + 1
        if il >= 12:
            intensity.append(float(line[:-1].split(',')[-2]))
            energy.append(float(line[:-1].split(',')[1]))
            q.append(float(line[:-1].split(',')[0]))
            _e = line[:-1].split(',')[-1]
            if _e == " 0":
                error.append(100000.)
            else:
                error.append(float(_e))
    intensity = np.array(intensity).reshape((Ne, Nq))
    error = np.array(error).reshape((Ne, Nq))
    energy = np.array(energy).reshape((Ne, Nq))
    q = np.array(q).reshape((Ne, Nq))
    dat = np.zeros((1, 3, Nq, Ne))
    dat[0, 0, :, :] = energy.T
    dat[0, 1, :, :] = intensity.T
    dat[0, 2, :, :] = error.T
    return dat


def gen_pkl(runno, head):
    txtfile = head + "runno" + runno + ".txt"
    if os.path.exists(txtfile):
       outfile = head + "run" + runno + "spectraorgtxt.pkl"
       dataset = get2ddata(txtfile)
       with open(outfile, 'wb') as f:
           pickle.dump(dataset, f, -1)


def run():
    if len(sys.argv) >= 2:
        runno = sys.argv[1]
    else:
        runno = "4174"
    head = "./"
    gen_pkl(runno, head)


run()
