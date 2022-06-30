#!/usr/bin/env python
import numpy as np
import sys
import pickle
import matplotlib.pyplot as plt
sys.path.append("/home/kazu/ktpro")
from qens_fit_class import qens_fit as qf
from get_qlist_nova_class_mpi import get_qlist as gq
from mpi4py import MPI



def eachdata(sfile, qmin, qmax, fig):
    prj = gq(pklfile=sfile)
    head = sfile.split(".")[0]
    spectrafile = head + "_" + "{:.2f}".format(qmin) + \
        "-" + "{:.2f}".format(qmax) + ".pkl"
    prj.read_pkl()
    prj.spect(qmin, qmax)
    prj.save_spectra(spectrafile)
    plt.plot(prj.ene, prj.spectra)
    plt.yscale('log')
    #plt.show()
    return(spectrafile)


def eachdata_justfilename(sfile, qmin, qmax):
    head = sfile.split(".")[0]
    spectrafile = head + "_" + "{:.2f}".format(qmin) + \
        "-" + "{:.2f}".format(qmax) + ".pkl"
    return(spectrafile)


def fit(devf, tf, variables, elim, qmin, qmax):
    prj = qf(devf, tf, elim, showplot=False, leastsq=False)
    prj.icorr()
    prj.preprocessh(doicorr=True)
    prj.bg = 0.
    fig = plt.figure()
    prj.optimize(variables=variables, figname='hist_fit_' +
                 "{:.2f}".format(qmin) + '-' + "{:.2f}".format(qmax) + '.png')


def run():
    prefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000025/"
    devallqf = prefix + "run6204s.pkl"
    tallqf = prefix + "run6205s.pkl"
    qmin = 0.55
    qmax = 0.70
    elim = [-0.03, 0.07]
    #variables = [0.5, 0.02, 0.2, 0.01, 0.1, 0.]
    variables = [0.5, 0.02, 0.02, 0.01]
    #variables = [0.5, 0.02, 0.01, 0.01]
    qmins = np.arange(0.2, 1.3, 0.1)
    qmaxs = qmins + 0.1
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        print(devallqf, tallqf, qmin)
    for qmin, qmax in zip(qmins, qmaxs):
        #savefile = "./q" + "{:.2f}".format(qmin) + "-" + \
        #        "{:.2f}".format(qmax) + ".pkl"
        fig = plt.figure()
        devf = eachdata(devallqf, qmin, qmax, fig)
        tf = eachdata(tallqf, qmin, qmax, fig)
        if rank == 0:
            print(qmin, qmax)
            plt.savefig('hist_run_' + '{:.2f}'.format(qmin) + '-' +
                        '{:.2f}'.format(qmax) + '.png')
            fit(devf, tf, variables, elim, qmin, qmax)


def fitrun():
    prefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000025/"
    devallqf = prefix + "run6204s.pkl"
    tallqf = prefix + "run6205s.pkl"
    qmin = 0.55
    qmax = 0.70
    elim = [-0.03, 0.07]
    #variables = [0.8, 0.0001, 0.24, 0.0002, 0.001, 1.2]
    variables = [0.19, 0.013, 0.25, 0.12]
    qmins = np.arange(0.9, 1.3, 0.1)
    qmaxs = qmins + 0.1
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    for qmin, qmax in zip(qmins, qmaxs):
        fig = plt.figure()
        devf = eachdata_justfilename(devallqf, qmin, qmax)
        tf = eachdata_justfilename(tallqf, qmin, qmax)
        if rank == 0:
            print(qmin, qmax)
            fit(devf, tf, variables, elim, qmin, qmax)


run()
