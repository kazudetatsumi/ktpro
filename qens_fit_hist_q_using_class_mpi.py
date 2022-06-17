#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append("/home/kazu/ktpro")
from qens_fit_class import qens_fit as qf
from get_qlist_nova_class import get_qlist as gq
from mpi4py import MPI


def eachdata(sfile, qmin, qmax, fig):
    prj = gq(pklfile=sfile)
    head = sfile.split(".")[0]
    spectrafile = head + "_" + "{:.1f}".format(qmin) + \
        "-" + "{:.1f}".format(qmax) + ".pkl"
    prj.read_pkl()
    prj.spect(qmin, qmax)
    prj.save_spectra(spectrafile)
    plt.plot(prj.ene, prj.spectra)
    plt.yscale('log')
    #plt.show()
    return(spectrafile)


def fit(devf, tf, variables, elim, qmin, qmax):
    prj = qf(devf, tf, elim, showplot=False, leastsq=False)
    prj.icorr()
    prj.preprocessh(doicorr=True)
    prj.bg = 0.
    fig = plt.figure()
    prj.optimize(variables=variables,
                 figname='hist_fit_' + "{:.1f}".format(qmin) + '-'
                 + "{:.1f}".format(qmax) + '.png')


def run():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    prefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000001/"
    devallqf = prefix + "run6204s.pkl"
    tallqf = prefix + "run6202s.pkl"
    qmin = 0.55
    qmax = 0.70
    elim = [-0.03, 0.07]
    variables = [0.5, 0.02, 0.2, 0.01, 0.1, 0.]
    # variables = [0.5, 0.02, 0.01, 0.01]
    qmins = np.arange(0.2, 1.3, 0.1)
    qmaxs = qmins + 0.1
    qmin = qmins[rank]
    qmax = qmaxs[rank]
    fig = plt.figure()
    devf = eachdata(devallqf, qmin, qmax, fig)
    tf = eachdata(tallqf, qmin, qmax, fig)
    plt.savefig('hist_run_' + "{:.1f}".format(qmin) + '-' +
                "{:.1f}".format(qmax) + '.png')
    for qidx in range(0, size):
        comm.Barrier()
        if qidx == rank:
            print("")
            print("Now fitting the spectra in the q range")
            print('between ' + "{:.1f}".format(qmin) + ' and '
                  + "{:.1f}".format(qmax) + " Angs-1")
            fit(devf, tf, variables, elim, qmin, qmax)
        comm.Barrier()


run()
