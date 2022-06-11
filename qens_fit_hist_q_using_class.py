#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append("/home/kazu/ktpro")
from qens_fit_class import qens_fit as qf
from get_qlist_nova_class import get_qlist as gq


def eachdata(sfile, qmin, qmax):
    prj = gq(pklfile=sfile)
    head = sfile.split(".")[0]
    spectrafile = head + "_" + str(qmin) + "-" + str(qmax) + ".pkl"
    prj.read_pkl()
    prj.spect(qmin, qmax)
    prj.save_spectra(spectrafile)
    #plt.plot(prj.ene, prj.spectra)
    #plt.show()
    return(spectrafile)


def fit(devf, tf, variables, elim):
    prj = qf(devf, tf, elim, showplot=True, leastsq=False)
    prj.icorr()
    prj.preprocessh(doicorr=True)
    prj.bg = 0.
    prj.optimize(variables=variables, figname='test.png')


def run():
    prefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000025/"
    devallqf = prefix + "run6204s.pkl"
    tallqf = prefix + "run6202s.pkl"
    qmin = 0.55
    qmax = 0.70
    elim = [-0.03, 0.07]
    variables = [0.5, 0.02, 0.2, 0.01, 0.1, 0.]
    #variables = [0.5, 0.02, 0.01, 0.01]
    qmins = np.arange(0.2, 0.9, 0.1)
    qmaxs = qmins + 0.1
    for qmin, qmax in zip(qmins, qmaxs):
        print(qmin, qmax)
        devf = eachdata(devallqf, qmin, qmax)
        tf = eachdata(tallqf, qmin, qmax)
        fit(devf, tf, variables, elim)


run()
