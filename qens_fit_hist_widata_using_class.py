#!/usr/bin/env python
import os
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class_hist_widata import runhistwithidata as rhwid


def testrun():
    np.set_printoptions(linewidth=120)
    outfile = "./outhistwidata.pkl"
    alpha = 0.5
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    prefix = "/home/kazu/desktop/210108/Tatsumi/pickles/0000001io/"
    idevf = prefix + "qens_kde_results_on_idata_6204.pkl"
    itf = prefix + "qens_kde_results_on_idata_6202.pkl"
    idfrw = prefix + "../../srlz/0000001io/run6204united_monispectra.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 30000
    binwidth1 = 0.0016
    binwidth2 = 0.0016
    if os.path.isfile(outfile):
        proj = rhwid(devf, tf, idevf, itf, outfile, idfrw, alpha, elim,
                     elimw, leastsq=True, numcycle=numcycle)
        proj.loadfile()
        proj.output()
        proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = rhwid(devf, tf, idevf, itf, outfile, idfrw, alpha, elim,
                     elimw, leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        proj.output()
        proj.savefile()


testrun()
