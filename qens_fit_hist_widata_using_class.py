#!/usr/bin/env python
import os
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class_hist_widata import runhistwithidata as rhwid


def testrun():
    np.set_printoptions(linewidth=120)
    prefix = "/home/kazu/desktop/210108/Tatsumi/winparam_exam/" +\
             "test_mpi_fort_de_0.000025/" +\
             "160_1_0000001io_Boxcar_simu/runtst/tmp/"
    outfile = prefix + "./outhistwidatanonneg.pkl"
    alpha = 0.5
    devf = prefix + "./qens_kde_o_divided_by_i_6204.pkl"
    tf = prefix + "./qens_kde_o_divided_by_i_6202.pkl"
    prefix = "/home/kazu/desktop/210108/Tatsumi/pickles/0000001io/"
    idevf = prefix + "qens_kde_results_on_idata_6204.pkl"
    itf = prefix + "qens_kde_results_on_idata_6202.pkl"
    idfrw = prefix + "../../srlz/0000001io/run6204united_monispectra.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 30000
    binwidth1 = 0.0010
    binwidth2 = 0.0004
    if os.path.isfile(outfile):
        proj = rhwid(devf, tf, idevf, itf, outfile, idfrw, alpha, elim,
                     elimw, leastsq=False, numcycle=numcycle)
        proj.loadfile()
        proj.output()
        proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = rhwid(devf, tf, idevf, itf, outfile, idfrw, alpha, elim,
                     elimw, leastsq=False, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        proj.output()
        proj.savefile()


testrun()
