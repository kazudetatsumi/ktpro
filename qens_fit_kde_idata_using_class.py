#!/usr/bin/env python
import numpy as np
import os
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class_kde_idata import runkdeidata as rki


def testrun():
    np.set_printoptions(linewidth=120)
    outfile = "./outkdeidata.pkl"
    alpha = 0.5
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    prefix = "/home/kazu/desktop/210108/Tatsumi/pickles/0000001io/"
    idevf = prefix + "qens_kde_results_on_idata_6204.pkl"
    itf = prefix + "qens_kde_results_on_idata_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 3
    binwidth1 = 0.0014
    binwidth2 = 0.0003
    if os.path.isfile(outfile):
        proj = rki(devf, tf, idevf, itf, outfile, alpha, elim, elimw,
                   leastsq=True, numcycle=numcycle)
        if proj.rank == 0:
            proj.loadfile()
            proj.output()
            proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = rki(devf, tf, idevf, itf, outfile, alpha, elim, elimw,
                   leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        if proj.rank == 0:
            proj.output()
            proj.savefile()


testrun()
