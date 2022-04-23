#!/usr/bin/env python
import os
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class_hist_noidata import runhistnoidata as rhn


def testrun():
    prefix = "/home/kazu/desktop/210108/Tatsumi/winparam_exam/" +\
             "test_mpi_fort_de_0.000025/160_1_0000001io_Boxcar_simu/runtst/tmp"
    outfile = prefix + "/outhistnoidata.pkl"
    outfile = "./outhistnoidata.pkl"
    alpha = 0.5
    devf = prefix + "/qens_kde_o_divided_by_i_6204.pkl"
    tf = prefix + "/qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 30
    binwidth1 = 0.0014
    binwdith2 = 0.00074
    binwdith2 = 0.0003
    np.set_printoptions(linewidth=120)
    if os.path.isfile(outfile):
        proj = rhn(devf, tf, outfile, alpha, elim, elimw,
                   leastsq=True, numcycle=numcycle)
        proj.loadfile()
        proj.output()
        proj.plot_distribution(binwidth1, binwdith2)
    else:
        np.random.seed(314)
        proj = rhn(devf, tf, outfile, alpha, elim, elimw,
                   leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        proj.output()
        proj.savefile()


testrun()
