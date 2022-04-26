#!/usr/bin/env python
import numpy as np
import os
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class_kde import runkdenoidata as rkn


def testrun():
    prefix = "/home/kazu/desktop/210108/Tatsumi/winparam_exam/" +\
             "test_mpi_fort_de_0.000025/160_1_0000001io_Boxcar_simu/runtst/tmp"
    outfile = prefix + "/outkde.pkl"
    #outfile = "./outkde.pkl"
    alpha = 0.5
    devf = prefix + "/qens_kde_o_divided_by_i_6204.pkl"
    tf = prefix + "/qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 10
    binwidth1 = 0.0014
    binwidth2 = 0.00074
    binwidth2 = 0.0003
    if os.path.isfile(outfile):
        proj = rkn(devf, tf, outfile, alpha, elim, elimw,
                   leastsq=True, numcycle=numcycle)
        if proj.rank == 0:
            proj.loadfile()
            proj.output()
            proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = rkn(devf, tf, outfile, alpha, elim, elimw,
                   leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        if proj.rank == 0:
            proj.output()
            proj.savefile()


testrun()
