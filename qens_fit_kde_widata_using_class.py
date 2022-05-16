#!/usr/bin/env python
import numpy as np
import os
import sys
from mpi4py import MPI
sys.path.append("/home/kazu/ktpro")
from qens_fit_class_kde_widata import runkdewithidata as rkw


def testrun():
    np.set_printoptions(linewidth=120)
    prefix = "/home/kazu/desktop/210108/Tatsumi/winparam_exam/" +\
             "test_mpi_fort_de_0.000025/160_1_0000001io_Boxcar_simu/runtst/tmp/"
    #outfile = prefix + "./outkdewidata_least_squares_alpha1_k3.pkl"
    outfile = prefix + "./outkdewidata_least_squares_alpha1_k3_test.pkl"
    alpha = 1.0
    devf = prefix + "./qens_kde_o_divided_by_i_6204.pkl"
    tf = prefix + "./qens_kde_o_divided_by_i_6202.pkl"
    prefix = "/home/kazu/desktop/210108/Tatsumi/pickles/0000001io/"
    idevf = prefix + "qens_kde_results_on_idata_6204.pkl"
    itf = prefix + "qens_kde_results_on_idata_6202.pkl"
    idfrw = prefix + "../../srlz/0000001io/run6204united_monispectra.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 10
    binwidth1 = 0.0010
    binwidth2 = 0.0005
    if os.path.isfile(outfile):
        proj = rkw(devf, tf, idevf, itf, outfile, idfrw, alpha, elim, elimw,
                   leastsq=False, numcycle=numcycle)
        if proj.rank == 0:
            proj.loadfile()
            proj.output()
            proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = rkw(devf, tf, idevf, itf, outfile, idfrw, alpha, elim, elimw,
                   leastsq=False, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        if proj.rank == 0:
            proj.output()
            print("WinFunc:", proj.WinFunc)
            print("M:", proj.M)
            print("winparam:", proj.winparam)
            proj.savefile()


testrun()
