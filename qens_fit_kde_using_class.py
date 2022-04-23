#!/usr/bin/env python
import numpy as np
import os
import sys
import pickle
import matplotlib.pyplot as plt
from ctypes import *
from mpi4py import MPI
sys.path.append("/home/kazu/ktpro")
from qens_class_fort_mpi import qens as qc
from qens_fit_class_hist import runhist as rh
lib = CDLL("/home/kazu/ktpro/ssvkernel_f90_mpi.so")
libssk = CDLL("/home/kazu/ktpro/sskernel_f90.so")
from qens_fit_class_kde import runkdenoidata as rkn


def testrun():
    prefix = "/home/kazu/desktop/210108/Tatsumi/winparam_exam/" +\
             "test_mpi_fort_de_0.000025/160_1_0000001io_Boxcar_simu/runtst/tmp"
    outfile = prefix + "/outkde.pkl"
    alpha = 0.5 
    devf = prefix + "/qens_kde_o_divided_by_i_6204.pkl"
    tf = prefix + "/qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 30000
    binwidth1 = 0.0016
    binwdith2 = 0.00074
    binwdith2 = 0.0016
    if os.path.isfile(outfile):
        proj = rkn(devf, tf, outfile, alpha, elim, elimw,
                   leastsq=True, numcycle=numcycle)
        proj.loadfile()
        proj.output()
        proj.plot_distribution(binwidth1, binwdith2)
    else:
        np.random.seed(314)
        proj = rkn(devf, tf, outfile, alpha, elim, elimw,
                   leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        proj.output()
        proj.savefile()


testrun()
