#!/usr/bin/env python
# This script two-lorentzian-fits kde results, and reconstructs the target
# function with the fitting parameter.
# The reconstructed ML function will be used to generate simulated data.
# example: qens_fit_kde_vs_hist_2lore_reconst_using_class.py 6202 6204 0875
# Kazuyoshi TATSUMI 2022/04/16
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
import qens_fit_class as qfc


def sqrun_kde_hist_2lore():
    np.set_printoptions(linewidth=120)
    if len(sys.argv) >= 2:
        runno = sys.argv[1]
    else:
        print("using default runno 6202")
        runno = "6202"
    if len(sys.argv) >= 3:
        runnod = sys.argv[2]
    else:
        print("using default runnod 6204")
        runnod = "6204"
    #if len(sys.argv) >= 4:
    #    frac = sys.argv[3]
    #else:
    #    print("using default frac """)
    #    frac = ""
    frac = ""

    head = "./"
    devf = head + "qens_sim_kde_"+runnod+".pkl"
    tf = head + "qens_sim_kde_"+runno+frac+".pkl"
    elim = [-0.03, 0.07]
    proj = qfc.qens_fit(devf, tf, elim, showplot=False)
    proj.preprocessnoi(doicorr=False)
    proj.optimize(variables=[1.58704786e-04, 2.67980295e-02, 6.02405238e-05,
                             6.88866588e-03, 2.21127501e-01, 5.02759930e-02],
                  figname="qens_kde_fit2.png")


sqrun_kde_hist_2lore()
