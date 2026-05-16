#!/usr/bin/env python
# This script do fitting analysis on kde results and then on histograms.
# For kde results, the background is optimized.
# For histograms, the background is fixed as the bkg_peak_ratio(kde) * peak(histogram)
# as suggested by Dr. Matsuura.
# example usage qens_fit_kde_vs_hist_using_class.py 6204 6202 0000025io 0875
# Kazuyoshi TATSUMI 2021/09/06
import numpy as np
import os
import sys
sys.path.append("/home/kazu/ktpro")
import qens_fit_class as qfc


def sqrun_kde_hist():
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
    if len(sys.argv) >= 4:
        dirname = sys.argv[3]
    else:
        print("using default dirname 000025io")
        dirname = "000025io"
    if len(sys.argv) >= 5:
        frac = sys.argv[4]
    else:
        print("using default frac """)
        frac = ""

    #head = "/home/kazu/desktop/210108/Tatsumi/pickles/"+dirname+"/"
    #note that we use the device function collected in the total measument time.
    head = "./"
    devf = head + "qens_kde_o_divided_by_i_"+runnod+".pkl"
    tf = head + "qens_kde_o_divided_by_i_"+runno+frac+".pkl"
    elim = [-0.03, 0.07]
    proj = qfc.qens_fit(devf, tf, elim, showplot=False, leastsq=False)
    proj.icorr()
    proj.preprocess(doicorr=True)
    #proj.optimize(figname="qens_kde_fit.png")
    proj.optimize(variables=[0.38705532, 0.01308485, 0.53720324, 2.19106007], figname=None)
    #head = "/home/kazu/desktop/210108/Tatsumi/srlz/"+dirname+"/"
    devf = head + "qens_hist_o_divided_by_i_"+runnod+".pkl"
    tf = head + "qens_hist_o_divided_by_i_"+runno+frac+".pkl"
    if os.path.isfile(devf) and os.path.isfile(tf):
        proj.devf = devf
        proj.tf = tf
        proj.icorr()
        proj.preprocessh(doicorr=True)
        #proj.optimize(variables=[1.46103037e-04, 1.23754329e-02, 5.20429443e-01],
        #proj.optimize(variables=[1.57991385e-04, 3.93973066e-03, 2.02717926e-01],
        #        figname="qens_hist_fit.png")
        proj.optimize(variables=[1.57991385e-04, 3.93973066e-03, 2.02717926e-01],
                figname=None)
        #proj.optimize(figname="qens_hist_fit_including_bg.png")
        proj.optimize(figname=None)


sqrun_kde_hist()

