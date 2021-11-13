#!/usr/bin/env python
# This script do fitting analysis on kde results and then on histograms.
# For kde results, the background is optimized.
# For histograms, the background is fixed as the bkg_peak_ratio(kde) * peak(histogram)
# as suggested by Dr. Matsuura.
# Kazuyoshi TATSUMI 2021/09/06
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
import qens_fit_class as qfc


def sqrun_kde_hist():
    if len(sys.argv) >= 2:
        dirname = sys.argv[1]
    else:
        print("using default dirname 000025io")
        dirname = "000025io"

    #head = "/home/kazu/desktop/210108/Tatsumi/pickles/"+dirname+"/"
    head = "./"
    devf = head + "qens_kde_o_divided_by_i_6204.pkl"
    tf = head + "qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    proj = qfc.qens_fit(devf, tf, elim, showplot=False)
    proj.icorr()
    proj.preprocess(doicorr=True)
    proj.optimize(figname="qens_kde_fit.png")

    head = "/home/kazu/desktop/210108/Tatsumi/srlz/"+dirname+"/"
    devf = head + "qens_hist_o_divided_by_i_6204.pkl"
    tf = head + "qens_hist_o_divided_by_i_6202.pkl"
    proj.devf = devf
    proj.tf = tf
    proj.icorr()
    proj.preprocessh(doicorr=True)
    proj.optimize(variables=[1.46103037e-04, 1.23754329e-02, 5.20429443e-01],
            figname="qens_hist_fit.png")
    proj.optimize(figname="qens_hist_fit_including_bg.png")


sqrun_kde_hist()

