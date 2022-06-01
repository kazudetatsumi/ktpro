#!/usr/bin/env python
# This script do fitting analysis on kde results and then on histograms.
# For kde results, the background is optimized.
# For histograms, the background is fixed as the bkg_peak_ratio(kde) * peak(histogram)
# as suggested by Dr. Matsuura.
# This is for the fitting analysis using 2 Lorentzian components. 
# example usage qens_fit_kde_vs_hist_2lore_using_class.py 6202 6204 0000025io 0875
# Kazuyoshi TATSUMI 2021/12/05
import numpy as np
import sys 
sys.path.append("/home/kazu/ktpro")
import qens_fit_class as qfc 


def sqrun_kde_hist_2lore():
    np.set_printoptions(linewidth=120)

    #if len(sys.argv) >= 2:
    #    dirname = sys.argv[1]
    #else:
    #    print("using default dirname 000025io")
    #    dirname = "0000025io"
    #if len(sys.argv) >= 3:
    #    frac = sys.argv[2]
    #else:
    #    print("using default frac """)
    #    frac = ""

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
    #proj.optimize(variables=[1.67780642e-04, 3.08144561e-02, 1.16049540e-04, 7.85885231e-03, 2.04665269e-01, 4.17453525e+00], figname="qens_kde_fit2.png")
    #proj.optimize(variables=[2.18704786e-04, 1.67980295e-02, 4.92405238e-05, 1.88866588e-03, 1.21127501e-01, 5.02759930e-02], figname="qens_kde_fit2.png")
    #proj.optimize(variables=[5.878804786e-01, 2.73880295e-02, 3.80705238e-01, 7.05866588e-03, 2.00127501e-01, 1.81659930e-02], figname="qens_kde_fit2.png")
    #proj.optimize(variables=[5.378804786e-01, 4.85880295e-03, 3.00705238e-01, 3.75866588e-04, 9.10127501e-04, 6.37659930e-06], figname="qens_kde_fit2.png")
    proj.optimize(variables=[6.452704786e-01, 4.58950295e-03, 3.77605238e-01, 3.06766588e-04, 2.27927501e-02, 1.25659930e-02], figname="qens_kde_fit2.png")
    #proj.optimize(variables=[3.852704786e-01, 1.35950295e-02, 4.47927501e-01, 2.25659930e-01], figname="qens_kde_fit2.png")


    head = "/home/kazu/desktop/210108/Tatsumi/srlz/"+dirname+"/"
    devf = head + "qens_hist_o_divided_by_i_"+runnod+".pkl"
    tf = head + "qens_hist_o_divided_by_i_"+runno+frac+".pkl"
    proj.devf = devf
    proj.tf = tf
    proj.icorr()
    proj.preprocessh(doicorr=True)
    #proj.optimize(variables=[1.46103037e-04, 1.23754329e-02, 5.20429443e-01],
    #        figname="qens_hist_fit.png")
    #proj.optimize(variables=[1.35639516e-05, 2.61566022e-02, 2.16304425e-06, 4.97615634e-03, 7.62240950e-02], figname="qens_hist_fit2_basefixed.png")
    #proj.optimize(variables=[2.64670711e-04, 2.67444797e-02, 4.57745873e-05, 5.06836899e-03, 1.45026317e-01, 3.12477175e-02], figname="qens_hist_fit2.png")
    proj.optimize(variables=[5.85639516e-01, 2.61566022e-02, 3.76304425e-01, 7.07615634e-03, 2.02240950e-01], figname="qens_hist_fit2_basefixed.png")
    #proj.optimize(variables=[5.84670711e-01, 4.76444797e-03, 2.99745873e-01, 3.76906899e-04, 9.18826317e-04, 6.38877175e-06], figname="qens_hist_fit2.png")
    #proj.optimize(variables=[5.878804786e-01, 2.73880295e-02, 3.80705238e-01, 7.05866588e-03, 2.00127501e-01, 1.81659930e-02], figname="qens_hist_fit2.png")
    proj.optimize(variables=[6.452704786e-01, 4.58950295e-03, 3.77605238e-01, 3.06766588e-04, 2.27927501e-02, 1.25659930e-02], figname="qens_hist_fit2.png")


sqrun_kde_hist_2lore()

