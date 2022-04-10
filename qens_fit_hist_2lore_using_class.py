#!/usr/bin/env python
# This script do a fitting analysis using 2 Lorentzian components
# ONLY on uncorrected and reconstructed histogram data.
# The data can be prepared by
# qens_hist_results_odata_divided_by_idata_reconst_using_class.py
# Usage: qens_fit_hist_2lore_using_class.py 6202 6204 0.000025 0125
# Example data are found at
# /home/kazu/desktop/210108/Tatsumi/winparm_exaam/test_reconstruct_python
# although the job shell file has such a line
# qens_fit_hist_2lore_using_class.py 6202 6204 0000025io 0.000025 0125
# where 0000025io is deleted at the present version.
# Kazuyoshi TATSUMI 2022/03/30
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
    if len(sys.argv) >= 4:
        denew = sys.argv[3]
    if len(sys.argv) >= 5:
        frac = sys.argv[4]
    else:
        print("using default frac """)
        frac = ""

    head = "./"
    devf = head + "qens_hist_o_divided_by_i_"+runnod+denew+".pkl"
    tf = head + "qens_hist_o_divided_by_i_"+runno+frac+denew+".pkl"
    elim = [-0.03, 0.07]
    proj = qfc.qens_fit(devf, tf, elim, showplot=False)
    proj.icorr()
    proj.preprocessh(doicorr=True)
    proj.optimize(variables=[2.64670711e-04, 2.67444797e-02, 4.57745873e-05,
                             5.06836899e-03, 1.45026317e-01, 3.12477175e-02],
                  figname="qens_hist_fit2.png")


sqrun_kde_hist_2lore()
