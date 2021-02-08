#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import compare_result_files_4w_extrapolat_class as cc
import numpy as np
import matplotlib.pyplot as plt



def run():
    infile1 = "/home/kazu/desktop/200204/" +\
            "fine/hourbyhour/ortho_opt_without_mask/result_only_extrapolate"
    infile2 = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
              "ortho_opt_without_mask/condparam09/result_only_extrapolate"
    label1 = r'$\alpha$=0'
    label2 = r'$\alpha$=0.9'
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    tcn = 3
    cn = 1
    ylabel = True
    mnj = 4
    mycc = cc.Compare(infile1, label1, infile2, label2, "17714 Ei42 Ei24",
                       stepsizes, tcn, cn, ylabel, mnj, dshifti=2)
    mycc.create_fig()
    mycc.plot_all_data()
    mycc.cn = 2
    mycc.infile1 = "/home/kazu/desktop/200522/" +\
                    "Ei42/veryfineq/result_only_extrapolate"
    mycc.infile2 = "/home/kazu/desktop/200522/" +\
                    "Ei42/veryfineq/condparam_09/result_only_extrapolate"
    mycc.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    mycc.dshifti = None
    mycc.dshiftf = -4
    mycc.ylabel = False
    mycc.get_all_data()
    mycc.plot_all_data()
    mycc.cn = 3
    mycc.infile1 = "/home/kazu/desktop/200522/" +\
                    "Ei24/fineq/result_only_extrapolate"
    mycc.infile2 = "/home/kazu/desktop/200522/" +\
                    "Ei24/fineq/condparam07/result_only_extrapolate"
    mycc.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    mycc.dshiftf = -6
    mycc.ylabel = False
    mycc.label2 = r'$\alpha$=0.7'
    mycc.mnj = 5
    mycc.get_all_data()
    mycc.plot_all_data()

run()
plt.show()
