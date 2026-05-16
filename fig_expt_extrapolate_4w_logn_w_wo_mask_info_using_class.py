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
    xlim = (1952.42556793021, 62110867.96516186)
    tcn = 3
    cn = 1
    ylabel = True
    mnj = 4
    vlims = [5.155200e+04, 5.189990e+05]
    tcountcommon = 123.0
    alpha = tcountcommon/vlims[1]
    xlim = (1.0/alpha, 62110867.96516186)
    title = "Extrapolated optimal bin-widths on experimental " +\
            "data 17714 Ei42 Ei24"
    mycc = cc.Compare(infile1, label1, infile2, label2, title,
                      stepsizes, tcn, cn, ylabel, mnj, dshifti=2)
    mycc.create_fig()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha)
    mycc.cn = 2
    mycc.infile1 = "/home/kazu/desktop/200522/" +\
                    "Ei42/veryfineq/result_only_extrapolate"
    mycc.infile2 = "/home/kazu/desktop/200522/" +\
                    "Ei42/veryfineq/condparam_09/result_only_extrapolate"
    vlims = [1.140130e+05, 5.906150e+05]
    tcountcommon = 483.0
    alpha = tcountcommon/vlims[1]
    mycc.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    mycc.dshifti = None
    mycc.dshiftf = -4
    mycc.ylabel = False
    mycc.get_all_data()
    mycc.plot_all_data(vlims=vlims, alpha=alpha)
    mycc.cn = 3
    mycc.infile1 = "/home/kazu/desktop/200522/" +\
                    "Ei24/fineq/result_only_extrapolate"
    mycc.infile2 = "/home/kazu/desktop/200522/" +\
                    "Ei24/fineq/condparam07/result_only_extrapolate"
    mycc.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    xlim = (654.519524723231, 3050611.069576253)
    vlims = [7.821000e+03, 9.685300e+04]
    tcountcommon = 1885.0
    alpha = tcountcommon/vlims[1]
    mycc.dshiftf = -6
    mycc.ylabel = False
    mycc.label2 = r'$\alpha$=0.7'
    mycc.mnj = 5
    mycc.get_all_data()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha)

run()
plt.savefig("expt_extraploate.eps")
plt.show()
