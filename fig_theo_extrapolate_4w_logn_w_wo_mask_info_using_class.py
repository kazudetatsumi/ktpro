#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import compare_result_files_4w_extrapolat_class as cc
import numpy as np
import matplotlib.pyplot as plt



def run():
    infile1 = "/home/kazu/desktop/200312/for_cu_new/old_orthotope" +\
              "/orthotope_ddscs_again2/expt_orthotope_bd/" + \
              "result_only_extrapolate"
    infile2 = "/home/kazu/desktop/200312/for_cu_new/old_orthotope" +\
              "/orthotope_ddscs_again2/expt_orthotope_bd/condparam09/" +\
              "result_only_extrapolate"
 
    label1 = r'$\alpha$=0'
    label2 = r'$\alpha$=0.9'
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    xlim = (1952.42556793021, 62110867.96516186)
    tcountcommon = 343.0
    tcountcommon_expt = 123.0
    tcount = 5.640260e+05
    alpha = tcountcommon/tcount
    vlims_expt = np.array([5.155200e+04, 5.189990e+05])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    vlims = vlims_expt/alpha*alpha_expt
    xlim = (1.0/alpha, xlim[1]*alpha_expt/alpha)
    tcn = 3
    cn = 1
    ylabel = True
    mnj = 4
    title = "Extrapolated optimal bin-widths on simulated data 17714 Ei42 Ei24"
    mycc = cc.Compare(infile1, label1, infile2, label2, title,
                      stepsizes, tcn, cn, ylabel, mnj, dshifti=2)
    mycc.create_fig()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha)

    mycc.cn = 2
    mycc.infile1 = "/home/kazu/desktop/200701/" +\
                   "orthotope_ddscs_kfki_debye/result_only_extrapolate"
    mycc.infile2 = "/home/kazu/desktop/200701/" +\
                   "orthotope_ddscs_kfki_debye/condparam09/result_only_extrapolate"
    mycc.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    xlim = (3956.076238453173, 27880800.646331012)
    #tcountcommon = 1424.0
    tcountcommon = 1610.0
    tcountcommon_expt = 483.0
    #tcount = 702056.0
    tcount = 797882.0
    alpha = tcountcommon/tcount
    vlims_expt = np.array([1.140130e+05, 5.906150e+05])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    xlim = (xlim[0]*alpha_expt/alpha, xlim[1]*alpha_expt/alpha)
    vlims = vlims_expt/alpha*alpha_expt
    mycc.dshifti = None
    mycc.dshiftf = -4
    mycc.ylabel = False
    mycc.get_all_data()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha)

    mycc.cn = 3
    mycc.infile1 = "/home/kazu/desktop/200903/morethan10meV/" +\
                   "samedataasfilled_ddscs_kfki_debye/result_only_extrapolate"
    mycc.infile2 = "/home/kazu/desktop/200903/morethan10meV/" +\
                   "samedataasfilled_ddscs_kfki_debye/condparam05/" +\
                   "result_only_extrapolate"
    mycc.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    xlim = (654.519524723231, 3050611.069576253)
    tcountcommon = 4324.0
    tcountcommon_expt = 1885.0
    tcount = 131514.0
    alpha = tcountcommon/tcount
    vlims_expt = np.array([7.821000e+03, 9.685300e+04])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    xlim = (xlim[0]*alpha_expt/alpha, xlim[1]*alpha_expt/alpha)
    vlims = vlims_expt/alpha*alpha_expt

    mycc.dshiftf = -6
    mycc.ylabel = False
    mycc.label2 = r'$\alpha$=0.5'
    mycc.mnj = 5
    mycc.get_all_data()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha)

run()
#plt.savefig("theo_extrapolate.pdf")
plt.show()
