#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import compare_result_files_4w_extrapolat_class as cc
import numpy as np
import matplotlib.pyplot as plt
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


def run():
    infile2 = "/home/kazu/desktop/200204/" +\
            "fine/hourbyhour/ortho_opt_without_mask/result_only_extrapolate"
    infile1 = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
              "ortho_opt_without_mask/condparam09/result_only_extrapolate"
    label2 = r'$\alpha$=0'
    label1 = r'$\alpha$=0.9'
    label1 = 'expt'
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
    title = "Extrapolated optimal bin-widths on experimental and simulated " +\
            "data 17714 Ei42 Ei24"
    mycc = cc.Compare(infile1, label1, infile2, label2, title,
                      stepsizes, tcn, cn, ylabel, mnj, dshifti=2)
    mycc.create_fig()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha, onlyone=True,
                       c='b',marker='x')

    mycc.infile2 = "/home/kazu/desktop/200312/for_cu_new/old_orthotope" +\
                   "/orthotope_ddscs_kfki_debye/" + \
                   "result_only_extrapolate"
    #mycc.infile1 = "/home/kazu/desktop/200312/for_cu_new/old_orthotope" +\
    #               "/orthotope_ddscs_kfki_debye/condparam09/" +\
    #               "result_only_extrapolate"
    mycc.infile1 = "/home/kazu/desktop/200312/for_cu_new/old_orthotope" +\
                   "/orthotope_again_gamma_corrected/condparam09/" +\
                   "result_only_extrapolate"
    xlim = (1952.42556793021, 62110867.96516186)
    #tcountcommon = 371.0
    tcountcommon = 15591.0
    tcountcommon_expt = 123.0
    #tcount = 626836.0
    tcount = 86103350
    alpha = tcountcommon/tcount
    vlims_expt = np.array([5.155200e+04, 5.189990e+05])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    vlims = vlims_expt/alpha*alpha_expt
    xlim = (1.0/alpha, xlim[1]*alpha_expt/alpha)
    mycc.label1 = 'simu'
    mycc.get_all_data()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha, onlyone=True,
                       c='r', marker='.')

    mycc.cn = 2
    mycc.infile2 = "/home/kazu/desktop/200522/" +\
                   "Ei42/veryfineq/result_only_extrapolate"
    mycc.infile1 = "/home/kazu/desktop/200522/" +\
                   "Ei42/veryfineq/condparam_09/result_only_extrapolate"
    vlims = [1.140130e+05, 5.906150e+05]
    tcountcommon = 483.0
    alpha = tcountcommon/vlims[1]
    mycc.label1 = 'expt'
    mycc.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    mycc.dshifti = None
    mycc.dshiftf = -4
    mycc.ylabel = False
    mycc.get_all_data()
    mycc.plot_all_data(vlims=vlims, alpha=alpha, onlyone=True, marker='x',
                       c='b')

    mycc.infile2 = "/home/kazu/desktop/200701/" +\
                   "orthotope_ddscs_kfki_debye/result_only_extrapolate"
    #mycc.infile1 = "/home/kazu/desktop/200701/" +\
    #               "orthotope_ddscs_kfki_debye/condparam09/result_only_extrapolate"
    mycc.infile1 = "/home/kazu/desktop/200701/" +\
                   "orthotope_gamma_corrected/condparam09/result_only_extrapolate"
    mycc.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    xlim = (3956.076238453173, 27880800.646331012)
    tcountcommon = 34470.0
    tcountcommon_expt = 483.0
    tcount = 53504590
    alpha = tcountcommon/tcount
    vlims_expt = np.array([1.140130e+05, 5.906150e+05])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    xlim = (xlim[0]*alpha_expt/alpha, xlim[1]*alpha_expt/alpha)
    vlims = vlims_expt/alpha*alpha_expt
    mycc.label1 = 'simu'
    mycc.dshifti = None
    mycc.dshiftf = -4
    mycc.ylabel = False
    mycc.get_all_data()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha, onlyone=True,
                       marker='.', c='r')

    mycc.cn = 3
    mycc.infile2 = "/home/kazu/desktop/200522/" +\
                   "Ei24/fineq/result_only_extrapolate"
    mycc.infile1 = "/home/kazu/desktop/200522/" +\
                   "Ei24/fineq/condparam07/result_only_extrapolate"
    mycc.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    xlim = (654.519524723231, 3050611.069576253)
    vlims = [7.821000e+03, 9.685300e+04]
    tcountcommon = 1885.0
    alpha = tcountcommon/vlims[1]
    mycc.label1 = 'expt'
    mycc.dshiftf = -6
    mycc.ylabel = False
    mycc.label2 = r'$\alpha$=0.7'
    mycc.mnj = 5
    mycc.get_all_data()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha, onlyone=True,
                       marker='x', c='b')
    mycc.infile2 = "/home/kazu/desktop/200903/morethan10meV/" +\
                   "samedataasfilled_ddscs_kfki_debye/result_only_extrapolate"
    mycc.infile1 = "/home/kazu/desktop/200903/morethan10meV/" +\
                   "samedataasfilled_ddscs_kfki_debye/condparam05/" +\
                   "result_only_extrapolate"
    mycc.infile1 = "/home/kazu/desktop/200903/morethan10meV/" +\
                   "samedataasfilled_gamma_corrected/condparam05/" +\
                   "result_only_extrapolate"
    mycc.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    xlim = (654.519524723231, 3050611.069576253)
    tcountcommon = 4370.0
    tcountcommon_expt = 1885.0
    tcount = 288168.0
    alpha = tcountcommon/tcount
    vlims_expt = np.array([7.821000e+03, 9.685300e+04])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    xlim = (xlim[0]*alpha_expt/alpha, xlim[1]*alpha_expt/alpha)
    vlims = vlims_expt/alpha*alpha_expt
    mycc.label1 = 'simu'
    mycc.dshiftf = -6
    mycc.ylabel = False
    mycc.label2 = r'$\alpha$=0.5'
    mycc.mnj = 5
    mycc.get_all_data()
    mycc.plot_all_data(xlim=xlim, vlims=vlims, alpha=alpha, onlyone=True,
                       marker='.', c='r')


run()
plt.savefig("expt_theo_noddscs_extraploate.pdf")
plt.show()
