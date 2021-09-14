#!/usr/bin/env python
# statistics on ise related data of ise_result
# Kazuyoshi TATSUMI 2020.
import sys
sys.path.append("/home/kazu/ktpro")
import compare_ise_ci_bin_class as cicbc
import numpy as np
import matplotlib.pyplot as plt


def run():
    head = "/home/kazu/desktop/200312/for_cu_new/old_orthotope/" +\
           "orthotope_ddscs_kfki_debye_rt/"
    isefn1 = "ise_withoutcond"
    isefn2 = "ise_condparam09"
    infiler1 = "result.txt_ise_lv"
    #infiler1 = "result.txt_ise_withcond"
    infiler2 = "condparam09/result.txt_vec"
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    label1 = r'$\alpha=0$'
    label2 = r'$\alpha=0.9$'
    tcn = 3
    cn = 0
    mnj = [6,12,6,2]
    tcount = 626836.0
    tcountcommon_expt = 123.0
    tcountcommon = 371.0
    alpha = tcountcommon/tcount
    vlims_expt = np.array([5.155200e+04, 5.189990e+05])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    vlims = vlims_expt/alpha*alpha_expt

    title = "bin-width differences from the searched bin-width" +\
            " for phantom data of 17714 Ei42 Ei24"

    ic = cicbc.ise_ci(head, isefn1, isefn2, label1, label2, tcn, cn, infiler1, infiler2, stepsizes, mnj)
    ic.create_fig(title)
    ic.plot_bin(alpha=alpha, vlims=vlims, rel=False)

    #ic.head = "/home/kazu/desktop/200701/orthotope_again_ddscs/"
    ic.head = "/home/kazu/desktop/200701/orthotope_ddscs_kfki_debye/"
    ic.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    ic.infiler1 = "result.txt_ise_lv"
    ic.cn = 1
    #tcount = 702056.0
    tcount = 796273.0
    tcountcommon_expt = 483.0
    #tcountcommon = 1424.0
    tcountcommon = 1669.0
    alpha = tcountcommon/tcount
    vlims_expt = np.array([1.140130e+05, 5.906150e+05])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    vlims = vlims_expt/alpha*alpha_expt
    ic.get_all_data()
    ic.mnj = [2,1,1,1]
    ic.plot_bin(ylabel=False, alpha=alpha, vlims=vlims, dels=[5, 7, 10, 11,15], rel=False)

    ic.head = "/home/kazu/desktop/200903/morethan10meV/" +\
              "samedataasfilled_ddscs_kfki_debye/"
    ic.isefn2 = "ise_condparam05"
    ic.infiler1 = "result.txt_ise_lv"
    ic.infiler2 = "condparam05/result.txt_vec"
    ic.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    ic.label2 = r'$\alpha=0.5$'
    ic.cn = 2
    ic.get_all_data()
    ic.mnj = [4,2,1,1]
    tcount = 131074.0
    tcountcommon_expt = 1885.0
    tcountcommon = 4342.0
    alpha = tcountcommon/tcount
    vlims_expt = np.array([7.821000e+03, 9.685300e+04])
    alpha_expt = tcountcommon_expt/vlims_expt[1]
    vlims = vlims_expt/alpha*alpha_expt
    ic.plot_bin(ylabel=False, alpha=alpha, vlims=vlims, dels=[9,11], rel=False )
    plt.savefig('binwidth_dif_rel.pdf')

    plt.show()


run()
