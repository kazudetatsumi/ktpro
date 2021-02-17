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
           "orthotope_ddscs_again2/expt_orthotope_bd/"
    isefn1 = "ise_withoutcond"
    isefn2 = "ise_condparam09"
    infiler1 = "result.txt_ise"
    infiler2 = "condparam09/result.txt_vec"
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    dV = 0.5*0.025*0.025*0.025
    label1 = r'$\alpha=0$'
    label2 = r'$\alpha=0.9$'
    tcn = 3
    cn = 0
    mnj = [1,2,2,2]
    title = "bin-width differences from the searched bin-width" +\
            " for phantom data of 17714 Ei42 Ei24"

    ic = cicbc.ise_ci(head, isefn1, isefn2, dV, label1, label2, tcn, cn, infiler1, infiler2, stepsizes, mnj)
    ic.create_fig(title)
    ic.plot_bin()

    ic.head = "/home/kazu/desktop/200701/orthotope_again_ddscs/"
    ic.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    ic.dV = 0.2*0.0125*0.025*0.05
    ic.cn = 1
    ic.infiler1 = "result.txt_ise"
    ic.infiler2 = "condparam09/result.txt_vec"
    ic.get_all_data()
    ic.mnj = [1,1,1,1]
    ic.plot_bin(ylabel=False)

    ic.head = "/home/kazu/desktop/200903/morethan10meV/" +\
              "samedataasfilled_again2_ddscs/"
    ic.isefn2 = "ise_condparam05"
    ic.infiler1 = "result.txt_ise"
    ic.infiler2 = "condparam05/result.txt_vec"
    ic.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    ic.dV = 0.08*0.01*0.01*0.04
    ic.label2 = r'$\alpha=0.5$'
    ic.cn = 2
    ic.get_all_data()
    ic.mnj = [7,2,1,4]
    ic.plot_bin(ylabel=False)

    plt.show()


run()
