#!/usr/bin/env python
# statistics on ise related data of ise_result
# Kazuyoshi TATSUMI 2020.
import sys
sys.path.append("/home/kazu/ktpro")
import compare_bin_diff_diff_class as cbddc
import numpy as np
import matplotlib.pyplot as plt

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


def run():
    head = "/home/kazu/desktop/200312/for_cu_new/old_orthotope/" +\
           "orthotope_ddscs_again2/expt_orthotope_bd/"
    isefn1 = "ise_withoutcond"
    isefn2 = "ise_condparam09"
    isefn3 = "ise_searched"
    infiler1 = "result.txt_vec"
    infiler2 = "condparam09/result.txt_vec"
    infiler3 = "result.txt_ise_lv"
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    label1 = r'$\alpha=0$'
    label2 = r'$\alpha=0.9$'
    tcn = 3
    cn = 0
    mnj = [1, 1, 1, 1]
    title = "mean differences in deviations of the optimized bin-widths" +\
            " from the searched bin-widths" +\
            " between the methods with and without using the mask info"

    ic = cbddc.SIseci(head, isefn1, isefn2, isefn3, label1, label2, tcn, cn,
                      infiler1, infiler2, infiler3, stepsizes, mnj)
    ic.create_fig(title)
    ic.plot_bin()

    ic.head = "/home/kazu/desktop/200701/orthotope_again_ddscs/"
    ic.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    ic.cn = 1
    ic.get_all_data()
    ic.plot_bin(ylabel=False)

    ic.head = "/home/kazu/desktop/200903/morethan10meV/" +\
              "samedataasfilled_again2_ddscs/"
    ic.isefn2 = "ise_condparam05"
    ic.infiler2 = "condparam05/result.txt_vec"
    ic.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    ic.label2 = r'$\alpha=0.5$'
    ic.cn = 2
    ic.get_all_data()
    ic.plot_bin(ylabel=False, mdel=0)
    #plt.savefig('binwidth_dif_dif.pdf')

    plt.show()


run()
