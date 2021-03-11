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
    tcn = 3
    cn = 1
    ylabel = True
    mnj = 4
    title = "Extrapolated optimal bin-widths on simulated data 17714 Ei42 Ei24"
    mycc = cc.Compare(infile1, label1, infile2, label2, title,
                      stepsizes, tcn, cn, ylabel, mnj, dshifti=2)
    mycc.create_fig()
    mycc.plot_all_data()
    mycc.cn = 2
    mycc.infile1 = "/home/kazu/desktop/200701/" +\
                   "orthotope_again_ddscs/result_only_extrapolate"
    mycc.infile2 = "/home/kazu/desktop/200701/" +\
                   "orthotope_again_ddscs/condparam09/result_only_extrapolate"
    mycc.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    mycc.dshifti = None
    mycc.dshiftf = -4
    mycc.ylabel = False
    mycc.get_all_data()
    mycc.plot_all_data()
    mycc.cn = 3
    mycc.infile1 = "/home/kazu/desktop/200903/morethan10meV/" +\
                   "samedataasfilled_again2_ddscs/result_only_extrapolate"
    mycc.infile2 = "/home/kazu/desktop/200903/morethan10meV/" +\
                   "samedataasfilled_again2_ddscs/condparam05/" +\
                   "result_only_extrapolate"
    mycc.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    mycc.dshiftf = -6
    mycc.ylabel = False
    mycc.label2 = r'$\alpha$=0.5'
    mycc.mnj = 5
    mycc.get_all_data()
    mycc.plot_all_data()

run()
plt.show()
