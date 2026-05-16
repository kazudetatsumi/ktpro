#!/usr/bin/env python
# This script plot the results of bin-width optimization, "reulst.txt_vec"
# to make Fig.1 for MLF annual report 2019.
# 2020 8/1 Kazuyoshi TATSUMI
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/home/kazu/ktpro')
import fig_expt_4w_logn_class as f4lc
import alpha_17714_Ei42_Ei24 as a1EE
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import LogFormatterSciNotation, ScalarFormatter

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['font.size'] = 12

plt.rcParams['font.family'] = 'Arial'




def run():
    num_clms = 3
    m = 2
    ulm = 0

    infile = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
             "ortho_opt_without_mask/condparam09/result.txt_vec"
    deltas = np.array([0.025, 0.025, 0.025, 0.5])
    clm = 1
    mnj = 4
    log = False
    projectset = f4lc.Plot_4w_Logn(infile, deltas, m, ulm, num_clms, mnj, log)
    #projectset.create_fig()
    projectset.fig = plt.figure(figsize=(5, 5))
    projectset.invplotter(clm, ylabel=True, alpha=a1EE.alpha('17714'))

    projectset.infile = "/home/kazu/desktop/200522/Ei42/veryfineq/" +\
             "condparam_09/result.txt_vec"
    projectset.deltas = np.array([0.0125, 0.025, 0.050, 0.2])
    clm = 2
    #projectset_Ei42 = Plot_4w_Logn(infile, deltas, m, ulm, num_clms)
    projectset.invplotter(clm, ylabel=False, alpha=a1EE.alpha('Ei42'))

    projectset.infile = "/home/kazu/desktop/200522/Ei24/fineq/" +\
             "/condparam07/result.txt_vec"
    projectset.deltas = np.array([0.01, 0.01, 0.04, 0.08])
    clm = 3
    projectset.mnj
    projectset.invplotter(clm, ylabel=False, alpha=a1EE.alpha('Ei24'))


run()
plt.savefig("fig_opt_ext_inv_binwidths.pdf")
plt.show()
