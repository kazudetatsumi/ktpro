#!/usr/bin/env python
# statistics on ise related data of ise_result
# Kazuyoshi TATSUMI 2020.
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import statsmodels.stats.api as sms
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)
plt.rcParams['font.size'] = 12

fig = plt.figure(figsize=(6, 9))
plt.rcParams['font.family'] = 'Arial'
fig.suptitle("integrated squared errors with respect to total counts")
head = "/home/kazu/desktop/200312/for_cu_new/old_orthotope/" +\
        "orthotope_ddscs_again2/expt_orthotope_bd/"


def get_characterized_data(isefn1, isefn2):
    for tryidx in range(1, 21):
        data1 = np.loadtxt(head + "try"+str(tryidx)+"/"+isefn1, delimiter=',')
        data2 = np.loadtxt(head + "try"+str(tryidx)+"/"+isefn2, delimiter=',')
        #data_f = np.loadtxt("try"+str(tryidx)+"/ise_filled", delimiter=',')
        if tryidx == 1:
            # the shape of the following arrays:
            # # of dataset, # of tcounts, # of prop
            #                              (tcount, ise, avepdfs, fracise)
            tdata1 = np.zeros((20, data1.shape[0], data1.shape[1]))
            tdata2 = np.zeros((20, data2.shape[0], data2.shape[1]))
        tdata1[tryidx-1, :, :] = data1
        tdata2[tryidx-1, :, :] = data2

    return(tdata1, tdata2)


def plot_cis(tdata1, tdata2, dV):
    tdata2[:, :, 0] = np.log10(tdata2[:, :, 0])
    tdata2[:, :, 1] = tdata2[:, :, 1]/dV
    tdata1[:, :, 1] = tdata1[:, :, 1]/dV
    cis = sms.CompareMeans(sms.DescrStatsW(tdata1[:, :, 1]),
                           sms.DescrStatsW(tdata2[:, :, 1])).tconfint_diff(
                                   usevar='unequal')
    gs = gridspec.GridSpec(2, 1, height_ratios=[3.2, 1])
    #ax = fig.add_subplot(2,1,1)
    ax = plt.subplot(gs[0])
    ax.scatter(np.tile(np.average(tdata2[:, :, 0], axis=0), 20).flatten() -
               0.030, tdata2[:, :, 1].flatten(),
               marker='_', s=30, c='k', lw=1, label=r'$\alpha=0.9$')
    ax.scatter(np.tile(np.average(tdata2[:, :, 0], axis=0), 20).flatten() +
               0.030, tdata1[:, :, 1].flatten(),
               marker='_', s=30, c='gray', lw=1, label=r'$\alpha=0.0$')
    ax.set_ylabel('integrated square errors ($rlu^{-3}meV^{-1}$)')
    miny = np.min([np.min(tdata1[:, :, 1]), np.min(tdata2[:, :, 1])])
    maxy = np.max([np.max(tdata1[:, :, 1]), np.max(tdata2[:, :, 1])])
    margin = (maxy - miny)*0.01
    ax.set_ylim(0, maxy+margin)
    plt.legend()
    ax = plt.subplot(gs[1])
    x = np.average(tdata2[:, :, 0], axis=0)
    y = (cis[0] + cis[1])/2.0
    yerr = np.abs(cis[0] - cis[1])/2.0
    ax.errorbar(x, y,  yerr=yerr,  capsize=3, fmt="none", ms="5", c='k',
                label=r"ISEs of $\alpha=0.0$  rel. to ISEs $\alpha=0.9$")
    ax.axhline(y=0, ls='--', c='k', lw=0.5)
    plt.legend()
    margin = (np.max(cis[1]) - np.min(cis[0]))*0.03
    ax.set_ylim(np.min(cis[0])-margin, np.max(cis[1])+margin)
    ax.set_ylabel('95 % CI ($rlu^{-3}meV^{-1}$)')
    ax.set_xlabel('log10(total count)')
    plt.subplots_adjust(wspace=0.4, hspace=0.0)
    plt.savefig("fig3_of_qbsf.pdf")
    plt.show()


def run():
    isefn1 = "ise_withoutcond"
    isefn2 = "ise_condparam09"
    dV = 0.5*0.025*0.025*0.025
    tdata1, tdata2 = get_characterized_data(isefn1, isefn2)
    plot_cis(tdata1, tdata2, dV)


run()
