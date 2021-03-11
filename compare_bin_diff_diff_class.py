#!/usr/bin/env python
# statistics on ise related data of ise_result
# Kazuyoshi TATSUMI 2020.
import sys
sys.path.append("/home/kazu/ktpro")
import mybootstrap
import compare_ise_ci_bin_class as cc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import LogFormatterSciNotation, ScalarFormatter


class SIseci(cc.ise_ci):
    def __init__(self, head, isefn1, isefn2, isefn3, label1, label2, tcn,
                 cn, infiler1, infiler2, infiler3, stepsizes, mnj):
        self.isefn3 = isefn3
        self.infiler3 = infiler3

        super(SIseci, self).__init__(head, isefn1, isefn2, label1,
                                     label2, tcn, cn, infiler1, infiler2,
                                     stepsizes, mnj)

    def get_all_data(self):
        self.tdatar1 = self.getdatar(self.infiler1)
        self.tdatar2 = self.getdatar(self.infiler2)
        self.tdatar3 = self.getdatar(self.infiler3)

    def diff_plotter(self, ax, diff, x, stepsize, marker, size, ec):
        ue, uc = np.unique(diff, return_counts=True)
        uc = uc/(diff.shape[0]*1.0)
        ue = ue*stepsize
        xse = np.ones_like(ue)*x
        ax.scatter(xse, ue, c=uc, ec=ec, lw=0.5, cmap='binary',
                   vmin=0.0, vmax=1.0, s=size, marker=marker)

    def plot_bin(self, ylabel=True, mdel=0):
        gs = gridspec.GridSpec(4, self.tcn, height_ratios=[1, 1, 1, 1])
        ys1 = np.array(self.tdatar1[:, mdel:, 1:], dtype=int)
        ys2 = np.array(self.tdatar2[:, mdel:, 1:], dtype=int)
        ys3 = np.array(self.tdatar3[:, mdel:, 1:], dtype=int)
        xs = np.average(self.tdatar1[:, mdel:, 0], axis=0)
        diff13 = np.abs(ys1 - ys3)
        diff23 = np.abs(ys2 - ys3)
        diff = (diff23 - diff13)*self.stepsizes
        wlist = ["qx", "qy", "qz", "w"]
        for widx, wlabel in enumerate(wlist):
            ax = plt.subplot(gs[widx, self.cn])
            y = []
            yerr = []
            for nidx in range(0, diff13.shape[1]):
                cl, cu = mybootstrap.bootstrap(diff[:, nidx, widx], 10000, np.mean, 0.95)
                yerr.append((cu - cl)/2.0)
                y.append((cl+cu)/2.0)
            y = np.array(y)
            yerr = np.array(yerr)
            ax.errorbar(xs, y, yerr=yerr, capsize=3, fmt="none", ms="5", c="k", elinewidth=3, capthick=3)
            if np.max(y+yerr) % self.stepsizes[widx] > 0.00001:
                ymax = (np.max(y+yerr)//self.stepsizes[widx] + 1)\
                        * self.stepsizes[widx]
            else:
                ymax = (np.max(y+yerr)//self.stepsizes[widx] + 0.10)\
                        * self.stepsizes[widx]
            if np.min(y-yerr) % self.stepsizes[widx] > 0.00001:
                ymin = (np.min(y-yerr)//self.stepsizes[widx])\
                        * self.stepsizes[widx]
            else:
                ymin = (np.min(y-yerr)//self.stepsizes[widx] - 1 - 0.10)\
                        * self.stepsizes[widx]
            ax.set_ylim(ymin, ymax)
            ax.yaxis.set_major_locator(MultipleLocator(
                                       self.stepsizes[widx]*self.mnj[widx]))
            ax.tick_params(top=True, right=True, direction='in', which='both', width=1.5)
            ax.tick_params(labelbottom=False)
            ax.axhline(y=0, ls='--', c='k', lw=0.5)
            ax.set_xscale('log')
            if ylabel and widx == 1:
                ax.set_ylabel('bin-width difference (rlu)')
            if ylabel and widx == 3:
                ax.set_ylabel('bin-width difference (meV)')
        plt.subplots_adjust(wspace=0.16, hspace=0.0)
        ax.tick_params(labelbottom=True)
        ax.set_xlabel('total count')


def samplerun():
    head = "/home/kazu/desktop/200312/for_cu_new/old_orthotope/" +\
           "orthotope_ddscs_again2/expt_orthotope_bd/"
    isefn1 = "ise_withoutcond"
    isefn2 = "ise_condparam09"
    isefn3 = "ise_searched"
    infiler1 = "result.txt_vec"
    infiler2 = "condparam09/result.txt_vec"
    infiler3 = "result.txt_ise"
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    label1 = r'$\alpha=0$'
    label2 = r'$\alpha=0.9$'
    tcn = 3
    cn = 0
    mnj = [1, 1, 1, 1]
    title = "mean differences in deviations of the optimized bin-widths from the searched bin-widths" +\
            " between the methods with and without using the mask info"

    ic = SIseci(head, isefn1, isefn2, isefn3, label1, label2, tcn, cn,
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
    ic.plot_bin(ylabel=False)

    plt.show()


#samplerun()
