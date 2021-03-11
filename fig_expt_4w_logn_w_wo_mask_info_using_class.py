#!/usr/bin/env python
import sys
sys.path.append("/home/kazu/ktpro")
import compare_result_files_4w_extrapolat_class as cc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import LogFormatterSciNotation, ScalarFormatter


class SCompare(cc.Compare):
    def __init__(self, infile1, label1, infile2, label2, dataname, stepsizes,
                 tcn, cn, ylabel, mnj, dshifti=None, dshiftf=None, log=True):
        self.log = log

        super(SCompare, self).__init__(infile1, label1, infile2, label2,
                                       dataname, stepsizes,
                                       tcn, cn, ylabel, mnj,
                                       dshifti=None, dshiftf=None)

    def create_fig(self):
        self.fig = plt.figure(figsize=(12, 9))
        self.fig.suptitle("Optimized bin-widths on experimental data " +
                          self.title)

    def plot_all_data(self):
        x = self.all_data[0, self.dshifti:self.dshiftf, 0]
        ys = self.all_data[:, self.dshifti:self.dshiftf, 1:]*self.stepsizes
        wlist = ["qx", "qy", "qz", "w"]

        for widx, wlabel in enumerate(wlist):
            ax = self.fig.add_subplot(4, self.tcn, self.tcn*widx+self.cn)
            ax.plot(x, ys[0, :, widx],
                    clip_on=False, linestyle="dotted", lw=1, label=self.label1,
                    marker="x", ms=8, mec='k', mfc='white', mew=1,  c='k')
            ax.plot(x, ys[1, :, widx],
                    clip_on=False, linestyle="dotted", lw=1, label=self.label2,
                    marker=".", ms=8, mec='k', mfc='w', mew=1, c='k')
            ax.set_ylim(0, np.max(ys[:, :, widx]) + self.stepsizes[widx])
            ax.yaxis.set_major_locator(MultipleLocator(self.stepsizes[widx]
                                                       * self.mnj))
            ax.yaxis.set_minor_locator(MultipleLocator(self.stepsizes[widx]))
            if self.log:
               ax.set_xscale('log')
               ax.text(np.max(x)*0.5, np.max(ys[:, :, widx])*0.9, wlabel)
               mergin = (np.max(np.log10(x))-np.min(np.log10(x)))*0.02
               ax.set_xlim(np.min(x)*10**(-mergin), np.max(x)*10**(mergin))
            else:
               ax.text(np.max(x)*0.7, np.max(ys[:, :, widx])*0.9, wlabel)
               mergin = (np.max(x)-np.min(x))*0.02
               ax.set_xlim(np.min(x)-mergin, np.max(x)+mergin)

            ax.xaxis.set_minor_formatter(plt.NullFormatter())
            ax.tick_params(top=True, right=True, direction='in', which='both')
            ax.tick_params(length=6, which='major')
            ax.tick_params(length=3, which='minor')
            ax.tick_params(labelbottom=False)
            if widx == 3:
                if self.log:
                    ax.set_xscale('log')
                ax.tick_params(labelbottom=True)
                ax.set_xlabel('total count')
            if widx == 1 and self.ylabel:
                ax.set_ylabel('bin width (rlu)')
            if widx == 3 and self.ylabel:
                ax.set_ylabel('bin width (meV)')
        plt.subplots_adjust(wspace=0.15, hspace=0.0)
        plt.legend(loc="best")


def run():
    infile1 = "/home/kazu/desktop/200204/" +\
            "fine/hourbyhour/ortho_opt_without_mask/result.txt_vec"
    infile2 = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
              "ortho_opt_without_mask/condparam09/result.txt_vec"
    label1 = r'$\alpha$=0'
    label2 = r'$\alpha$=0.9'
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    tcn = 3
    cn = 1
    ylabel = True
    mnj = 4
    mycc = SCompare(infile1, label1, infile2, label2, "17714 Ei42 Ei24",
                    stepsizes, tcn, cn, ylabel, mnj, log=True)
    mycc.create_fig()
    mycc.plot_all_data()
    mycc.cn = 2
    mycc.infile1 = "/home/kazu/desktop/200522/" +\
                   "Ei42/veryfineq/result.txt_vec"
    mycc.infile2 = "/home/kazu/desktop/200522/" +\
                   "Ei42/veryfineq/condparam_09/result.txt_vec"
    mycc.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    mycc.dshifti = None
    mycc.ylabel = False
    mycc.get_all_data()
    mycc.plot_all_data()
    mycc.cn = 3
    mycc.infile1 = "/home/kazu/desktop/200522/" +\
                   "Ei24/fineq/result.txt_vec"
    mycc.infile2 = "/home/kazu/desktop/200522/" +\
                   "Ei24/fineq/condparam07/result.txt_vec"
    mycc.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    mycc.ylabel = False
    mycc.label2 = r'$\alpha$=0.7'
    mycc.mnj = 5
    mycc.get_all_data()
    mycc.plot_all_data()

run()
plt.show()
