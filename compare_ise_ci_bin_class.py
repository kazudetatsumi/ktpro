#!/usr/bin/env python
# statistics on ise related data of ise_result
# Kazuyoshi TATSUMI 2020.
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
import mybootstrap
from matplotlib import pyplot as plt
from matplotlib import gridspec
#import statsmodels.stats.api as sms
import re
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import LogFormatterSciNotation, ScalarFormatter

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


class ise_ci:

    def __init__(self, head, isefn1, isefn2, label1, label2, tcn, cn,
                 infiler1, infiler2, stepsizes, mnj):
        self.head = head
        self.isefn1 = isefn1
        self.isefn2 = isefn2
        self.label1 = label1
        self.label2 = label2
        self.tcn = tcn
        self.cn = cn
        self.infiler1 = infiler1
        self.infiler2 = infiler2
        self.stepsizes = stepsizes
        self.mnj = mnj
        self.get_all_data()

    def create_fig(self, title=None):
        self.fig = plt.figure(figsize=(16, 12))
        if title is None:
            self.fig.suptitle("integrated squared errors with respect to total" +
                         "counts for phantom data of 17714 Ei42 Ei24")
        else:
            self.fig.suptitle(title)

    def getdata(self, isefn):
        data = [np.loadtxt(self.head + "try" + str(tryidx) + "/" +
                           isefn, delimiter=',')
                for tryidx in range(1, 21)]
        tdata = np.stack(data, axis=0)
        return tdata

    def getdatar_sub(self, fulinfile):
        nqx = []
        nqy = []
        nqz = []
        nw = []
        nn = []
        with open(fulinfile, 'r') as f:
            for line in f:
                if re.compile(r'[A-z]{2,}').search(line) and len(nqx) > 0:
                    break
                if not re.compile(r'[A-z]{2,}').search(line):
                    values = line.split()
                    if re.compile(r'[A-z]').search(values[0]):
                        nqx.append(float(values[1]))
                        nqy.append(float(values[2]))
                        nqz.append(float(values[3]))
                        nw.append(float(values[4]))
                        nn.append(float(values[0]))
                    else:
                        nqx.append(float(values[0]))
                        nqy.append(float(values[1]))
                        nqz.append(float(values[2]))
                        nw.append(float(values[3]))
                        nn.append(float(values[4]))
        return np.transpose(np.array([nn, nqx, nqy, nqz, nw]))

    def getdatar(self, infile):
        data = [self.getdatar_sub(self.head + "try" + str(tryidx) + "/" + infile)
                for tryidx in range(1, 21)]
        tdata = np.stack(data, axis=0)
        return np.array(tdata, dtype='int32')

    def get_all_data(self):
        self.tdata1 = self.getdata(self.isefn1)
        self.tdata2 = self.getdata(self.isefn2)
        self.tdatar1 = self.getdatar(self.infiler1)
        self.tdatar2 = self.getdatar(self.infiler2)
        self.maxdiff = np.max(np.max(np.abs(self.tdatar1 - self.tdatar2),
                                     axis=2), axis=0)

    def plot_cis(self, shift=0.025, ylabel=True, mdel=0):
        y1 = self.tdata1[:, mdel:, 1]/np.prod(self.stepsizes)
        y2 = self.tdata2[:, mdel:, 1]/np.prod(self.stepsizes)
        x12 = np.tile(np.average(np.concatenate([self.tdata1[:, mdel:, 0],
                                                self.tdata2[:, mdel:, 0]], axis=0),
                                 axis=0),
                      20).flatten()
        #cis = sms.CompareMeans(sms.DescrStatsW(y1),
        #                       sms.DescrStatsW(y2)).tconfint_diff(
        #                       usevar='unequal')
        diff = y2 - y1
        cls = []
        cus = []
        for idx in range(0, y1.shape[1]):
            cl, cu = mybootstrap.bootstrap(diff[:, idx], 10000, np.mean, 0.95)
            cls.append(cl)
            cus.append(cu)
        cls = np.array(cls)
        cus = np.array(cus)
        yerr = cus - cls
        y = (cls+cus)/2.0
        gs = gridspec.GridSpec(2, self.tcn, height_ratios=[3.2, 1])
        ax = plt.subplot(gs[0, self.cn])
        ax.scatter(x12*10.0**(shift), y1.flatten(), marker='_', s=2, c='red',
                   lw=1, label=self.label1)
        ax.scatter(x12*10.0**(-shift), y2.flatten(), marker='_', s=2, c='blue',
                   lw=1, label=self.label2)
        ax.set_xscale('log')
        if ylabel:
            ax.set_ylabel('$\overline{ISE}$ ($rlu^{-3}meV^{-1}$)')
        miny = np.min([np.min(y1), np.min(y2)])
        maxy = np.max([np.max(y1), np.max(y2)])
        margin = (maxy - miny)*0.05
        ax.set_ylim(0, maxy+margin)
        ax.set_xlim(np.min(x12)*10.0**(-2.5*shift), np.max(x12)*10.0**(2.5*shift))
        ax.tick_params(top=True, right=True, direction='in', which='both')
        plt.legend()
        ax = plt.subplot(gs[1, self.cn])
        x = np.average(np.concatenate([self.tdata1[:, mdel:, 0],
                                       self.tdata2[:, mdel:, 0]], axis=0), axis=0)
        ax.errorbar(x, y,  yerr=yerr,  capsize=3, fmt="none", ms="5", c='k')
        ax.axhline(y=0, ls='--', c='k', lw=0.5)
        ax.set_xscale('log')
        ax.set_xlim(np.min(x)*10.0**(-2.5*shift), np.max(x)*10.0**(2.5*shift))
        if ylabel:
            ax.set_ylabel('mean difference \n ($rlu^{-3}meV^{-1}$)')
        ax.set_xlabel('total count')
        ax.tick_params(top=True, right=True, direction='in', which='both')

    def plot_bin(self, ylabel=True, mdel=0):
        gs = gridspec.GridSpec(4, self.tcn, height_ratios=[1, 1, 1, 1])
        ys1 = np.array(self.tdatar1[:, :, 1:], dtype=int)
        ys2 = np.array(self.tdatar2[:, :, 1:], dtype=int)
        xs = np.average(self.tdatar1[:, :, 0], axis=0)
        diff = ys1 - ys2
        wlist = ["qx", "qy", "qz", "w"]
        for widx, wlabel in enumerate(wlist):
            ax = plt.subplot(gs[widx, self.cn])
            for nidx in range(mdel, diff.shape[1]):
                ue, uc = np.unique(diff[:, nidx, widx], return_counts=True)
                uc = uc/(diff.shape[0]*1.0)
                ue = ue*self.stepsizes[widx]
                xse = np.ones_like(ue)*xs[nidx]
                c = ax.scatter(xse, ue, c=uc, ec='k', lw=1.5, cmap='binary',
                               vmin=0.0, vmax=1.0, s=60)

            ax.yaxis.set_major_locator(MultipleLocator(
                                       self.stepsizes[widx]*self.mnj[widx]))
            if self.mnj[widx] != 1:
                ax.yaxis.set_minor_locator(MultipleLocator(self.stepsizes[widx]))
            ax.tick_params(top=True, right=True, direction='in', which='both', width=2)
            ax.tick_params(labelbottom=True, labelleft=True)
            ax.axhline(y=0, ls='--', c='k', lw=0.5)
            ax.set_xscale('log')
            if ylabel and widx == 1:
                ax.set_ylabel('bin-width difference (rlu)')
            if ylabel and widx == 3:
                ax.set_ylabel('bin-width difference (meV)')
            self.fig.colorbar(c, ax=ax)
        plt.subplots_adjust(wspace=0.16, hspace=0.16)
        ax.tick_params(labelbottom=True)
        ax.set_xlabel('total count')


def samplerunbin():
    head = "/home/kazu/desktop/200312/for_cu_new/old_orthotope/" +\
           "orthotope_ddscs_again2/expt_orthotope_bd/"
    isefn1 = "ise_withoutcond"
    isefn2 = "ise_condparam09"
    infiler1 = "result.txt_ise"
    infiler2 = "condparam09/result.txt_vec"
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    label1 = r'$\alpha=0$'
    label2 = r'$\alpha=0.9$'
    tcn = 3
    cn = 0
    mnj = [1, 2, 2, 2]
    title = "bin-width differences from the searched bin-width" +\
            " for phantom data of 17714 Ei42 Ei24"

    ic = ise_ci(head, isefn1, isefn2, label1, label2, tcn, cn, infiler1,
                infiler2, stepsizes, mnj)
    ic.create_fig(title)
    ic.plot_bin()

    ic.head = "/home/kazu/desktop/200701/orthotope_again_ddscs/"
    ic.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    ic.cn = 1
    ic.infiler1 = "result.txt_ise"
    ic.infiler2 = "condparam09/result.txt_vec"
    ic.get_all_data()
    ic.mnj = [1, 1, 1, 1]
    ic.plot_bin(ylabel=False)

    ic.head = "/home/kazu/desktop/200903/morethan10meV/" +\
              "samedataasfilled_again2_ddscs/"
    ic.isefn2 = "ise_condparam05"
    ic.infiler1 = "result.txt_ise"
    ic.infiler2 = "condparam05/result.txt_vec"
    ic.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    ic.label2 = r'$\alpha=0.5$'
    ic.cn = 2
    ic.get_all_data()
    ic.mnj = [7, 2, 1, 4]
    ic.plot_bin(ylabel=False)

    plt.show()

def samplerunise():
    head = "/home/kazu/desktop/200312/for_cu_new/old_orthotope/" +\
           "orthotope_ddscs_again2/expt_orthotope_bd/"
    isefn1 = "ise_searched"
    isefn2 = "ise_condparam09"
    infiler1 = "result.txt_vec"
    infiler2 = "condparam09/result.txt_vec"
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    label1 = r'$\alpha=0$'
    label2 = r'$\alpha=0.9$'
    tcn = 3
    cn = 0
    mnj = 4
    ic = ise_ci(head, isefn1, isefn2, label1, label2, tcn, cn,
                infiler1, infiler2, stepsizes, mnj)
    ic.create_fig()
    ic.plot_cis(shift=0.034)

    ic.head = "/home/kazu/desktop/200701/orthotope_again_ddscs/"
    ic.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    ic.cn = 1
    ic.get_all_data()
    ic.plot_cis(shift=0.006, ylabel=False)

    ic.head = "/home/kazu/desktop/200903/morethan10meV/" +\
              "samedataasfilled_again2_ddscs/"
    ic.isefn2 = "ise_condparam05"
    ic.infiler2 = "condparam05/result.txt_vec"
    ic.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    ic.label2 = r'$\alpha=0.5$'
    ic.cn = 2
    ic.mnj = 5
    ic.get_all_data()
    ic.plot_cis(shift=0.01, ylabel=False)

    plt.show()

#samplerunbin()
#samplerunise()
