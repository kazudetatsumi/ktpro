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

params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)


class ise_ci:

    def __init__(self, head, isefn1, isefn2, dV, label1, label2, tcn, cn,
                 infiler1, infiler2, stepsizes, mnj):
        self.head = head
        self.isefn1 = isefn1
        self.isefn2 = isefn2
        self.dV = dV
        self.label1 = label1
        self.label2 = label2
        self.tcn = tcn
        self.cn = cn
        self.infiler1 = infiler1
        self.infiler2 = infiler2
        self.stepsizes = stepsizes
        self.mnj = mnj
        self.get_all_data()

    def create_fig(self):
        fig = plt.figure(figsize=(16, 8))
        fig.suptitle("integrated squared errors with respect to total counts" +
                     " for phantom data of 17714 Ei42 Ei24")

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
        return tdata

    def get_all_data(self):
        self.tdata1 = self.getdata(self.isefn1)
        self.tdata2 = self.getdata(self.isefn2)
        self.tdatar1 = self.getdatar(self.infiler1)
        self.tdatar2 = self.getdatar(self.infiler2)
        self.maxdiff = np.max(np.max(np.abs(self.tdatar1 - self.tdatar2), axis=2), axis=0)
        #print(self.maxdiff)

    def plot_cis(self, shift=0.025, ylabel=True):
        y1 = self.tdata1[:, :, 1]/self.dV
        y2 = self.tdata2[:, :, 1]/self.dV
        x12 = np.tile(np.average(np.concatenate([self.tdata1[:, :, 0],
                                                self.tdata2[:, :, 0]], axis=0),
                                 axis=0),
                      20).flatten()
        #cis = sms.CompareMeans(sms.DescrStatsW(y1),
        #                       sms.DescrStatsW(y2)).tconfint_diff(
        #                       usevar='unequal')
        diff = y1 - y2
        cls = []
        cus = []
        for idx in range(0, y1.shape[1]):
            cl, cu = mybootstrap.bootstrap(diff[:, idx], 10000, np.mean, 0.95)
            cls.append(cl)
            cus.append(cu)
        cis = (np.array(cls), np.array(cus))
        #gs = gridspec.GridSpec(6, self.tcn, height_ratios=[3.2, 1, 2, 2, 2, 2])
        gs = gridspec.GridSpec(2, self.tcn, height_ratios=[3.2, 1])
        ax = plt.subplot(gs[0, self.cn])
        ax.scatter(x12*10.0**(shift), y1.flatten(), marker='_', s=2, c='red',
                   lw=1, label=self.label1)
        ax.scatter(x12*10.0**(-shift), y2.flatten(), marker='_', s=2, c='blue',
                   lw=1, label=self.label2)
        ax.set_xscale('log')
        if ylabel:
            ax.set_ylabel('integrated square errors ($rlu^{-3}meV^{-1}$)')
        miny = np.min([np.min(y1), np.min(y2)])
        maxy = np.max([np.max(y1), np.max(y2)])
        margin = (maxy - miny)*0.01
        ax.set_ylim(0, maxy+margin)
        ax.set_xlim(np.min(x12)*10.0**(-2.5*shift), np.max(x12)*10.0**(2.5*shift))
        ax.tick_params(top=True, right=True, direction='in', which='both')
        plt.legend()
        ax = plt.subplot(gs[1, self.cn])
        x = np.average(np.concatenate([self.tdata1[:, :, 0], self.tdata2[:, :, 0]], axis=0), axis=0)
        y = (cis[0] + cis[1])/2.0
        #yerr = np.abs(cis[0] - cis[1])/2.0
        yerr = np.stack([cis[0], cis[1]], axis=0)
        ax.errorbar(x, y,  yerr=yerr,  capsize=3, fmt="none", ms="5", c='k')
        ax.axhline(y=0, ls='--', c='k', lw=0.5)
        ax.set_xscale('log')
        margin = (np.max(cis[1]) - np.min(cis[0]))*0.03
        #ax.set_ylim(np.min(cis[0])-margin, np.max(cis[1])+margin)
        ax.set_xlim(np.min(x)*10.0**(-2.5*shift), np.max(x)*10.0**(2.5*shift))
        if ylabel:
            ax.set_ylabel('mean difference \n ($rlu^{-3}meV^{-1}$)')
        ax.set_xlabel('total count')
        ax.tick_params(top=True, right=True, direction='in', which='both')

        ys = np.stack([self.tdatar1[0, :, 1:], self.tdatar2[0, :, 1:]], axis=0)*self.stepsizes
        ys1 = np.array(self.tdatar1[:, :, 1:], dtype=int)
        ys2 = np.array(self.tdatar2[:, :, 1:], dtype=int)
        ysmf1 = np.zeros((ys1.shape[1], ys1.shape[2]))
        ysmf2 = np.zeros((ys2.shape[1], ys2.shape[2]))
        for idx1 in range(0, ys1.shape[1]):
            for idx2 in range(0, ys2.shape[2]):
                ysmf1[idx1, idx2] = np.argmax(np.bincount(ys1[:, idx1, idx2]))
                ysmf2[idx1, idx2] = np.argmax(np.bincount(ys2[:, idx1, idx2]))
        wlist = ["qx", "qy", "qz", "w"]
        ysmf1 = ysmf1*self.stepsizes
        ysmf2 = ysmf2*self.stepsizes
        #for widx, wlabel in enumerate(wlist):
        #    ax = plt.subplot(gs[2+widx, self.cn])
        #    #ax.plot(x, ys[0, :, widx],
        #    #        clip_on=False, linestyle="dotted", lw=0.8, label=self.label1,
        #    #        marker="x", ms=5, mec='k', mfc='white', mew=0.8,  c='k')
        #    #ax.plot(x, ys[1, :, widx],
        #    #        clip_on=False, linestyle="dotted", lw=0.8, label=self.label2,
        #    #        marker=".", ms=5, mec='k', mfc='w', mew=0.8, c='k')
        #    ax.scatter(x, ysmf1[:,widx], marker='x', s=20, c='k',
        #               lw=1, label=self.label1)
        #    ax.scatter(x, ysmf2[:,widx], marker='.', s=20, c='gray',
        #               lw=1, label=self.label2)
        #    ax.set_ylim(0, np.max(ys[:, :, widx]) + self.stepsizes[widx])
        #    ax.set_xlim(np.min(x)*10.0**(-2.5*shift), np.max(x)*10.0**(2.5*shift))
        #    ax.yaxis.set_major_locator(MultipleLocator(self.stepsizes[widx]
        #                                               * self.mnj))
        #    ax.yaxis.set_minor_locator(MultipleLocator(self.stepsizes[widx]))
        #    ax.tick_params(top=True, right=True, direction='in', which='both')
        #    ax.tick_params(length=6, which='major')
        #    ax.tick_params(length=3, which='minor')
        #    ax.tick_params(labelbottom=False)
        #    ax.set_xscale('log')
        plt.subplots_adjust(wspace=0.1, hspace=0.0)

def samplerun():
    head = "/home/kazu/desktop/200312/for_cu_new/old_orthotope/" +\
           "orthotope_ddscs_again2/expt_orthotope_bd/"
    isefn1 = "ise_withoutcond"
    isefn2 = "ise_condparam09"
    infiler1 = "result.txt_vec"
    infiler2 = "condparam09/result.txt_vec"
    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
    dV = 0.5*0.025*0.025*0.025
    label1 = r'$\alpha=0$'
    label2 = r'$\alpha=0.9$'
    tcn = 3
    cn = 0
    mnj = 4
    ic = ise_ci(head, isefn1, isefn2, dV, label1, label2, tcn, cn, infiler1, infiler2, stepsizes, mnj)
    ic.create_fig()
    ic.plot_cis(shift=0.034)

    ic.head = "/home/kazu/desktop/200701/orthotope_again_ddscs/"
    ic.stepsizes = np.array([0.0125, 0.025, 0.05, 0.2])
    ic.dV = 0.2*0.0125*0.025*0.05
    ic.cn = 1
    ic.get_all_data()
    ic.plot_cis(shift=0.022)

    ic.head = "/home/kazu/desktop/200903/morethan10meV/" +\
              "samedataasfilled_again2_ddscs/"
    ic.isefn2 = "ise_condparam05"
    ic.infiler2 = "condparam05/result.txt_vec"
    ic.stepsizes = np.array([0.01, 0.01, 0.04, 0.08])
    ic.dV = 0.08*0.01*0.01*0.04
    ic.label2 = r'$\alpha=0.5$'
    ic.cn = 2
    ic.mnj = 5
    ic.get_all_data()
    ic.plot_cis(shift=0.01)

    plt.show()
    #tdata1, tdata2 = get_characterized_data(isefn1, isefn2)
    #plot_cis(tdata1, tdata2, dV)


#samplerun()
