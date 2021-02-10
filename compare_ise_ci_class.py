#!/usr/bin/env python
# statistics on ise related data of ise_result
# Kazuyoshi TATSUMI 2020.
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
import statsmodels.stats.api as sms
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)


class ise_ci:

    def __init__(self, head, isefn1, isefn2, dV, label1, label2, tcn, cn):
        self.head = head
        self.isefn1 = isefn1
        self.isefn2 = isefn2
        self.dV = dV
        self.label1 = label1
        self.label2 = label2
        self.tcn = tcn
        self.cn = cn
        self.get_all_data()

    def create_fig(self):
        fig = plt.figure(figsize=(16, 9))
        fig.suptitle("integrated squared errors with respect to total counts" +
                     " for phantom data of 17714 Ei42 Ei24")

    def getdata(self, isefn):
        for tryidx in range(1, 21):
            data = np.loadtxt(self.head + "try"+str(tryidx)+"/" + isefn,
                              delimiter=',')
            if tryidx == 1:
                tdata = np.zeros((20, data.shape[0], data.shape[1]))
            tdata[tryidx-1, :, :] = data
        return tdata

    def get_all_data(self):
        self.tdata1 = self.getdata(self.isefn1)
        self.tdata2 = self.getdata(self.isefn2)

    def plot_cis(self, shift=0.025):
        y1 = self.tdata1[:, :, 1]/self.dV
        y2 = self.tdata2[:, :, 1]/self.dV
        x12 = np.tile(np.average(np.concatenate([self.tdata1[:, :, 0],
                                                self.tdata2[:, :, 0]], axis=0),
                                 axis=0),
                      20).flatten()
        cis = sms.CompareMeans(sms.DescrStatsW(y1),
                               sms.DescrStatsW(y2)).tconfint_diff(
                               usevar='unequal')
        gs = gridspec.GridSpec(2, self.tcn, height_ratios=[3.2, 1])
        ax = plt.subplot(gs[0, self.cn])
        ax.scatter(x12*10.0**(shift), y1.flatten(), marker='_', s=20, c='k',
                   lw=1, label=self.label1)
        ax.scatter(x12*10.0**(-shift), y2.flatten(), marker='_', s=20, c='gray',
                   lw=1, label=self.label2)
        ax.set_xscale('log')
        ax.set_ylabel('integrated square errors ($rlu^{-3}meV^{-1}$)')
        miny = np.min([np.min(y1), np.min(y2)])
        maxy = np.max([np.max(y1), np.max(y2)])
        margin = (maxy - miny)*0.01
        ax.set_ylim(0, maxy+margin)
        ax.set_xlim(np.min(x12)*10.0**(-2.5*shift), np.max(x12)*10.0**(2.5*shift))
        plt.legend()
        ax = plt.subplot(gs[1, self.cn])
        x = np.average(self.tdata2[:, :, 0], axis=0)
        y = (cis[0] + cis[1])/2.0
        yerr = np.abs(cis[0] - cis[1])/2.0
        ax.errorbar(x, y,  yerr=yerr,  capsize=3, fmt="none", ms="5", c='k')
        ax.axhline(y=0, ls='--', c='k', lw=0.5)
        ax.set_xscale('log')
        margin = (np.max(cis[1]) - np.min(cis[0]))*0.03
        ax.set_ylim(np.min(cis[0])-margin, np.max(cis[1])+margin)
        ax.set_ylabel('mean difference ($rlu^{-3}meV^{-1}$)')
        ax.set_xlabel('total count')
        plt.subplots_adjust(wspace=0.2, hspace=0.0)


def run():
    head = "/home/kazu/desktop/200312/for_cu_new/old_orthotope/" +\
           "orthotope_ddscs_again2/expt_orthotope_bd/"
    isefn1 = "ise_withoutcond"
    isefn2 = "ise_condparam09"
    dV = 0.5*0.025*0.025*0.025
    label1 = r'$\alpha=0$'
    label2 = r'$\alpha=0.9$'
    tcn = 3
    cn = 0
    ic = ise_ci(head, isefn1, isefn2, dV, label1, label2, tcn, cn)
    ic.create_fig()
    ic.plot_cis(shift=0.034)

    ic.head = "/home/kazu/desktop/200701/orthotope_again_ddscs/"
    ic.dV = 0.2*0.0125*0.025*0.05
    ic.cn = 1
    ic.get_all_data()
    ic.plot_cis(shift=0.022)

    ic.head = "/home/kazu/desktop/200903/morethan10meV/" +\
              "samedataasfilled_again2_ddscs/"
    ic.isefn2 = "ise_condparam05"
    ic.dV = 0.08*0.01*0.01*0.04
    ic.label2 = r'$\alpha=0.5$'
    ic.cn = 2
    ic.get_all_data()
    ic.plot_cis(shift=0.01)

    plt.show()
    #tdata1, tdata2 = get_characterized_data(isefn1, isefn2)
    #plot_cis(tdata1, tdata2, dV)


run()
