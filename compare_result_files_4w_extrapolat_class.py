#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import LogFormatterSciNotation, ScalarFormatter
import numpy as np
import re


class CustomTicker(LogFormatterSciNotation):
    def __call__(self, x, pos=None):
        return "{x:g}".format(x=x)


class Compare:

    def __init__(self, infile1, label1, infile2, label2, title, stepsizes,
                 tcn, cn, ylabel, mnj, dshifti=None, dshiftf=None):
        self.infile1 = infile1
        self.infile2 = infile2
        self.label1 = label1
        self.label2 = label2
        self.title = title
        self.stepsizes = stepsizes
        self.tcn = tcn
        self.cn = cn
        self.ylabel = ylabel
        self.mnj = mnj
        self.dshifti = dshifti
        self.dshiftf = dshiftf
        self.get_all_data()
        #self.plot_all_data()

    def getdata(self, infile):
        nqx = []
        nqy = []
        nqz = []
        nw = []
        nn = []
        with open(infile, 'r') as f:
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

    def get_all_data(self):
        data1 = self.getdata(self.infile1)
        data2 = self.getdata(self.infile2)
        self.all_data = np.zeros((2, data1.shape[0], data1.shape[1]))
        self.all_data[0, :, :] = data1
        self.all_data[1, :, :] = data2

    def create_fig(self):
        self.fig = plt.figure(figsize=(12, 9))
        self.fig.suptitle(self.title)

    def plot_all_data(self, log=True, xlim=None, vlims=None, alpha=1.0):
        x = self.all_data[0, self.dshifti:self.dshiftf, 0]*alpha
        ys = self.all_data[:, self.dshifti:self.dshiftf, 1:]*self.stepsizes
        wlist = ["qx", "qy", "qz", "w"]

        for widx, wlabel in enumerate(wlist):
            ax = self.fig.add_subplot(4, self.tcn, self.tcn*widx+self.cn)
            ax.plot(x, ys[0, :, widx],
                    clip_on=False, linestyle="dotted", lw=0.8, label=self.label1,
                    marker="x", ms=5, mec='k', mfc='white', mew=0.8,  c='k')
            ax.plot(x, ys[1, :, widx],
                    clip_on=False, linestyle="dotted", lw=0.8, label=self.label2,
                    marker=".", ms=5, mec='k', mfc='w', mew=0.8, c='k')

            ax.set_ylim(0, np.max(ys[:, :, widx]) + self.stepsizes[widx])
            if log:
                ax.text(np.max(x)*0.5, np.max(ys[:, :, widx])*0.9, wlabel)
                mergin = (np.max(np.log10(x))-np.min(np.log10(x)))*0.02
                if xlim is not None:
                    ax.set_xlim((xlim[0]*alpha, xlim[1]*alpha))
                else:
                    ax.set_xlim(np.min(x)*10**(-mergin), np.max(x)*10**(mergin))
                #print(np.min(x)*10**(-mergin),  np.max(x)*10**(mergin))
                ax.set_xscale('log')
            else:
                ax.text((np.max(x)-np.min(x))*0.7, np.max(ys[:, :, widx])*0.9,
                        wlabel)
                mergin = (np.max(x)-np.min(x))*0.02
                ax.set_xlim(np.min(x)-mergin, np.max(x)+mergin)

            if vlims is not None:
                ax.axvline(x=vlims[0]*alpha, color='gray', lw=0.5)
                ax.axvline(x=vlims[1]*alpha, color='gray', lw=0.5)

            ax.yaxis.set_major_locator(MultipleLocator(self.stepsizes[widx]
                                                       * self.mnj))
            ax.yaxis.set_minor_locator(MultipleLocator(self.stepsizes[widx]))
            ax.xaxis.set_minor_formatter(plt.NullFormatter())
            ax.tick_params(top=True, right=True, direction='in', which='both')
            ax.tick_params(length=6, which='major')
            ax.tick_params(length=3, which='minor')
            ax.tick_params(labelbottom=False)
            if widx == 3:
                ax.tick_params(labelbottom=True)
                if abs(alpha - 1.0) > 0.0001:
                    ax.set_xlabel('count in common space')
                else:
                    ax.set_xlabel('total count')
            if widx == 1 and self.ylabel:
                ax.set_ylabel('bin width (rlu)')
            if widx == 3 and self.ylabel:
                ax.set_ylabel('bin width (meV)')
        plt.subplots_adjust(wspace=0.15, hspace=0.0)
        plt.legend(loc="best")
        #plt.show()


#def samplerun():
#    infile1 = "/home/kazu/desktop/200204/" +\
#            "fine/hourbyhour/ortho_opt_without_mask/result_only_extrapolate"
#    infile2 = "/home/kazu/desktop/200204/fine/hourbyhour/" +\
#              "ortho_opt_without_mask/condparam09/result_only_extrapolate"
#    label1 = r'$\alpha$=0'
#    label2 = r'$\alpha$=0.9'
#    stepsizes = np.array([0.025, 0.025, 0.025, 0.5])
#    projectset = Compare(infile1, label1, infile2, label2, "17714", stepsizes)

#samplerun()
