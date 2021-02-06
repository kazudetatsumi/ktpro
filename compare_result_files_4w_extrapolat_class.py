#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.ticker import LogFormatterSciNotation
import numpy as np


#class CustomTicker(LogFormatterSciNotation):
#    def __call__(self, x, pos=None):
#        return "{x:g}".format(x=x)


class Compare:

    def __init__(self, infile1, label1, infile2, label2, dataname, stepsizes,
                 dshifti=None, dshiftf=None):
        self.infile1 = infile1
        self.infile2 = infile2
        self.label1 = label1
        self.label2 = label2
        self.dataname = dataname
        self.stepsizes = stepsizes
        self.dshifti = dshifti
        self.dshiftf = dshiftf
        self.all_data()
        self.plot_all_data()

    def getdata(self, infile):
        nqx = []
        nqy = []
        nqz = []
        nw = []
        nn = []
        with open(infile, 'r') as f:
            for line in f:
                values = line.split()
                nqx.append(float(values[0]))
                nqy.append(float(values[1]))
                nqz.append(float(values[2]))
                nw.append(float(values[3]))
                nn.append(float(values[4]))
        return np.transpose(np.array([nn, nqx, nqy, nqz, nw]))

    def all_data(self):
        data1 = self.getdata(self.infile1)
        data2 = self.getdata(self.infile2)
        self.all_data = np.zeros((2, data1.shape[0], data1.shape[1]))
        self.all_data[0, :, :] = data1
        self.all_data[1, :, :] = data2

    def plot_all_data(self):
        fig = plt.figure(figsize=(6, 9))
        fig.suptitle("extrapolated bin-widths of 4D INS orthotope data "
                     + self.dataname)
        #x = np.log10(self.all_data[0, 2:, 0])
        x = self.all_data[0, self.dshifti:self.dshiftf, 0]
        ys = self.all_data[:, self.dshifti:self.dshiftf, 1:]*self.stepsizes
        wlist = ["qx", "qy", "qz", "w"]

        for widx, wlabel in enumerate(wlist):
            ax = fig.add_subplot(4, 1, widx+1)
            ax.plot(x, ys[1, :, widx],
                    clip_on=False, linestyle="dotted", lw=1, label=self.label2,
                    marker=".", ms=4, mec='k', mfc='k', c='k')
            ax.plot(x, ys[0, :, widx],
                    clip_on=False, linestyle="dotted", lw=1, label=self.label1,
                    marker="x", ms=4, mec='k', mfc='white', c='k')
            ax.text(np.max(x)*0.5, np.max(ys[:, :, widx])*0.9, wlabel)
            ax.tick_params(labelbottom=False)
            ax.set_ylim(0, np.max(ys[:, :, widx])*1.06)
            mergin = (np.max(np.log10(x))-np.min(np.log10(x)))*0.01
            ax.set_xlim(np.min(x)*10**(-mergin), np.max(x)*10**(mergin))
            ax.yaxis.set_major_locator(MultipleLocator(self.stepsizes[widx]*4.0))
            ax.yaxis.set_minor_locator(MultipleLocator(self.stepsizes[widx]))
            ax.xaxis.set_major_locator(MultipleLocator(1))
            ax.tick_params(top=True, right=True, direction='in', which='both')
            ax.tick_params(length=6, which='major')
            ax.tick_params(length=3, which='minor')
            ax.set_xscale('log')
            #ax.xaxis.set_major_formatter(CustomTicker())
            if widx == 3:
                ax.tick_params(labelbottom=True)
                ax.set_xlabel('total count')
            if widx == 1:
                ax.set_ylabel('bin width (rlu)')
            if widx == 3:
                ax.set_ylabel('bin width (meV)')
        plt.subplots_adjust(wspace=0.4, hspace=0.0)
        plt.legend(loc="center right")
        plt.show()


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
