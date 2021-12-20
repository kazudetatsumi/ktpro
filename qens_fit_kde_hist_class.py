#!/usr/bin/env python
# This script do fitting analysis on kde results and then on histograms.
# For kde results, the background is optimized.
# For histograms, the background is fixed as
# the bkg_peak_ratio(kde) * peak(histogram)
# as suggested by Dr. Matsuura.
# Kazuyoshi TATSUMI 2021/09/06
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
import qens_fit_class as qfc
import matplotlib.pyplot as plt


class sqrun_kde_hist(qfc.qens_fit):
    def __init__(self, elim=[-0.03, 0.10]):
        self.elim = elim
        self.plot = False

    def fitalldata(self, head, dirname, fsts, outf):
        khead = head + "pickles/" + dirname
        hhead = head + "srlz/" + dirname
        outfile = khead + "kh_fit_results.txt"
        self.data = np.zeros((len(fsts), 5))
        for ifx, fs in enumerate(fsts):
            #print("fracevents: ", fs)
            if fs == "1":
                sfs = "_"
            else:
                values = fs.split(".")
                sfs = "_" + values[0] + values[1] + "_"
            self.kdevf = khead +\
                "qens_run6204united_kde_results_on_sdata_qsel.pkl"
            self.ktf = khead + "qens_run6202united" + sfs +\
                "kde_results_on_sdata_qsel.pkl"
            self.hdevf = hhead + "run6204united_sspectra.pkl"
            self.htf = hhead + "run6202united" + sfs + "sspectra.pkl"
            self.data[ifx, 0] = float(fs)
            self.quiet = True
            self.data[ifx, 1:] = self.kde_hist()
        np.savetxt(outfile, self.data, fmt='%.5e')

    def plotter(self, dataname=None, isend=True, twocomp=False,
                ymin=None, ymax=None,
                c=['pink', 'red', 'lightblue', 'blue'],
                loc=None, bbox_to_anchor=None):
        if dataname is None:
            if twocomp:
                plt.errorbar(self.data[:, 0], self.data[:, 1],
                             yerr=self.data[:, 2], marker="o", label='kde1',
                             ms=2, elinewidth=1, lw=0, capsize=3, c=c[0])
                plt.errorbar(self.data[:, 0], self.data[:, 3],
                             yerr=self.data[:, 4], marker="o", label='kde2',
                             ms=2, elinewidth=1, lw=0, capsize=3, c=c[1])
                plt.errorbar(self.data[:, 0], self.data[:, 5],
                             yerr=self.data[:, 6], marker="o", label='hist1',
                             ms=2, elinewidth=1, lw=0, capsize=3, c=c[2])
                plt.errorbar(self.data[:, 0], self.data[:, 7],
                             yerr=self.data[:, 8], marker="o", label='hist2',
                             ms=2, elinewidth=1, lw=0, capsize=3, c=c[3])
                if ymin is not None and ymax is not None:
                    plt.ylim([ymin, ymax])
            else:
                plt.errorbar(self.data[:, 0], self.data[:, 1],
                             yerr=self.data[:, 2], marker="o", label='kde',
                             ms=2, elinewidth=1, lw=0, capsize=3)
                plt.errorbar(self.data[:, 0], self.data[:, 3],
                             yerr=self.data[:, 4], marker="o", label='hist',
                             ms=2, elinewidth=1, lw=0, capsize=3)
                if ymin is not None and ymax is not None:
                    plt.ylim([ymin, ymax])
        else:
            plt.errorbar(self.data[:, 0], self.data[:, 1],
                         yerr=self.data[:, 2], marker="o", label=dataname,
                         ms=2, elinewidth=1, lw=0, capsize=3)
        if isend:
            plt.xlabel('Fraction of event data used')
            plt.ylabel('Gamma (meV)')
            plt.tick_params(direction='in', right=True, top=True)
            if loc is not None and bbox_to_anchor is not None:
                plt.legend(loc=loc, bbox_to_anchor=bbox_to_anchor)
            plt.show()

def samplerun():
    head = "/home/kazu/desktop/210108/Tatsumi/"
    dirname = "0000025new/"
    fsts = ["0.375", "0.5", "0.625", "0.75", "0.875", "1"]
    outf = "kh_fit_results.txt"
    prj = sqrun_kde_hist()
    prj.fitalldata(head, dirname, fsts, outf)
    prj.plotter()


#samplerun()
