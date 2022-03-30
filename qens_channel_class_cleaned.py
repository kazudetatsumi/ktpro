# Adaptive kernel density estimation on qens data:
# Presently, this script can prepare the input for the adaptive kernel density
# estimation from txt files containing a histogram per a detector element.
# These txt files are to be deduced from the raw data by a script written by
# BL02 staffs and stored in a directory.
# This script can apply the adaptive kernel density estimation on the prepared
# data and show the plots of the results, i.e., the estimated density and the
# estimated destribution of the kernel widths.
# This script uses the python methods of the adaptive kenerl denssity
# estimation distributed at the site linked by Dr. Shimazaki's  HP.
# The adaptive kenel density estimation was developped by Dr. Shimazaki.
# Kazuyoshi TATSUMI 2021/04/03
import numpy as np
import os
import re
import pickle
import matplotlib.pylab as plt
import sys
sys.path.append("/home/kazu/desktop/210108/AdaptiveKDE/adaptivekde")
import ssvkernel
import sskernel
sys.path.append("/home/kazu/ktpro")
import qens_class as qcorg

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


class qens_channel(qcorg.qens):
    def __init__(self, datadir, save_file, odata=True, sfile=None, qsel=False,
                 optsm=False, winparam=1, M=80, WinFunc='Boxcar',
                 figname='qens_out.png', showplot=True):
        super().__init__(datadir, save_file, odata=odata, sfile=sfile, qsel=qsel,
                       optsm=optsm, winparam=winparam, M=M, WinFunc=WinFunc,
                       figname=figname, showplot=showplot)

    def run_ssvkernel(self):
        if self.optsm:
        ## tin and tin_real are modified now, so as to be used commonly for
        ## different data sets.
            #self.tin = np.linspace(0.0, 1000.0, 2001)
            tinmax = 10.**int(np.log10(self.selected_spectra.shape[0])+1.)
            self.tin = np.linspace(0.0, tinmax, int(tinmax)*2+1)
            print("Check parameters of horizontal axis")
            print("de=", self.de, "selected_energy[0]=",
                  self.selected_energy[0], "num channels=", self.tin.shape[0])
            self.tin_real = self.tin*self.de + self.selected_energy[0]
        else:
            #self.y = ssvkernel.ssvkernel(np.array(testx))
            T = (np.max(self.xvec) - np.min(self.xvec))
            dx = np.sort(np.diff(np.sort(self.xvec)))
            dt_samp = dx[np.nonzero(dx)][0]
            # To relieve the effect of the finite counts near the upper and lower
            # bounds of the measurement range
            # we extend the range by some fraction of the original range
            #mergin = T*0.2
            #tin = np.linspace(np.min(self.xvec) - mergin, np.max(self.xvec) + mergin,
            #                  int(min(np.ceil(T*1.4 / dt_samp), 1e3)))
            # tin is the space where densities are inferred.

            self.tin = np.linspace(np.min(self.xvec), np.max(self.xvec),
                              int(min(np.ceil(T / dt_samp), 1e3)))
            self.tin_real = np.linspace(np.min(self.xvec_real), np.max(self.xvec_real),
                                   int(min(np.ceil(T / dt_samp), 1e3)))
        self.y = ssvkernel.ssvkernel(self.xvec, self.tin, M=self.M,
                                     winparam=self.winparam,
                                     WinFunc=self.WinFunc)
        self.y_ = sskernel.sskernel(self.xvec, self.tin)

        scf = (np.min(self.xvec_real) - np.max(self.xvec_real)) /\
              (np.min(self.xvec) - np.max(self.xvec))

        if self.optsm:
            self.optsmear()
            snorms = self.sspectra/np.sum(self.sspectra)/self.sde

        norms = self.selected_spectra/np.sum(self.selected_spectra)/self.de
        fig = plt.figure(figsize=(6, 8))
        ax = fig.add_subplot(3, 1, 1)
        ax.bar(self.selected_energy, norms, width=self.de, label='expt data')
        if self.optsm:
            #ax.bar(self.senergy, snorms, width=self.sde, label='expt sdata')
            ax.plot(self.tin_real, self.yck/self.de, c='r', label='yck')
            ax.plot(self.tin_real, self.y[0]/self.de, c='k', label='ssvkernel')
        else:
            ax.plot(self.tin_real, self.y[0]/self.de, c='r', label='ssvkernel')
            ax.plot(self.tin_real, self.y_[0]/self.de, c='k', label='sskernel')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False, width=1.5)
        ax.set_ylabel('density')
        ax.set_xlabel('energy ($\mu eV$)')
        ax.set_yticks([0, 50, 100])
        plt.legend()
        ax = fig.add_subplot(3, 1, 2)
        ax.set_ylabel('density')
        if self.optsm:
            ax.bar(self.senergy, snorms, width=self.sde, label='expt sdata',
                   bottom=0.0001)
            ax.plot(self.tin_real, self.yck/self.de, c='r', label='yck')
            ax.plot(self.tin_real, self.y[0]/self.de, c='k', label='ssvkernel')
        else:
            ax.bar(self.selected_energy, norms, width=self.de,
                   label='expt data', bottom=0.0001)
            ax.plot(self.tin_real, self.y_[0]/self.de, c='k', label='sskernel')
            ax.plot(self.tin_real, self.y[0]/self.de, c='r', label='ssvkernel')
        ax.set_yscale('log')
        ax.set_ylim(0.0001, np.max(self.y_[0])/self.de)
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False, width=1.5)
        tmprange = ax.get_xlim()
        plt.legend()
        ax = fig.add_subplot(3, 1, 3)
        ax.plot(self.tin_real, self.y[2]*scf, c='r', label='ssvkernel')
        ax.plot(self.tin_real, np.zeros_like(self.tin_real)+self.y_[2]*scf,
                c='k', label='sskernel')
        ax.set_ylabel('band-width')
        ax.set_xlabel('energy (meV)')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=True, width=1.5)
        ax.set_xlim(tmprange)
        ax.set_ylim(0, np.max(self.y[2]*scf))
        ax.set_yticks([0, 0.003, 0.006])
        plt.subplots_adjust(hspace=0.0)
        plt.legend()
        plt.savefig(self.figname)
        if self.showplot:
            plt.show()

def samplerun():
    datadir = "/home/kazu/desktop/210108/Tatsumi/run6202_containers/"
    #dataset_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_dir = datadir + "../"
    save_file = dataset_dir + "/qens.pkl"
    proj = qens_channel(datadir, save_file)
    proj.select_spectra()
    proj.add_shift()
    proj.run_ssvkernel()


#samplerun()
