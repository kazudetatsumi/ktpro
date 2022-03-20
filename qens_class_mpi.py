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
from mpi4py import MPI
sys.path.append("/home/kazu/desktop/210108/AdaptiveKDE/adaptivekde")
## if you use fotran library of ssvk or ssk, uncomment the corresponding lines
## and comment out the import lines of the  python library.
#import ssvkernel_fort as ssvkernel
#import sskernel_fort as sskernel 
import ssvkernel_mpi
import sskernel

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


class qens:
    def __init__(self, datadir, save_file, odata=True, qsel=False,
                 winparam=1, M=80, WinFunc='Boxcar',
                 figname='qens_out.png', showplot=True):
        self.datadir = datadir
        self.save_file = save_file
        self.odata = odata
        self.qsel = qsel               # in def self.select_energy()
        self.winparam = winparam
        self.M = M
        self.WinFunc = WinFunc
        self.figname = figname
        self.showplot = showplot
        with open(self.save_file, 'rb') as f:
            self.dataset = pickle.load(f, encoding='latin1')

    def select_spectra(self):
        spectra = self.dataset['spectra']
        if self.odata:  # case outgoing beam
            if self.qsel:   # case spectra are already integrated over a
                            # specific q range.
                self.selected_spectra = self.dataset['spectra']
                self.selected_energy = self.dataset['energy']
            else:       # case dataset are distributed over 2-D PSD elements.
                dp = self.dataset['detector_position']
                mask = np.where((dp[:, 0] >= 10) & (dp[:, 1] <= 65))[0]
                self.selected_spectra = np.sum(spectra[mask, 1, :], axis=0)
                self.selected_energy = spectra[0, 0, :]
        else:
            spectra[0, 0, :] = spectra[0, 0, :] - 2.085
            mergin = 0.001
            xlim = np.array([-0.10025 - mergin, 0.14975 + mergin])
            mask = np.where((spectra[0, 0, :] >= xlim[0]) &
                            (spectra[0, 0, :] <= xlim[1]))[0]
            self.selected_spectra = spectra[0, 1, mask]
            self.selected_energy = spectra[0, 0, mask]
        #print(self.selected_energy[np.argmax(self.selected_spectra)])
        self.de = self.selected_energy[1] - self.selected_energy[0]
        self.get_xvec()

    def get_xvec(self):
        # To keep the accuracy, kde is executed on the channel numbers in the
        # histogram data.
        # The actual energies were retrieved by "_real" variables.
        print("CHK", self.selected_spectra.shape[0])
        self.xvec = np.array([idx for idx in
                             range(0, self.selected_spectra.shape[0]) for
                             num_repeat in
                             range(0, int(self.selected_spectra[idx]))
                              ], dtype=float)
        self.xvec_real = np.array([self.selected_energy[idx] for idx in
                                  range(0, self.selected_spectra.shape[0]) for
                                  num_repeat in
                                  range(0, int(self.selected_spectra[idx]))
                                   ], dtype=float)

    def add_shift(self):
        np.random.seed(314)
        self.xvecorg = np.array(self.xvec)
        self.shift = np.random.uniform(-0.5, 0.5, size=self.xvec.shape[0])
        self.xvec += self.shift
        self.xvec_real += self.shift*self.de
        #print(self.xvec[0:30])

    def run_ssvkernel(self):
        self.tin = np.arange(self.selected_energy.shape[0])
        self.tin_real = np.linspace(self.selected_energy[0],
                                    self.selected_energy[-1], num=800)
                                    #self.selected_energy[-1], num=8000)
                                    #self.selected_energy[-1], num=66700)
                                    #self.selected_energy[-1], num=200000)
        comm = MPI.COMM_WORLD
        self.y = ssvkernel_mpi.ssvkernel(self.xvec_real, self.tin_real,
                                         M=self.M, winparam=self.winparam,
                                         WinFunc=self.WinFunc)
        self.y_ = sskernel.sskernel(comm, self.xvec_real, self.tin_real)

    def plotter(self):
        #scf = (np.min(self.xvec_real) - np.max(self.xvec_real)) /\
        #      (np.min(self.xvec) - np.max(self.xvec))

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
            #ax.plot(self.tin_real, self.y[0]/self.de, c='r', label='ssvkernel')
            #ax.plot(self.tin_real, self.y_[0]/self.de, c='k', label='sskernel')
            ax.plot(self.tin_real, self.y[0], c='r', label='ssvkernel')
            ax.plot(self.tin_real, self.y_[0], c='k', label='sskernel')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False, width=1.5)
        ax.set_ylabel('density')
        ax.set_xlabel('energy ($\mu eV$)')
        #ax.set_yticks([0, 35, 70])
        #ax.set_ylim([0, 70])
        ax.set_ylim(0., np.max(self.y_[0])*1.2)
        plt.legend()
        ax = fig.add_subplot(3, 1, 2)
        ax.set_ylabel('density')
        ax.bar(self.selected_energy, norms, width=self.de,
               label='expt data', bottom=0.0001)
        ax.plot(self.tin_real, self.y_[0], c='k', label='sskernel')
        ax.plot(self.tin_real, self.y[0], c='r', label='ssvkernel')
        ax.set_yscale('log')
        ax.set_ylim(0.0001, np.max(self.y_[0])/self.de)
        ax.set_ylim(0.0001, np.max(self.y_[0]))
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False, width=1.5)
        tmprange = ax.get_xlim()
        plt.legend()
        ax = fig.add_subplot(3, 1, 3)
        #ax.plot(self.tin_real, self.y[2]*scf, c='r', label='ssvkernel')
        #ax.plot(self.tin_real, np.zeros_like(self.tin_real)+self.y_[2]*scf,
        ax.plot(self.tin_real, self.y[2], c='r', label='ssvkernel')
        ax.plot(self.tin_real, np.zeros_like(self.tin_real)+self.y_[2],
                c='k', label='sskernel')
        ax.set_ylabel('band-width')
        ax.set_xlabel('energy (meV)')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=True, width=1.5)
        ax.set_xlim(tmprange)
        #ax.set_ylim(0, np.max(self.y[2]*scf))
        #ax.set_ylim(0, 5)
        #ax.set_yticks([0, 0.003, 0.006])
        plt.subplots_adjust(hspace=0.0)
        plt.legend()
        plt.savefig(self.figname)
        if self.showplot:
            plt.show()

    def save_output(self,  output_file):
        dataset = {}
        dataset['y_ssvk'] = self.y
        dataset['y_ssk'] = self.y_
        dataset['tin_real'] = self.tin_real
        dataset['xlim'] = np.array([np.min(self.xvec_real),
                                   np.max(self.xvec_real)])
        with open(output_file, 'wb') as f:
            pickle.dump(dataset, f, -1)

    def save_outputs(self,  output_file):
        dataset = {}
        dataset['ys_ssvk'] = self.ys
        dataset['tin_real'] = self.tin_real
        #dataset['xlim'] = np.array([np.min(self.xvec_real),
        #                           np.max(self.xvec_real)])
        with open(output_file, 'wb') as f:
            pickle.dump(dataset, f, -1)


def samplerun():
    datadir = "/home/kazu/desktop/210108/Tatsumi/run6202_containers/"
    #dataset_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_dir = datadir + "../"
    save_file = dataset_dir + "/qens.pkl"
    proj = qens(datadir, save_file)
    proj.select_spectra()
    proj.add_shift()
    proj.run_ssvkernel()


#samplerun()
