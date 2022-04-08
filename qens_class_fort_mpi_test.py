# Adaptive kernel density estimation on qens data:
# Presently, this script can prepare the input for the adaptive kernel density
# estimation from txt files containing a histogram per a detector element.
# These txt files are to be deduced from the raw data by a script written by
# BL02 staffs and stored in a directory.
# This script can apply the adaptive kernel density estimation on the prepared
# data and show the plots of the results, i.e., the estimated density and the
# estimated distribution of the kernel widths.
# This script uses the python methods of the adaptive kenerl density
# estimation distributed at the site linked by Dr. Shimazaki's  HP.
# The adaptive kenel density estimation was developed by Dr. Shimazaki.
# Kazuyoshi TATSUMI 2021/04/03
# MPI parallelized costfunction, which is very messy, but works.
# Kazuyoshi TATSUMI 2022/03/17
import numpy as np
import os
import re
from ctypes import *
import pickle
import matplotlib
matplotlib.use('agg')
import matplotlib.pylab as plt
from mpi4py import MPI
import sys
sys.path.append("/home/kazu/desktop/210108/AdaptiveKDE/adaptivekde")
## ssvkernel compatibility between python and fortran versions is now destroyed.
## This class has an alternative method using mpi.
lib = CDLL("/home/kazu/ktpro/ssvkernel_f90_mpi_test.so")
## Either of sskernel_fort (fortran ver.) or sskernel (python ver.) can be set by
## uncommenting the corresponding line below.
import sskernel_fort as sskernel 
#import sskernel

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
        #print("CHK", self.selected_spectra.shape[0])
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
        self.xvecorg = np.array(self.xvec)
        self.shift = np.random.uniform(-0.5, 0.5, size=self.xvec.shape[0])
        self.xvec += self.shift
        self.xvec_real += self.shift*self.de
        #print(self.xvec[0:30])

    def add_shift_de(self):
        self.xvecorg = np.array(self.xvec)
        self.shift = np.random.uniform(0., 1., size=self.xvec.shape[0])
        self.xvec += self.shift
        self.xvec_real += self.shift*self.de
        #print(self.xvec[0:30])

    def run_ssvkernel(self):
        self.tin = np.arange(self.selected_energy.shape[0])
        #print("Check parameters of horizontal axis")
        #print("de=", self.de, "selected_energy[0]=",
        #      self.selected_energy[0], "num channels=", self.tin.shape[0])
        self.tin_real = np.linspace(self.selected_energy[0],
                                    #self.selected_energy[-1],
                                    #num=self.selected_spectra.shape[0])
                                    #self.selected_energy[-1], num=8000)
                                    #self.selected_energy[-1], num=80000)
                                    #self.selected_energy[-1], num=66700)
                                    #self.selected_energy[-1], num=200000)
                                    self.selected_energy[-1], num=2000000)
        #print(self.tin_real[0:10])
        print('number of tin_real elements=', self.tin_real.shape[0])

        if self.WinFunc=='Boxcar':
            WinFuncNo=1
        if self.WinFunc=='Gauss':
            WinFuncNo=2
        if self.WinFunc=='Cauchy':
            WinFuncNo=3

        self.y = self.calc_ssvkernel_f90(WinFuncNo)
        self.y_ = sskernel.sskernel(self.xvec_real, self.tin_real)

    def calc_ssvkernel_f90(self, WinFuncNo):
        lib.ssvk.restype = c_void_p
        lib.ssvk.argtypes = [ 
                            POINTER(c_int32),
                            POINTER(c_int),
                            POINTER(c_double),
                            POINTER(c_int),
                            POINTER(c_int),
                            POINTER(c_int),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)
                            ]
        xsize = self.xvec_real.shape[0]
        tinsize = self.tin_real.shape[0]
        yopt = np.zeros((tinsize))
        optw = np.zeros((tinsize))
        comm = MPI.COMM_WORLD
        comm = comm.py2f()

        lib.ssvk(
                c_int32(comm),
                c_int(self.M),
                c_double(self.winparam),
                c_int(xsize),
                c_int(tinsize),
                c_int(WinFuncNo),
                self.xvec_real,
                self.tin_real,
                optw,
                yopt
                )
        return yopt, self.tin_real, optw

    def plotter(self):
        #norms = self.selected_spectra/np.sum(self.selected_spectra)/self.de
        rank = MPI.COMM_WORLD.Get_rank()
        size = MPI.COMM_WORLD.Get_size()
        if rank == 0:
            fig = plt.figure(figsize=(6, 8))
            ax = fig.add_subplot(3, 1, 1)
            #ax.bar(self.selected_energy, norms, width=self.de, label='expt data')
            ax.plot(self.tin_real, self.y[0], c='r', label='ssvkernel')
            ax.plot(self.tin_real, self.y_[0], c='k', label='sskernel')
            ax.tick_params(top=True, right=True, direction='in', which='both',
                           labelbottom=False, width=1.5)
            ax.set_ylabel('density')
            ax.set_xlabel('energy ($\mu eV$)')
            ax.set_ylim(0., np.max(self.y_[0])*1.2)
            plt.legend()
            ax = fig.add_subplot(3, 1, 2)
            ax.set_ylabel('density')
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
            ax.plot(self.tin_real, self.y[2], c='r', label='ssvkernel')
            ax.plot(self.tin_real, np.zeros_like(self.tin_real)+self.y_[2],
                    c='k', label='sskernel')
            ax.set_ylabel('band-width')
            ax.set_xlabel('energy (meV)')
            ax.tick_params(top=True, right=True, direction='in', which='both',
                           labelbottom=True, width=1.5)
            ax.set_xlim(tmprange)
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
