#!/usr/bin/env python
import numpy as np
import scipy.optimize as so
import scipy.signal as ss
import sys
import re
import pickle
import matplotlib.pyplot as plt
from ctypes import *
from mpi4py import MPI
sys.path.append("/home/kazu/ktpro")
from qens_class_fort_mpi import qens as qc
from qens_fit_class_hist import runhist as rh
lib = CDLL("/home/kazu/ktpro/ssvkernel_f90_mpi.so")
libssk = CDLL("/home/kazu/ktpro/sskernel_f90.so")


class runkdenoidata(rh, qc):
    def __init__(self, devf, tf, elim, elimw, numcycle=100):
        self.elim = elim
        self.elimw = elimw
        self.devf = devf
        self.tf = tf
        self.elim = elim
        self.numcycle = numcycle
        self.rank = MPI.COMM_WORLD.Get_rank()

    def get_xmlyd(self):
        x, yd, yt = self.preprocess()
        out = self.optimize(x, yd, yt,
                            variables=[2.18704786e-04, 1.67980295e-02,
                                       4.92405238e-05, 1.88866588e-03,
                                       1.21127501e-01, 5.02759930e-02])
        if self.rank == 0:
            print(out)
        #self.ml = self.reconstruct(self.x, self.yd, out)
        self.ml = self.reconstruct(x, yd, out)
        self.yd = yd
        self.x = x
        #print(self.x[0], self.x[-1])

    def preprocess(self):
        self.icorr()
        xd, self.yd = self.get_data(self.devf)
        self.x, yt = self.get_data(self.tf)
        xdl, ydl = self.limit(xd, self.yd, self.elimw)
        xtl, ytl = self.limit(self.x, yt, self.elimw)
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        return xtl, ydlc, ytlc

    def kde(self, x, y):
        self.WinFunc = 'Boxcar'
        self.M = 160
        self.winparam = 1
        self.selected_spectra = y
        self.selected_energy = x
        self.de = self.selected_energy[1] - self.selected_energy[0]
        self.get_xvec()
        self.add_shift_de()
        self.run_ssvkernel()

    def cycle(self):
        self.outall = np.zeros((self.numcycle, 6))
        for cyidx in range(0, self.numcycle):
            simt = np.zeros(self.ml.shape)
            simd = np.zeros(self.yd.shape)
            if self.rank == 0:
                simd, simt = self.generate_data()
            simt = MPI.COMM_WORLD.bcast(simt)
            simd = MPI.COMM_WORLD.bcast(simd)
            self.kde(self.x, simd)
            simyd = self.y[0]
            #plt.plot(self.y[1], self.y[0])
            self.kde(self.x, simt)
            simyt = self.y[0]
            #plt.plot(self.y[1], self.y[0])
            #plt.yscale('log')
            #plt.show()
            out = self.optimize(self.y[1], simyd, simyt,
                                variables=[1.73704786e-05, 2.66580295e-02,
                                           9.96405238e-06, 7.00766588e-03,
                                           2.00077501e-01, 1.78759930e-01])
            if out[0] < 0 and out[1] < 0:
                print("negative-negative")
                out[0] = out[0]*(-1.)
                out[1] = out[1]*(-1.)
            if out[2] < 0 and out[3] < 0:
                print("negative-negative")
                out[2] = out[2]*(-1.)
                out[3] = out[3]*(-1.)
            if out[1] < out[3]:
                print("exchange")
                tmpout = out[1]
                tmpout2 = out[0]
                out[1] = out[3]
                out[3] = tmpout
                out[0] = out[2]
                out[2] = tmpout2
            if self.rank == 0:
                print(cyidx, out)
            self.outall[cyidx, :] = out
        if self.rank == 0:
            print(np.average(self.outall[:, 1]), np.std(self.outall[:, 1]))
            print(np.average(self.outall[:, 3]), np.std(self.outall[:, 3]))
        mask = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                        & (self.outall[:, 2] > 0) & (self.outall[:, 3] > 0)
                        & (self.outall[:, 4] > 0) & (self.outall[:, 5] > 0))
        self.outnonneg = self.outall[mask]
        if self.rank == 0:
            print(np.average(self.outnonneg[:, 1]), "+/-",
                  np.std(self.outnonneg[:, 1]))
            print(np.average(self.outnonneg[:, 3]), "+/-",
                  np.std(self.outnonneg[:, 3]))
            print(self.outnonneg.shape[0], "/", self.numcycle)

    def correction(self, x, yd, yt):
        x = x + 2.085
        yd = yd / (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        yt = yt / (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        return yd, yt

    def optimize(self, x, yd, yt,
                 variables=[6.e-6, 2.e-2, 1.e-6, 4.e-3, 7.e-3, 3.e-1]):
        out = so.leastsq(self.res, variables, args=(x, yd, yt), full_output=1,
                         epsfcn=0.0001)
        return out[0]

    def res(self, coeffs, x, d, t):
        [alpha1, gamma1, alpha2, gamma2,  delta, base] = coeffs
        y = alpha1*self.convloreorg(d, gamma1, x)\
            + alpha2*self.convloreorg(d, gamma2, x)\
            + delta*d + base
        # A smaller energy range is set for the squre differences,
        # because y involves convolution and this setting is preferable
        # to decrease the edge effect.
        xl, dif = self.limit(x, t-y, self.elim)
        return dif

    def reconstruct(self, x, yd, out):
        _alpha, _gamma, _alpha2, _gamma2,  _delta, _base = out
        return _alpha*self.convloreorg(yd, _gamma, x)\
            + _alpha2*self.convloreorg(yd, _gamma2, x)\
            + _delta*yd + _base

    def limit(self, x, y, elim):
        mask = np.where((x > elim[0]) & (x < elim[1]))
        return x[mask], y[mask]

    def generate_data(self):
        return np.random.poisson(self.yd/np.sum(self.yd)*59146.)*1.,\
               np.random.poisson(self.ml/np.sum(self.ml)*18944.)*1.

    def run_ssvkernel(self):
        self.tin = np.arange(self.selected_energy.shape[0])
        self.tin_real = np.linspace(self.selected_energy[0],
                                    self.selected_energy[-1], num=4800)
        #print('number of tin_real elements=', self.tin_real.shape[0])

        if self.WinFunc=='Boxcar':
            WinFuncNo=1
        if self.WinFunc=='Gauss':
            WinFuncNo=2
        if self.WinFunc=='Cauchy':
            WinFuncNo=3

        self.y = self.calc_ssvkernel_f90(WinFuncNo)
        self.y_ = self.calc_sskernel_f90()

    def calc_ssvkernel_f90(self, WinFuncNo):
        lib.ssvk.restype = c_void_p
        lib.ssvk.argtypes = [
                            POINTER(c_int32),
                            POINTER(c_int),
                            POINTER(c_double),
                            POINTER(c_int),
                            POINTER(c_int),
                            POINTER(c_int),
                            POINTER(c_int),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2)
                            ]
        xsize = self.xvec_real.shape[0]
        tinsize = self.tin_real.shape[0]
        yopt = np.zeros((tinsize))
        optw = np.zeros((tinsize))
        nb = 1
        yb = np.zeros((nb, tinsize))
        comm = MPI.COMM_WORLD
        comm = comm.py2f()
        MPI.COMM_WORLD.barrier()
        lib.ssvk(
                c_int32(comm),
                c_int(self.M),
                c_double(self.winparam),
                c_int(xsize),
                c_int(tinsize),
                c_int(WinFuncNo),
                c_int(nb),
                self.xvec_real,
                self.tin_real,
                optw,
                yopt,
                yb
                )
        MPI.COMM_WORLD.barrier()
        return yopt, self.tin_real, optw, yb

    def calc_sskernel_f90(self):
        libssk.ssk.restype = c_void_p
        libssk.ssk.argtypes = [
                            POINTER(c_double),
                            POINTER(c_int),
                            POINTER(c_int),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1),
                            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1)
                            ]
        xsize = self.xvec_real.shape[0]
        tinsize = self.tin_real.shape[0]
        yopt = np.zeros((tinsize))
        optw = c_double()
        MPI.COMM_WORLD.barrier()
        libssk.ssk(
                   byref(optw),
                   c_int(xsize),
                   c_int(tinsize),
                   self.xvec_real,
                   self.tin_real,
                   yopt
                   )
        MPI.COMM_WORLD.barrier()
        return yopt, self.tin_real, optw



def testrun():
    np.set_printoptions(linewidth=120)
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    proj = runkdenoidata(devf, tf, elim, elimw, numcycle=30)
    proj.get_xmlyd()
    proj.cycle()


testrun()
