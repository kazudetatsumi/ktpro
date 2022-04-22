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
                            #variables=[2.18704786e-04, 1.67980295e-02,
                            #           4.92405238e-05, 1.88866588e-03,
                            #           1.21127501e-01, 5.02759930e-02])
                            variables=[0.59704786e-00, 2.67980295e-02,
                                       3.82405238e-01, 7.88866588e-03,
                                       0.21127501e+00, 1.82759930e-02])
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

    def kde_baloon(self, x, y):
        self.selected_energy = x
        self.selected_spectra = y
        self.de = self.selected_energy[1] - self.selected_energy[0]
        self.get_xvec()
        self.add_shift_de()
        self.hist()
        return self.baloon_estimator()

    def cycle(self):
        self.outall = np.zeros((self.numcycle, 6))
        for cyidx in range(0, self.numcycle):
            simt = np.zeros(self.ml.shape)
            simd = np.zeros(self.yd.shape)
            if self.rank == 0:
                simd, simt = self.generate_data()
            simt = MPI.COMM_WORLD.bcast(simt)
            simd = MPI.COMM_WORLD.bcast(simd)
            ## kde for dev func.
            self.kde(self.x, simd)
            self.dt = self.y[1][1]-self.y[1][0]
            simyd = self.y[0]
            self.kde(self.x, simt)
            simyt = self.y[0]
            ## baloon for dev func, with the same bandwidths as those of target
            #self.dt = self.y[1][1] - self.y[1][0]
            #simyd = self.kde_baloon(self.x, simd)
            simyd = simyd/np.sum(simyd)/self.dt
            simyt = simyt/np.sum(simyt)/self.dt*100.
            out = self.optimize(self.y[1], simyd, simyt,
                                #variables=[1.73704786e-05, 2.66580295e-02,
                                #           9.96405238e-06, 7.00766588e-03,
                                #           2.00077501e-01, 1.78759930e-01])
                                #variables=[0.59704786e-00, 2.67980295e-02,
                                #           3.82405238e-01, 7.88866588e-03,
                                #           0.21127501e+00, 1.82759930e-02])
                                variables=[5.4e+01, 2.65e-02,
                                           3.7e+01, 7.0e-03,
                                           1.8e+01, 1.4e+01])
            if out[0] < 0 and out[1] < 0:
                #print("negative-negative")
                out[0] = out[0]*(-1.)
                out[1] = out[1]*(-1.)
            if out[2] < 0 and out[3] < 0:
                #print("negative-negative")
                out[2] = out[2]*(-1.)
                out[3] = out[3]*(-1.)
            if out[1] < out[3]:
                #print("exchange")
                tmpout = out[1]
                tmpout2 = out[0]
                out[1] = out[3]
                out[3] = tmpout
                out[0] = out[2]
                out[2] = tmpout2
            #if self.rank == 0:
                #print(cyidx, out)
            self.outall[cyidx, :] = out

    def generate_data(self):
        return np.random.poisson(self.yd/np.sum(self.yd)*59146.*0.5)*1.,\
               np.random.poisson(self.ml/np.sum(self.ml)*18944.*0.5)*1.

    def run_ssvkernel(self):
        self.tin = np.arange(self.selected_energy.shape[0])
        self.tin_real = np.linspace(self.selected_energy[0],
                                    self.selected_energy[-1], num=480)
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
        #MPI.COMM_WORLD.barrier()
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
        #MPI.COMM_WORLD.barrier()
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
        #MPI.COMM_WORLD.barrier()
        libssk.ssk(
                   byref(optw),
                   c_int(xsize),
                   c_int(tinsize),
                   self.xvec_real,
                   self.tin_real,
                   yopt
                   )
        #MPI.COMM_WORLD.barrier()
        return yopt, self.tin_real, optw

    def baloon_estimator(self):
        y_hist_nz = self.y_hist[self.y_hist > 0]
        tin_nz = self.y[1][self.y_hist > 0]
        yv = np.zeros((self.y[1].shape[0]))
        for xchidx in range(self.y[1].shape[0]):
            yv[xchidx] = np.sum(y_hist_nz * self.dt * self.Gauss(self.y[1][xchidx]-tin_nz,
                                                            self.y[2][xchidx]))
        return yv * np.sum(self.y_hist) / np.sum(yv * self.dt)

    def Gauss(self, x, w):
        return 1. / (2. * np.pi)**2 / w * np.exp(-x**2 / 2. / w**2)

    def hist(self):
        thist = np.concatenate((self.y[1], (self.y[1][-1]+self.dt)[np.newaxis]))
        self.y_hist = np.histogram(self.xvec_real, thist-self.dt/2.)[0] / self.dt


def testrun():
    np.random.seed(314)
    np.set_printoptions(linewidth=120)
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    proj = runkdenoidata(devf, tf, elim, elimw, numcycle=3000)
    proj.get_xmlyd()
    proj.cycle()
    if proj.rank == 0:
        proj.output()


testrun()
