#!/usr/bin/env python
import numpy as np
import scipy.optimize as so
import scipy.signal as ss
import sys
import re
import pickle
from ctypes import *
from mpi4py import MPI
sys.path.append("/home/kazu/ktpro")
lib = CDLL("/home/kazu/ktpro/ssvkernel_f90_mpi.so")
libssk = CDLL("/home/kazu/ktpro/sskernel_f90.so")


class runkdenoidata():
    def __init__(self, devf, tf, elim, elimw, numcycle=100):
        self.elim = elim
        self.elimw = elimw
        self.devf = devf
        self.tf = tf
        self.elim = elim
        self.numcycle = numcycle

    def get_xmlyd(self):
        x, yd, yt = self.preprocess()
        out = self.optimize(x, yd, yt,
                            variables=[2.18704786e-04, 1.67980295e-02,
                                       4.92405238e-05, 1.88866588e-03,
                                       1.21127501e-01, 5.02759930e-02])
        self.ml = self.reconstruct(self.x, self.yd, out)

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
        self.selected_spectra = x
        self.selected_energy = y
        self.de = self.selected_energy[1] - self.selected_energy[0]
        self.get_xvec()
        self.add_shift_de()
        self.run_ssvkernel()

    def cycle(self):
        self.outall = np.zeros((self.numcycle, 6))
        for cyidx in range(0, self.numcycle):
            simd, simt = self.generate_data()
            self.kde(self.x, simd)
            simyd = self.y[0]
            self.kde(self.x, simt)
            simyt = self.y[0]
            out = self.optimize(self.x, simyd, simyt,
                                variables=[6.11704786e-06, 2.51980295e-02,
                                           1.55405238e-06, 4.28866588e-03,
                                           7.97127501e-03, 3.52759930e-01])
            self.outall[cyidx, :] = out
        print(np.average(self.outall[:, 1]), np.std(self.outall[:, 1]))
        print(np.average(self.outall[:, 3]), np.std(self.outall[:, 3]))
        mask = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                        & (self.outall[:, 2] > 0) & (self.outall[:, 3] > 0)
                        & (self.outall[:, 4] > 0) & (self.outall[:, 5] > 0))
        self.outnonneg = self.outall[mask]
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
        y = alpha1*self.convlore(d, gamma1, x)\
            + alpha2*self.convlore(d, gamma2, x)\
            + delta*d + base
        # A smaller energy range is set for the squre differences,
        # because y involves convolution and this setting is preferable
        # to decrease the edge effect.
        xl, dif = self.limit(x, t-y, self.elim)
        return dif

    def reconstruct(self, x, yd, out):
        _alpha, _gamma, _alpha2, _gamma2,  _delta, _base = out
        return self.convlore(_alpha*yd, _gamma, x)\
            + self.convlore(_alpha2*yd, _gamma2, x)\
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
                                    self.selected_energy[-1], num=8000)
        print('number of tin_real elements=', self.tin_real.shape[0])

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
        libssk.ssk(
                   byref(optw),
                   c_int(xsize),
                   c_int(tinsize),
                   self.xvec_real,
                   self.tin_real,
                   yopt
                   )
        return yopt, self.tin_real, optw

    def get_xvec(self):
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

    def add_shift_de(self):
        rank = MPI.COMM_WORLD.Get_rank()
        if rank == 0:
            self.shift = np.random.uniform(0., 1., size=self.xvec.shape[0])
            self.xvec_real += self.shift*self.de
        self.xvec_real = MPI.COMM_WORLD.bcast(self.xvec_real)

    def res_icorr(self, coeffs, x, t):
        [k0, k1, k2, k3] = coeffs
        y = k0 + k1*x + k2*x**2 + k3*x**3
        return t - y

    def get_icorrdata(self, icorrfile):
        x = []
        y = []
        for line in open(icorrfile):
            if not re.compile('[A-z[]').search(line):
                x.append(float(line.split()[0]))
                y.append(float(line.split()[1]))
        x = np.array(x)
        y = np.array(y)
        return x, y

    def icorr(self):
        x, y = self.get_icorrdata("/home/kazu/desktop/210108/Tatsumi/" +
                                  "Ilambda_correction.txt")
        variables = [10, 10, 10, 10]
        out = so.leastsq(self.res_icorr, variables,
                         args=(x, y), full_output=1,
                         epsfcn=0.0001)
        [_k0, _k1, _k2, _k3] = out[0]
        self.k = out[0]

    def get_data(self, infile):
        with open(infile, 'rb') as f:
            data = pickle.load(f, encoding='latin1')
        return data['xo'], data['yr_ssvk']

    def fun_lore(self, x, gamma):
        return 1/np.pi*gamma/(x**2 + gamma**2)

    def convlore(self, f, gamma, x):
        ex = np.linspace(np.min(x)-(np.max(x)-np.min(x))*0.5,
                         np.max(x)+(np.max(x)-np.min(x))*0.5, x.shape[0]*2+1)
        win = self.fun_lore(ex - ex[int(x.shape[0])], gamma)
        _convd = ss.convolve(f, win, mode='same', method='fft')
        return _convd



def testrun():
    np.set_printoptions(linewidth=120)
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.06, 0.10]
    proj = runkdenoidata(devf, tf, elim, elimw, numcycle=3)
    proj.get_xmlyd()
    proj.cycle()


testrun()
