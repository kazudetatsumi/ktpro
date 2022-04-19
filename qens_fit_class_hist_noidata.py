#!/usr/bin/env python
import numpy as np
import scipy.optimize as so
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class import qens_fit as qf


class runhistnoidata(qf):
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
        self.ml = self.reconstruct(x, yd, out)
        self.yd = yd
        self.x = x

    def preprocess(self):
        self.icorr()
        xd, yd = self.get_data(self.devf)
        xt, yt = self.get_data(self.tf)
        xdl, ydl = self.limit(xd, yd, self.elimw)
        xtl, ytl = self.limit(xt, yt, self.elimw)
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        return xtl, ydlc, ytlc

    def cycle(self):
        self.outall = np.zeros((self.numcycle, 6))
        for cyidx in range(0, self.numcycle):
            simd, simt = self.generate_data()
            out = self.optimize(self.x, simd, simt,
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


def testrun():
    np.set_printoptions(linewidth=120)
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.06, 0.10]
    proj = runhistnoidata(devf, tf, elim, elimw, numcycle=999)
    proj.get_xmlyd()
    proj.cycle()


#testrun()
