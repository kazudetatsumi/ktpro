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
        # normalization
        dx = x[1] - x[0]
        yd = yd / np.sum(yd) / dx
        yt = yt / np.sum(yt) / dx
        out = self.optimize(x, yd, yt,
#                            variables=[2.18704786e-04, 1.67980295e-02,
#                                       4.92405238e-05, 1.88866588e-03,
#                                       1.21127501e-01, 5.02759930e-02])
                            #variables=[0.59704786e-00, 2.67980295e-02,
                            #           3.82405238e-01, 7.88866588e-03,
                            #           0.21127501e+00, 1.82759930e-02])
                            variables=[0.54e-00, 2.7e-02,
                                       3.5e-01, 7.0e-03,
                                       0.18e+00, 1.9e-02])
        print(out)
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
            # normalization
            # dx = self.x[1] - self.x[0]
            # simd = simd / np.sum(simd) / dx
            # simt = simt / np.sum(simt) / dx * 100.
            # rebinning
            tin_real, simdr = self.rebin_generated_samples(self.x, simd, num=480)
            tin_real, simtr = self.rebin_generated_samples(self.x, simt, num=480)
            dx = tin_real[1] - tin_real[0]
            simdr = simdr / np.sum(simdr) / dx
            simtr = simtr / np.sum(simtr) / dx * 100.
            out = self.optimize(tin_real, simdr, simtr,
            # out = self.optimize(self.x, simd, simt,
#                                variables=[1.73704786e-05, 2.66580295e-02,
#                                           9.96405238e-06, 7.00766588e-03,
#                                           2.00077501e-01, 1.78759930e-01])
                                #variables=[0.18704786e-00, 2.67980295e-02,
                                #           1.02405238e-01, 6.48866588e-03,
                                #           0.06127501e+00, 0.13759930e-00])
                                #variables=[0.7e-00, 2.6e-02,
                                #           2.3e-01, 7.0e-03,
                                #           0.17e+00, 0.38e-00])
                                variables=[6.2e+01, 2.7e-02,
                                           2.5e+01, 7.0e-03,
                                           2.0e+01, 4.0e+01])
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
            print(cyidx, out)
            self.outall[cyidx, :] = out
        print(np.average(self.outall[:, 1]), np.std(self.outall[:, 1]))
        print(np.average(self.outall[:, 3]), np.std(self.outall[:, 3]))
        mask = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                        & (self.outall[:, 2] > 0) & (self.outall[:, 3] > 0)
                        & (self.outall[:, 4] > 0) & (self.outall[:, 5] > 0))
        self.outnonneg = self.outall[mask]
        maskwobg = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                        & (self.outall[:, 2] > 0) & (self.outall[:, 3] > 0)
                        & (self.outall[:, 4] > 0))
        self.outnonnegwobg = self.outall[maskwobg]
        print(np.average(self.outnonneg[:, 1]), "+/-",
              np.std(self.outnonneg[:, 1]))
        print(np.average(self.outnonneg[:, 3]), "+/-",
              np.std(self.outnonneg[:, 3]))
        print(self.outnonneg.shape[0], "/", self.numcycle)
        print(np.average(self.outnonnegwobg[:, 1]), "+/-",
              np.std(self.outnonnegwobg[:, 1]))
        print(np.average(self.outnonnegwobg[:, 3]), "+/-",
              np.std(self.outnonnegwobg[:, 3]))
        print(self.outnonnegwobg.shape[0], "/", self.numcycle)


    def correction(self, x, yd, yt):
        x = x + 2.085
        yd = yd / (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        yt = yt / (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        return yd, yt

    def optimize(self, x, yd, yt,
                 variables=[6.e-6, 2.e-2, 1.e-6, 4.e-3, 7.e-3, 3.e-1]):
        # leastsq
        # out = so.leastsq(self.res, variables, args=(x, yd, yt), full_output=1,
        #                  epsfcn=0.0001)
        # return out[0]
        # least_squares
        bounds = (0, np.inf)
        #out = so.least_squares(self.res, variables, bounds=bounds, args=(x, yd, yt))
        out = so.least_squares(self.res, variables,  args=(x, yd, yt))
        print("status:", out.status)
        return out.x

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
        #_base = _base * 9.
        return _alpha*self.convlore(yd, _gamma, x)\
            + _alpha2*self.convlore(yd, _gamma2, x)\
            + _delta*yd + _base

    def limit(self, x, y, elim):
        mask = np.where((x > elim[0]) & (x < elim[1]))
        return x[mask], y[mask]

    def generate_data(self):
        return np.random.poisson(self.yd/np.sum(self.yd)*59146.*0.5)*1.,\
               np.random.poisson(self.ml/np.sum(self.ml)*18944.*0.5)*1.


def testrun():
    np.random.seed(314)
    np.set_printoptions(linewidth=120)
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    proj = runhistnoidata(devf, tf, elim, elimw, numcycle=300)
    proj.get_xmlyd()
    proj.cycle()


#testrun()
