#!/usr/bin/env python
import numpy as np
import os
import sys
from ctypes import *
from mpi4py import MPI
sys.path.append("/home/kazu/ktpro")
from qens_class_fort_mpi import qens as qc
#from qens_class import qens as qc
from qens_fit_class_hist_noidata  import runhistnoidata as rhn


class runkdenoidata(rhn, qc):
    def __init__(self, devf, tf, outfile, alpha, elim, elimw, numcycle=100,
                 leastsq=True):
        self.elim = elim
        self.elimw = elimw
        self.devf = devf
        self.tf = tf
        self.outfile = outfile
        self.alpha = alpha
        self.elim = elim
        self.numcycle = numcycle
        self.leastsq = leastsq
        self.rank = MPI.COMM_WORLD.Get_rank()

    def get_xmlyd(self):
        x, yd, yt = self.preprocess()
        _out = self.optimize(x, yd, yt,
                             # variables=[2.18704786e-04, 1.67980295e-02,
                             #           4.92405238e-05, 1.88866588e-03,
                             #           1.21127501e-01, 5.02759930e-02])
                             variables=[0.59704786e-00, 2.67980295e-02,
                                        3.82405238e-01, 7.88866588e-03,
                                        0.21127501e+00, 1.82759930e-02])
        if self.rank == 0:
            print(_out[0])
        # self.ml = self.reconstruct(self.x, self.yd, out)
        self.ml = self.reconstruct(x, yd, _out[0])
        self.yd = yd
        self.x = x
        # print(self.x[0], self.x[-1])

    def preprocess(self):
        self.icorr()
        xd, self.yd = self.get_data(self.devf)
        self.x, yt = self.get_data(self.tf)
        xdl, ydl = self.limit2(xd, self.yd, self.elimw)
        xtl, ytl = self.limit2(self.x, yt, self.elimw)
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        return xtl, ydlc, ytlc

    def kde(self, x, y, M=160, winparam=1, num=800, WinFunc='Boxcar',
            isnovariablebw=False):
        self.WinFunc = WinFunc
        self.M = M
        self.winparam = winparam
        self.selected_spectra = y
        self.selected_energy = x
        self.de = self.selected_energy[1] - self.selected_energy[0]
        self.get_xvec()
        self.add_shift_de()
        self.run_ssvkernel(num=num, isnovariablebw=isnovariablebw)

    def kde_baloon(self, x, y):
        self.selected_energy = x
        self.selected_spectra = y
        self.de = self.selected_energy[1] - self.selected_energy[0]
        self.get_xvec()
        self.add_shift_de()
        self.hist()
        return self.baloon_estimator()

    def cycle(self):
        self.outall = []
        for cyidx in range(0, self.numcycle):
            simt = np.zeros(self.ml.shape)
            simd = np.zeros(self.yd.shape)
            if self.rank == 0:
                simd, simt = self.generate_data()
            simt = MPI.COMM_WORLD.bcast(simt)
            simd = MPI.COMM_WORLD.bcast(simd)
            # kde for dev func.
            self.kde(self.x, simd)
            self.dt = self.y[1][1]-self.y[1][0]
            simyd = self.y[0]
            self.kde(self.x, simt)
            simyt = self.y[0]
            # baloon for dev func, with the same bandwidths as those of target
            # self.dt = self.y[1][1] - self.y[1][0]
            # simyd = self.kde_baloon(self.x, simd)
            simyd = simyd/np.sum(simyd)/self.dt
            simyt = simyt/np.sum(simyt)/self.dt*100.
            _out = self.optimize(self.y[1], simyd, simyt,
                                 # variables=[1.73704786e-05, 2.66580295e-02,
                                 #            9.96405238e-06, 7.00766588e-03,
                                 #            2.00077501e-01, 1.78759930e-01])
                                 # variables=[0.59704786e-00, 2.67980295e-02,
                                 #            3.82405238e-01, 7.88866588e-03,
                                 #            0.21127501e+00, 1.82759930e-02])
                                 variables=[5.4e+01, 2.65e-02,
                                            3.7e+01, 7.0e-03,
                                            1.8e+01, 1.4e+01])
            self.check_out(cyidx, _out)

    def run_ssvkernel_notused(self):
        self.tin = np.arange(self.selected_energy.shape[0])
        self.tin_real = np.linspace(self.selected_energy[0],
                                    #self.selected_energy[-1], num=480)
                                    self.selected_energy[-1], num=8000)
        #print('number of tin_real elements=', self.tin_real.shape[0])

        if self.WinFunc == 'Boxcar':
            WinFuncNo = 1
        if self.WinFunc == 'Gauss':
            WinFuncNo = 2
        if self.WinFunc == 'Cauchy':
            WinFuncNo = 3

        self.y = self.calc_ssvkernel_f90(WinFuncNo)

    def baloon_estimator(self):
        y_hist_nz = self.y_hist[self.y_hist > 0]
        tin_nz = self.y[1][self.y_hist > 0]
        yv = np.zeros((self.y[1].shape[0]))
        for xchidx in range(self.y[1].shape[0]):
            yv[xchidx] = np.sum(y_hist_nz * self.dt *
                                self.Gauss(self.y[1][xchidx]-tin_nz,
                                           self.y[2][xchidx]))
        return yv * np.sum(self.y_hist) / np.sum(yv * self.dt)

    def Gauss(self, x, w):
        return 1. / (2. * np.pi)**2 / w * np.exp(-x**2 / 2. / w**2)

    def hist(self):
        thist = np.concatenate((self.y[1], (self.y[1][-1]+self.dt)[np.newaxis])
                               )
        self.y_hist = np.histogram(self.xvec_real,
                                   thist-self.dt/2.)[0] / self.dt


def testrun():
    outfile = './outkde.pkl'
    alpha = 0.5
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 10
    binwidth1 = 0.0016
    binwidth2 = 0.0016
    if os.path.isfile(outfile):
        proj = runkdenoidata(devf, tf, outfile, alpha, elim, elimw,
                             numcycle=numcycle, leastsq=True)
        proj.loadfile()
        proj.output()
        proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = runkdenoidata(devf, tf, outfile, alpha, elim, elimw,
                             numcycle=numcycle, leastsq=True)
        proj.get_xmlyd()
        proj.cycle()
        if proj.rank == 0:
            proj.output()
            proj.savefile()


#testrun()
