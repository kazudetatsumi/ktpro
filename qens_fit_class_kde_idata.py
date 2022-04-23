#!/usr/bin/env python
import numpy as np
import os
import sys
from mpi4py import MPI
sys.path.append("/home/kazu/ktpro")
from qens_class_fort_mpi import qens as qc
from qens_fit_class_kde import runkdenoidata as rkn


class runkdeidata(rkn, qc):
    def __init__(self, devf, tf, idevf, itf, outfile, alpha, elim, elimw,
                 numcycle=100, leastsq=True):
        self.elim = elim
        self.elimw = elimw
        self.devf = devf
        self.tf = tf
        self.idevf = idevf
        self.itf = itf
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
                             #            4.92405238e-05, 1.88866588e-03,
                             #            1.21127501e-01, 5.02759930e-02])
                             variables=[0.59704786e-00, 2.67980295e-02,
                                        3.82405238e-01, 7.88866588e-03,
                                        0.21127501e+00, 1.82759930e-02])
        if self.rank == 0:
            print(_out[0])
        self.ml = self.reconstruct(x, yd, _out[0])
        self.yd, self.ml = self.decorrection(x, yd, self.ml)
        self.ml = self.ml*np.interp(x, self.ixt, self.iyt)
        self.yd = self.yd*np.interp(x, self.ixd, self.iyd)
        self.x = x

    def preprocess(self):
        self.icorr()
        xd, self.yd = self.get_data(self.devf)
        self.x, yt = self.get_data(self.tf)
        self.ixd, self.iyd = self.get_idata(self.idevf)
        self.ixt, self.iyt = self.get_idata(self.itf)
        xdl, ydl = self.limit(xd, self.yd, self.elimw)
        xtl, ytl = self.limit(self.x, yt, self.elimw)
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        return xtl, ydlc, ytlc

    def cycle(self):
        self.outall = []
        for cyidx in range(0, self.numcycle):
            simt = np.zeros(self.ml.shape)
            simd = np.zeros(self.yd.shape)
            if self.rank == 0:
                simd, simt = self.generate_data()
            simt = MPI.COMM_WORLD.bcast(simt)
            simd = MPI.COMM_WORLD.bcast(simd)
            self.kde(self.x, simd)
            self.dt = self.y[1][1]-self.y[1][0]
            if cyidx == 0:
                iyd = np.interp(self.y[1], self.ixd, self.iyd)
            simyd = self.y[0]
            simydi = simyd/iyd
            self.kde(self.x, simt)
            if cyidx == 0:
                iyt = np.interp(self.y[1], self.ixt, self.iyt)
            simyt = self.y[0]
            simyti = simyt/iyt
            simydc, simytc = self.correction(self.y[1], simydi, simyti)
            # check spectral shapes
            # plt.plot(self.y[1], simyd/10.0, label='raw')
            # plt.plot(self.y[1], simydi, label='divi')
            # plt.plot(self.y[1], simydc, label='divi+correction')
            # plt.yscale('log')
            # plt.legend()
            # plt.show()
            simydc = simydc/np.sum(simydc)/self.dt
            simytc = simytc/np.sum(simytc)/self.dt*100.
            _out = self.optimize(self.y[1], simydc, simytc,
                                 # variables=[1.73704786e-05, 2.66580295e-02,
                                 #            9.96405238e-06, 7.00766588e-03,
                                 #            2.00077501e-01, 1.78759930e-01])
                                 # variables=[0.59704786e-00, 2.67980295e-02,
                                 #           3.82405238e-01, 7.88866588e-03,
                                 #            0.21127501e+00, 1.82759930e-02])
                                 variables=[4.8e+01, 2.65e-02,
                                            3.2e+01, 7.0e-03,
                                            1.9e+01, 1.1e+01])
            self.check_out(cyidx, _out)


def testrun():
    np.set_printoptions(linewidth=120)
    outfile = "./outkdeidata.pkl"
    alpha = 0.5
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    prefix = "/home/kazu/desktop/210108/Tatsumi/pickles/0000001io/"
    idevf = prefix + "qens_kde_results_on_idata_6204.pkl"
    itf = prefix + "qens_kde_results_on_idata_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 10
    binwidth1 = 0.0016
    binwidth2 = 0.0016
    if os.path.isfile(outfile):
        proj = runkdeidata(devf, tf, idevf, itf, outfile, alpha, elim, elimw,
                           leastsq=True, numcycle=numcycle)
        if proj.rank == 0:
            proj.loadfile()
            proj.output()
            proj.plot_distribution(binwidth1, binwidth2)
    else:
        np.random.seed(314)
        proj = runkdeidata(devf, tf, idevf, itf, outfile, alpha, elim, elimw,
                           leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        if proj.rank == 0:
            proj.output()
            proj.savefile()


#testrun()
