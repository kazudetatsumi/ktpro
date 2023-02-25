#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
#import os
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from mpi4py import MPI
from qens_kde_resampled import qens_kde_resampled as qkr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


class qens_balloon_resamples(qkr):
    def __init__(self, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1):
        self.runNos = runNos
        self.Nb = Nb
        self.gammas = np.zeros((Nb, 2))
        self.DefineFiles()
        self.elim = [-0.03, 0.07]
        self.rank = MPI.COMM_WORLD.Get_rank()
        self.leastsq = False

    def getrsspectra(self, rsfile, inb=0):
        qkr.__init__(self, pklfile=rsfile)
        return self.spectrab[inb, 0, :], self.spectrab[inb, 1, :]

    def CalcBandW(self, orgfile, inb=0):
        qkr.__init__(self, pklfile=orgfile)
        self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :], num=6400)
        #plt.plot(self.y[1], self.y[0])
        #plt.show()
        return self.y

    def balloon(self, ky, sy):
        bw = np.interp(sy[0], ky[1], ky[2])
        idx = sy[1].nonzero()
        sy_nz = sy[1][idx]
        t_nz = sy[0][idx]
        yv = np.zeros_like(sy[0])
        dt = min(np.diff(sy[0]))
        for k in range(yv.shape[0]):
            yv[k] = np.sum(sy_nz * dt * self.Gauss(sy[0][k]-t_nz, bw[k]))
        yv = yv * np.max(ky[0]) / np.max(yv)
        #plt.plot(sy[0], yv)
        #plt.yscale('log')
        #plt.show()
        return yv

    def Gauss(self, x, w):
        return 1.0/((2.0*np.pi)**0.5*w)*np.exp(-(x/w)**2/2.0)

    def eachrunno(self, fidx, inb):
        sy = self.getrsspectra(self.rsfiles[fidx], inb)
        syb = self.balloon(self.kys[fidx], sy)
        return sy[0], syb, sy[1]
        #return sy[0], sy[1], sy[1]

    def DefineFiles(self):
        self.rsfiles = []
        self.orgfiles = []
        for runno in self.runNos:
            self.rsfiles.append("./run" + str(runno) + "spectrab.pkl")
            self.orgfiles.append("./run" + str(runno) + "spectraorg.pkl")

    def DoQf(self, inb):
        xt, yt, yth = self.eachrunno(0, inb)
        xd, yd, ydh = self.eachrunno(1, inb)
        self.icorr()
        xtl, ytl = self.limit2(xt, yt, self.elim)
        xdl, ydl = self.limit2(xd, yd, self.elim)
        if inb == 0:
            print(np.sum(xtl - xdl), 'check x_tf - x_df')
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        self.bg = 0.
        #self.optimize(xdl, ydlc, ytlc, variables=[0.655, 0.0129, 0.200, 0.00208])
        self.check_out(inb, self.optimize(xdl, ydlc, ytlc,
                       variables=[0.655, 0.0129, 0.200, 0.00208]))

    def run(self):
        self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles]
        self.outall = []
        for inb in range(self.Nb):
            self.DoQf(inb)
        self.output()



def testrun():
    Nb = 48
    elim = [-0.03, 0.07]
    prj = qens_balloon_resamples(runNos=[6202, 6204], elim=elim, Nb=Nb)
    #print(qens_balloon_resamples.__mro__)
    prj.run()


testrun()
