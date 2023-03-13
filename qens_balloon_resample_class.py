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
#import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('TkAgg')


class qens_balloon_resamples(qkr):
    def __init__(self, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1,
                 ishist=False, num=6400, rsmodifier="b", orgmodifier="org",
                 prefix="./", variables=[0.655, 0.0129, 0.200, 0.00208],
                 quiet=False):
        self.runNos = runNos
        self.Nb = Nb
        self.gammas = np.zeros((Nb, 2))
        self.elim = [-0.03, 0.07]
        self.ishist = ishist
        self.num = num
        self.rsmodifier = rsmodifier
        self.orgmodifier = orgmodifier
        self.prefix = prefix
        self.variables = variables
        self.quiet = quiet
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.leastsq = False
        self.DefineFiles()

    def getrsspectra(self, rsfile, inb=0):
        qkr.__init__(self, pklfile=rsfile)
        return self.spectrab[inb, 0, :], self.spectrab[inb, 1, :]

    def CalcBandW(self, orgfile, inb=0):
        qkr.__init__(self, pklfile=orgfile)
        if self.ishist:
            if self.rank == 0 and not self.quiet:
                print("skipping KDE because ishist", self.ishist)
            self.y = "dummy"
        else:
            self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :],
                     num=self.num)
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

    def DefineFiles(self):
        self.rsfiles = []
        self.orgfiles = []
        for runno in self.runNos:
            self.rsfiles.append(self.prefix + "run" + str(runno) + "spectra" +
                                self.rsmodifier + ".pkl")
            self.orgfiles.append(self.prefix + "run" + str(runno) + "spectra" +
                                 self.orgmodifier + ".pkl")
        if self.rank == 0:
            if not self.quiet:
                print(self.prefix)
            print("")
            print("rsfiles:", [rsf.split(self.prefix)[1] for rsf in self.rsfiles])
            print("orgfiles:", [orgf.split(self.prefix)[1] for orgf in self.orgfiles])

    def eachrunno(self, fidx, inb):
        sy = self.getrsspectra(self.rsfiles[fidx], inb)
        if self.ishist:
            return sy[0], sy[1], sy[1]
        else:
            syb = self.balloon(self.kys[fidx], sy)
            return sy[0], syb, sy[1]

    def DoQf(self, inb):
        xt, yt, yth = self.eachrunno(0, inb)
        xd, yd, ydh = self.eachrunno(1, inb)
        self.icorr()
        xtl, ytl = self.limit2(xt, yt, self.elim)
        xdl, ydl = self.limit2(xd, yd, self.elim)
        if inb == 0 and self.rank == 0:
            if np.sum(xtl - xdl) > 0.000001:
                print('WARNING, check x_tf - x_df')
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        self.bg = 0.
        ####self.optimize(xdl, ydlc, ytlc, variables=[0.655, 0.0129, 0.200, 0.00208])
        self.check_out(inb, self.optimize(xdl, ydlc, ytlc,
        #               variables=[0.655, 0.0129, 0.200, 0.00208]))
                       #variables=[0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]))
                                          variables=self.variables))

    def run(self):
        self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles
                    ]
        self.outall = []
        #for inb in range(self.Nb):
        ##MPI Distributing the work to each rank
        ### in the case of Nb=1, we still want self.outall
        #if self.Nb == 1:
        #    if self.rank == 0:
        #        self.DoQf(0)
        #else:
        #    for inb in range(self.rank*self.Nb//self.size,
        #                     (self.rank+1)*self.Nb//self.size):
        #        self.DoQf(inb)
        #    ##MPI Gathering the result from each rank
        #    outallt = np.zeros((self.size*np.array(self.outall).size))
        #    self.comm.Allgather(np.array(self.outall).flatten(), outallt)
        #    self.outall = outallt.reshape((self.Nb, -1))
        #    if self.Nb > 1 and self.rank == 0:
        #        self.output()
        for inb in range(self.rank*self.Nb//self.size,
                         (self.rank+1)*self.Nb//self.size):
            self.DoQf(inb)
        ##MPI Gathering the result from each rank
        outallt = np.zeros((self.size*np.array(self.outall).size))
        self.comm.Allgather(np.array(self.outall).flatten(), outallt)
        self.outall = outallt.reshape((self.Nb, -1))
        if self.Nb > 1 and self.rank == 0:
            self.output()

    def run_eachkde(self):
        self.outall = []
        #for inb in range(self.Nb):
        ##MPI Distributing the work to each rank
        ### in the case of Nb=1, we still want self.outall
        #if self.Nb == 1:
        #    if self.rank == 0:
        #        self.DoQf(0)
        #else:
        #    for inb in range(self.rank*self.Nb//self.size,
        #                     (self.rank+1)*self.Nb//self.size):
        #        self.DoQf(inb)
        #    ##MPI Gathering the result from each rank
        #    outallt = np.zeros((self.size*np.array(self.outall).size))
        #    self.comm.Allgather(np.array(self.outall).flatten(), outallt)
        #    self.outall = outallt.reshape((self.Nb, -1))
        #    if self.Nb > 1 and self.rank == 0:
        #        self.output()
        #for inb in range(self.rank*self.Nb//self.size,
        #                 (self.rank+1)*self.Nb//self.size):
        for inb in range(self.Nb):
            self.kys = [self.CalcBandW(orgfile, inb=inb) for orgfile in self.orgfiles
                        ]
            self.DoQf(inb)
        if self.Nb > 1 and self.rank == 0:
            self.output()
        ##MPI Gathering the result from each rank
        #outallt = np.zeros((self.size*np.array(self.outall).size))
        #self.comm.Allgather(np.array(self.outall).flatten(), outallt)
        #self.outall = outallt.reshape((self.Nb, -1))


def testrun():
    np.set_printoptions(suppress=True)
    Nb = 1
    elim = [-0.03, 0.07]
    ishist = True
    num = 6400
    rsmodifier = "org"
    variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
    variables = [0.655, 0.0129, 0.200, 0.00208]
    prefix = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/"
    prj = qens_balloon_resamples(runNos=[6202, 6204], elim=elim, Nb=Nb,
                                 ishist=ishist, num=num, rsmodifier=rsmodifier,
                                 prefix=prefix, variables=variables)
    #print(qens_balloon_resamples.__mro__)
    prj.run()
    prj.ishist = False
    prj.run()


#testrun()
