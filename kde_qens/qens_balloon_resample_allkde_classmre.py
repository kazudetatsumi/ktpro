#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
# mre stands for multiple qs, ebinned data and channel specific errors
# for least square fitting.
# CI_of_intensities is overriding the original in the 
# qens_balloon_resample_class, to generate rebinned resampled histogram data
# sets.
from mpi4py import MPI
import numpy as np
import scipy.optimize as so
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_classmre import qens_balloon_resamples as sqkr


class qens_balloon_resamples(sqkr):
    def __init__(self, qidx, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1,
                 ishist=False, num=6400, rsmodifier="b", orgmodifier="org",
                 prefixes=["./", "./"], variables=[0.655, 0.0129, 0.200, 0.00208],
                 quiet=False):
        self.qidx = qidx
        self.runNos = runNos
        self.Nb = Nb
        self.gammas = np.zeros((Nb, 2))
        self.elim = elim
        self.ishist = ishist
        self.num = num
        self.rsmodifier = rsmodifier
        self.orgmodifier = orgmodifier
        self.prefixes = prefixes
        self.variables = variables
        self.quiet = quiet
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.leastsq = False
        self.DefineFiles()

    def getrsspectra(self, rsfile, inb=0):
        #super(sqkr, self).__init__(pklfile=rsfile)
        self.pklfile = rsfile
        self.load_pkl()
        print("getrsspectra: chkm slicing spectrab at qidx 0, 1, 2")
        return self.spectrab[inb, 0],\
            self.spectrab[inb, 1],\
            self.spectrab[inb, 2]

    def DefineFiles(self):
        self.rsfiles = []
        self.orgfiles = []
        for runno, prefix in zip(self.runNos, self.prefixes):
            self.rsfiles.append(prefix + "run" + str(runno) + "spectra" +
                                self.rsmodifier + ".pkl." + str(self.qidx))
            self.orgfiles.append(prefix + "run" + str(runno) + "spectra" +
                                 self.orgmodifier + ".pkl." + str(self.qidx))
        if self.rank == 0:
            for prefix in self.prefixes:
                if not self.quiet:
                    print(prefix)
                print("")
                print("rsfiles:", self.rsfiles)
                print("orgfiles:", self.orgfiles)

    def CI_of_intensities(self):
        # Overriding the original method in qens_balloon_resample_class.
        # to rebin the intensities and generate rebinned allhist files.
        if not self.ishist:
            self.kys = [self.CalcBandW(self.orgfiles[0], inb=0)]
        ytrlcs = []
        for inb in range(self.rank*self.Nb//self.size,
                         (self.rank+1)*self.Nb//self.size):
            xt, yt, dummy = self.eachrunno(0, inb)
            xtr, ytr = self.rebin(xt, yt)
            self.icorr()
            xtrl, ytrl = self.limit2(xtr, ytr, self.elim)
            dummy, ytrlc = self.correction(xtrl, ytrl, ytrl)
            ytrlcs.append(ytrlc)
        outallt = np.array(self.comm.gather(np.array(ytrlcs).flatten(),
                           root=0))
        if self.rank == 0:
            self.outall = outallt.reshape((self.Nb, -1))

