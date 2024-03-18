#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
# Modified init and DefineFiles for multi qs.
# Kazuyoshi TATSUMI 2023/06/12

import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from mpi4py import MPI
from qens_balloon_resample_class import Sqens_balloon_resamples as qkr


class Sqens_balloon_resamples(qkr):
    def __init__(self, qidx, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1,
                 ishist=False, num=6400, M=160, winparam=1, WinFunc='Boxcar',
                 rsmodifier="b", orgmodifier="org", prefixes=["./", "./"],
                 variables=[0.655, 0.0129, 0.200, 0.00208], quiet=False,
                 ispltchk=False, isnovariablebw=False):
        self.qidx = qidx
        self.runNos = runNos
        self.Nb = Nb
        self.gammas = np.zeros((Nb, 2))
        self.elim = elim
        self.ishist = ishist
        self.num = num
        self.M = M
        self.winparam = winparam
        self.WinFunc = WinFunc
        self.rsmodifier = rsmodifier
        self.orgmodifier = orgmodifier
        self.prefixes = prefixes
        self.variables = variables
        self.quiet = quiet
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.leastsq = False
        self.ispltchk = ispltchk
        self.isnovariablebw = isnovariablebw
        self.DefineFiles()

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


def testrun():
    np.set_printoptions(suppress=True)
    Nb = 1
    elim = [-0.03, 0.07]
    ishist = True
    num = 6400
    rsmodifier = "org"
    variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
    variables = [0.655, 0.0129, 0.200, 0.00208]
    preprefix = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/"
    # prefix = preprefix + "100/"
    # prefix = preprefix + "0125/back/test5/6208/nonmpi_test/"
    # prefix = preprefix + "0125/back/test5/6208/whole/"
    prefix = preprefix + "test/0125/back/test5/6208/nonmpi_test/dnapca03/"
    # prj = qens_balloon_resamples(runNos=[6202, 6204], elim=elim, Nb=Nb,
    prj = Sqens_balloon_resamples(runNos=[6207, 6204], elim=elim, Nb=Nb,
                                  ishist=ishist, num=num,
                                  rsmodifier=rsmodifier,
                                  prefixes=[prefix, prefix],
                                  variables=variables)
    # print(qens_balloon_resamples.__mro__)
    prj.run()
    prj.ishist = False
    prj.run()


# testrun()
