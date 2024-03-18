#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
# Rebinning and weights at each bin are added in order to reproduce the
# results of JPSJ  (2021) T. Yamada et al..
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_classme import qens_balloon_resamples as qbr
np.set_printoptions(suppress=True)


class qens_org_classm(qbr):
    def __init__(self, qidx, qsize, elim=[-0.03, 0.07], rsmodifier="org",
                 ishist=True, orgmodifier="org", variables=[0.655, 0.0129,
                 0.200, 0.00208], prefixes=["./", "./"], num=6400, M=160,
                 winparam=1, WinFunc='Boxcar', quiet=True, io=False,
                 runNos=[6206, 6204], ispltchk=False, isnovariablebw=False):
        self.qsize = qsize
        Nb = 1
        if ishist:
            self.outfile = "outhist" + str(runNos[0]) + "me.pkl"
        else:
            if io:
                if isnovariablebw:
                    self.outfile = "outkdeionovbw" + str(runNos[0]) + "me.pkl"
                else:
                    self.outfile = "outkdeio" + str(runNos[0]) + "me.pkl"
            else:
                if isnovariablebw:
                    self.outfile = "outkdenovbw" + str(runNos[0]) + "me.pkl"
                else:
                    self.outfile = "outkde" + str(runNos[0]) + "me.pkl"
        qbr.__init__(self, qidx, runNos=runNos, elim=elim, Nb=Nb,
                     ishist=ishist, num=num, M=M, winparam=winparam,
                     WinFunc=WinFunc, rsmodifier=rsmodifier,
                     orgmodifier=orgmodifier, prefixes=prefixes,
                     variables=variables, quiet=quiet, ispltchk=ispltchk,
                     isnovariablebw=isnovariablebw)

    def run_for_mqs(self):
        for qidx in range(self.qsize):
            print("CHECK qidx:", qidx)
            self.qidx = qidx
            if "kdeio" in self.outfile:
                self.run_eachkde_io()
            else:
                self.run_eachkde()
            if qidx == 0:
                outall = self.outall
            else:
                outall = np.concatenate((outall, self.outall), axis=0)
        if self.rank == 0:
            self.outall = outall
            self.savefile()


def samplerun():
    runNos = [6206, 6204]
    ishist = True
    qidx = 0
    prj = qens_org_classm(qidx, runNos=runNos, ishist=ishist)
    prj.run_for_mqs()


#samplerun()
