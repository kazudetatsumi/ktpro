#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_classm import qens_balloon_resamples as qbr
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
            self.outfile = "outhist" + str(runNos[0]) + "m.pkl"
        else:
            if io:
                if isnovariablebw:
                    self.outfile = "outkdeionovbw" + str(runNos[0]) + "m.pkl"
                else:
                    self.outfile = "outkdeio" + str(runNos[0]) + "m.pkl"
            else:
                if isnovariablebw:
                    if prefixes[0] == prefixes[1]:
                        self.outfile = "outkdenovbw" + str(runNos[0]) + "m.pkl"
                    else:
                        self.outfile = "outkdenovbwprefixes" + str(runNos[0]) + "m.pkl"
                else:
                    self.outfile = "outkde" + str(runNos[0]) + "m.pkl"
        qbr.__init__(self, qidx, runNos=runNos, elim=elim, Nb=Nb,
                     ishist=ishist, num=num, M=M, WinFunc=WinFunc,
                     winparam=winparam, rsmodifier=rsmodifier,
                     orgmodifier=orgmodifier, prefixes=prefixes,
                     variables=variables, quiet=quiet,
                     ispltchk=ispltchk, isnovariablebw=isnovariablebw)

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

    def combine_qs(self):
        _out = self.get_out()
        self.outfile += ".bak"
        __out = self.get_out()
        self.outall = np.concatenate((__out, _out), axis=0)
        self.outfile = self.outfile.split("bak")[0] + "combined"
        self.savefile()

    def get_out(self):
        self.loadfile()
        return np.array(self.outall)


def samplerun():
    runNos = [6206, 6204]
    ishist = True
    qidx = 0
    prj = qens_org_classm(qidx, runNos=runNos, ishist=ishist)
    prj.run_for_mqs()


#samplerun()
