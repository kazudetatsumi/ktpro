#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
# Rebinning and weights at each bin are added in order to reproduce the results of
# JPSJ  (2021) T. Yamada et al.. 
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_classmre import qens_balloon_resamples as qbr
np.set_printoptions(suppress=True)


class qens_org_classm(qbr):
    def __init__(self, qidx, qsize, elim=[-0.03, 0.07], rsmodifier="orge", ishist=True,
                 orgmodifier="orge", variables=[0.655, 0.0129, 0.200, 0.00208],
                 prefix="./", num=6400, quiet=True, runNos=[6206, 6204]):
        self.qsize = qsize
        Nb = 1
        if ishist:
            self.outfile = "outhist" + str(runNos[0]) + "m.pkl"
        else:
            self.outfile = "outkde" + str(runNos[0]) + "m.pkl"
        qbr.__init__(self, qidx, runNos=runNos, elim=elim, Nb=Nb, ishist=ishist,
                     num=num, rsmodifier=rsmodifier,
                     orgmodifier=orgmodifier, prefix=prefix,
                     variables=variables, quiet=quiet)

    def run_for_mqs(self):
        for qidx in range(0, self.qsize):
            print("CHECK qidx:", qidx)
            self.qidx = qidx
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
