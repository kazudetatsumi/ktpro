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
    def __init__(self, qidx, elim=[-0.03, 0.07], rsmodifier="org", ishist=True,
                 orgmodifier="org", variables=[0.655, 0.0129, 0.200, 0.00208],
                 prefix="./", num=6400, quiet=True, runNos=[6206, 6204]):
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
        for qidx in range(0, 6):
            self.qidx = qidx
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
