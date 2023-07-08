#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
import os
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_classm2 import Sqens_balloon_resamples as qbr


class Sqbr(qbr):
    def __init__(self, qidx, outfile, Nb, runNos=[6202, 6204],
                 elim=[-0.03, 0.07], binwidth1=0.0005, binwidth2=0.0005,
                 binwidth=0.00025):
        ishist = False
        super().__init__(qidx, runNos=runNos, elim=elim, Nb=Nb,
                         ishist=ishist)
        self.CI_of_intensities_io()
        if self.rank == 0:
            self.outfile = outfile
            self.savefile()


def testrun():
    prefix = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/"
    outfile = prefix + "100/q0.9-1.0/test/testm/runno6202/outhist.pkl"
    runNos = [6202, 6204]
    Nb = 20
    elim = [-0.03, 0.07]
    # binwidth1 = 0.0005
    # binwidth2 = 0.0005
    binwidth = 0.00025
    qidx = 0
    Sqbr(qidx, outfile, Nb, runNos=runNos, elim=elim, binwidth=binwidth)


# testrun()
