#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
import os
import sys
sys.path.append("/home/kazu/ktpro")
#from qens_balloon_resample_classm2r import Sqens_balloon_resamples as qbr
from qens_balloon_resample_classmrk import qens_balloon_resamples as qbr


class Sqbr(qbr):
    def __init__(self, qidx, outfile, Nb, num=6400, runNos=[6202, 6204],
                 elim=[-0.03, 0.07], binwidth1=0.0005, binwidth2=0.0005,
                 binwidth=0.00025, io=False):
        if "2comps" in outfile:
            variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
        else:
            variables = [0.655, 0.0129, 0.200, 0.00208]
        if 'outhist' in outfile:
            ishist = True
        else:
            ishist = False
        if os.path.isfile(outfile):
            super().__init__(qidx, num=num, runNos=runNos, elim=elim, Nb=Nb)
            self.outfile = outfile
            self.loadfile()
            self.output()
            if self.rank == 0:
                if len(variables) == 4:
                    self.plot_distribution_single(binwidth)
                    print("energy step:", binwidth)
                else:
                    self.plot_distribution(binwidth1, binwidth2)
                    print("energy step:", binwidth1, binwidth2)
        elif io:
            super().__init__(qidx, runNos=runNos, elim=elim, Nb=Nb,
                             ishist=ishist, variables=variables)
            print("IO")
            self.run_io()
            if self.rank == 0:
                self.outfile = outfile
                self.savefile()
        elif not io:
            super().__init__(qidx, num=num, runNos=runNos, elim=elim, Nb=Nb,
                             ishist=ishist, variables=variables)
            #self.getbins()
            self.run()
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
