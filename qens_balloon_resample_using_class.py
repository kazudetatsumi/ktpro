#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
import os
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import qens_balloon_resamples as qbr


def testrun():
    Nb = 701
    elim = [-0.03, 0.07]
    #outfile = "outkde_2comps.pkl"
    #outfile = "outkde.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outhist_2comps.pkl"
    outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outkde_2comps.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outkde.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outhist.pkl"
    #outfile = "./outhist.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outkde.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outhist.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outkde_2comps.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outhist_2comps.pkl"
    binwidth1 = 0.001
    binwidth2 = 0.001
    binwidth = 0.0007
    print(outfile)
    if os.path.isfile(outfile):
        prj = qbr(runNos=[6202, 6204], elim=elim, Nb=Nb)
        prj.outfile = outfile
        prj.loadfile()
        prj.output()
        #prj.plot_distribution_single(binwidth)
        prj.plot_distribution(binwidth1, binwidth2)
    else:
        prj = qbr(runNos=[6202, 6204], elim=elim, Nb=Nb)
        prj.run()
        #prj.plot_distribution_single(binwidth)
        prj.plot_distribution(binwidth1, binwidth2)
        prj.outfile = outfile
        prj.savefile()


testrun()
