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
    Nb = 2800
    Nb = 100
    elim = [-0.03, 0.07]
    #outfile = "outkde_2comps.pkl"
    #outfile = "outkde.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outhist_2comps.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outkde_2comps.pkl"
    outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outkde.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outhist.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/outkdetmp_2comps.pkl"
    #outfile = "./outhist.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outkde.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outhist.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outhisttmp.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outkdetmp.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outkde_2comps.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/outhist_2comps.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/latterhalf_boot/div10/outhist.pkl"
    #outfile = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/latterhalf_boot/div10/outhkde.pkl"
    binwidth1 = 0.0005
    binwidth2 = 0.0005
    binwidth = 0.000375
    #binwidth = 0.002
    print(outfile)
    if "2comps" in outfile:
        variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
    else:
        variables = [0.655, 0.0129, 0.200, 0.00208]
    if 'outhist' in outfile:
        ishist = True
    else:
        ishist = False
    #print(outfile)
    if os.path.isfile(outfile):
        prj = qbr(runNos=[6202, 6204], elim=elim, Nb=Nb)
        prj.outfile = outfile
        prj.loadfile()
        prj.output()
        if prj.rank == 0:
            if len(variables) == 4:
                prj.plot_distribution_single(binwidth)
                print("energy step:", binwidth)
            else:
                prj.plot_distribution(binwidth1, binwidth2)
                print("energy step:", binwidth1, binwidth2)
    else:
        prj = qbr(runNos=[6202, 6204], elim=elim, Nb=Nb, ishist=ishist, variables=variables)
        prj.run()
        if prj.rank == 0:
            if len(variables) == 4:
                prj.plot_distribution_single(binwidth)
            else:
                prj.plot_distribution(binwidth1, binwidth2)
            prj.outfile = outfile
            prj.savefile()


testrun()
