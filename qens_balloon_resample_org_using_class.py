#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import qens_balloon_resamples as qbr


def showkdeorhist(ishist):
    if ishist:
        print('Entering Results on HIST')
    else:
        print('Entering Results on KDE')


def testrun():
    np.set_printoptions(suppress=True)
    Nb = 1
    elim = [-0.03, 0.07]
    ishist = True
    num = 6400
    rsmodifier = "org"
    orgmodifier = "org"
    variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
    variables = [0.655, 0.0129, 0.200, 0.00208]
    prefix = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/100/"
    print("prefix:", prefix)
    print("number of bins of tin in kde:", num)
    quiet = True
    prj = qbr(runNos=[6202, 6204], elim=elim, Nb=Nb, ishist=ishist, num=num,
              rsmodifier=rsmodifier, orgmodifier=orgmodifier, prefix=prefix,
              variables=variables, quiet=quiet)
    if prj.rank == 0:
        print("prefix:", prefix)
        print("number of bins of tin in kde:", num)
        showkdeorhist(prj.ishist)
    prj.run_eachkde()
    prj.ishist = False
    if prj.rank == 0:
        showkdeorhist(prj.ishist)
    prj.run_eachkde()
    prj.variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
    prj.ishist = True
    if prj.rank == 0:
        showkdeorhist(prj.ishist)
        prj.outall = []
        prj.DoQf(0)
        prj.ishist = False
        showkdeorhist(prj.ishist)
        prj.DoQf(0)


def testruns():
    np.set_printoptions(suppress=True)
    Nb = 1
    elim = [-0.03, 0.07]
    ishist = True
    num = 6400
    variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
    variables = [0.655, 0.0129, 0.200, 0.00208]
    prefix = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/test/0125/back/test5/orgs/"
    #print("")
    #print("prefix:", prefix)
    #print("number of bins of tin in kde:", num)
    quiet = True
    for i in range(1,4):
        rsmodifier = "org" + str(i)
        orgmodifier = "org" + str(i)
        prj = qbr(runNos=[6202, 6204], elim=elim, Nb=Nb, ishist=ishist, num=num,
                  rsmodifier=rsmodifier, orgmodifier=orgmodifier, prefix=prefix,
                  variables=variables, quiet=quiet)
        if prj.rank == 0:
            if i == 1:
                print("prefix:", prefix)
                print("number of bins of tin in kde:", num)
            showkdeorhist(prj.ishist)
        prj.run()
        prj.ishist = False
        if prj.rank == 0:
            showkdeorhist(prj.ishist)
        prj.run()
        prj.variables = [0.8, 0.01, 0.24, 0.0002, 0.001, 1.2]
        prj.ishist = True
        if prj.rank == 0:
            showkdeorhist(prj.ishist)
            prj.outall = []
            prj.DoQf(0)
            prj.ishist = False
            showkdeorhist(prj.ishist)
            prj.DoQf(0)


testrun()
#testruns()
