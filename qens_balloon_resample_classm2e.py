#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
# Modified init and DefineFiles for multi qs.
# Kazuyoshi TATSUMI 2023/06/12
import pickle
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from mpi4py import MPI
from qens_balloon_resample_class import Sqens_balloon_resamples as qkr


class Sqens_balloon_resamples(qkr):
    def __init__(self, qidx, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1,
                 ishist=False, num=6400, rsmodifier="b", orgmodifier="org",
                 prefix="./", variables=[0.655, 0.0129, 0.200, 0.00208],
                 quiet=False, ispltchk=False):
        self.qidx = qidx
        self.runNos = runNos
        self.Nb = Nb
        self.gammas = np.zeros((Nb, 2))
        self.elim = elim
        self.ishist = ishist
        self.num = num
        self.rsmodifier = rsmodifier
        self.orgmodifier = orgmodifier
        self.prefix = prefix
        self.variables = variables
        self.quiet = quiet
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.leastsq = False
        self.ispltchk = ispltchk
        self.DefineFiles()

    def DefineFiles(self):
        self.rsfiles = []
        self.orgfiles = []
        for runno in self.runNos:
            self.rsfiles.append(self.prefix + "run" + str(runno) + "spectra" +
                                self.rsmodifier + ".pkl." + str(self.qidx))
            self.orgfiles.append(self.prefix + "run" + str(runno) + "spectra" +
                                 self.orgmodifier + ".pkl." + str(self.qidx))
        if self.rank == 0:
            if not self.quiet:
                print(self.prefix)
            print("")
            print("rsfiles:", [rsf.split(self.prefix)[1] for rsf in
                               self.rsfiles])
            print("orgfiles:", [orgf.split(self.prefix)[1] for orgf in
                                self.orgfiles])

    def geterrorbars(self):
        with open("outallkde.pkl."+str(self.qidx), 'rb') as f:
            dat = pickle.load(f)['out']
        return np.std(dat, axis=0)

    def geterrorbarsio(self):
        with open("outallkdeio.pkl."+str(self.qidx), 'rb') as f:
            dat = pickle.load(f)['out']
        return np.std(dat, axis=0)

    def res(self, coeffs, x, d, t):
        if len(coeffs) == 6:
            [alpha1, gamma1, alpha2, gamma2,  delta, base] = coeffs
            y = alpha1*self.convlore(d, gamma1, x)\
                + alpha2*self.convlore(d, gamma2, x)\
                + delta*d + base
        if len(coeffs) == 5:
            [alpha1, gamma1, alpha2, gamma2,  delta] = coeffs
            y = alpha1*self.convlore(d, gamma1, x)\
                + alpha2*self.convlore(d, gamma2, x)\
                + delta*d + self.bg
        if len(coeffs) == 4:
            [alpha, gamma, delta, base] = coeffs
            y = alpha*self.convlore(d, gamma, x)\
                + delta*d + base
        if len(coeffs) == 3:
            [alpha, gamma, delta] = coeffs
            y = alpha*self.convlore(d, gamma, x)\
                + delta*d + self.bg
        xl, dif = self.limit2(x, (t-y)/self.etl, self.elim)
        #xl, dif = self.limit2(x, t-y, self.elim)
        return dif

    def run(self):
        self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles
                    ]
        self.etl = self.geterrorbars()
        self.outall = []
        for inb in range(self.rank*self.Nb//self.size,
                         (self.rank+1)*self.Nb//self.size):
            self.DoQf(inb)
        ##MPI Gathering the result from each rank
        outallt = np.zeros((self.size*np.array(self.outall).size))
        self.comm.Allgather(np.array(self.outall).flatten(), outallt)
        self.outall = outallt.reshape((self.Nb, -1))
        if self.Nb > 1 and self.rank == 0:
            self.output()

    def run_io(self):
        self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles
                    ]
        self.check_idata()
        self.etl = self.geterrorbarsio()
        self.outall = []
        for inb in range(self.rank*self.Nb//self.size,
                         (self.rank+1)*self.Nb//self.size):
            self.kyos = [self.eachrunno_io(fidx, inb) for fidx in range(2)]
            self.kyios = [self.io(kyo, kyi) for kyo, kyi in
                          zip(self.kyos, self.kyis)]
            self.DoQfio(inb)
        ##MPI Gathering the result from each rank
        outallt = np.zeros((self.size*np.array(self.outall).size))
        self.comm.Allgather(np.array(self.outall).flatten(), outallt)
        self.outall = outallt.reshape((self.Nb, -1))
        if self.Nb > 1 and self.rank == 0:
            self.output()


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
                                  rsmodifier=rsmodifier, prefix=prefix,
                                  variables=variables)
    # print(qens_balloon_resamples.__mro__)
    prj.run()
    prj.ishist = False
    prj.run()


# testrun()
