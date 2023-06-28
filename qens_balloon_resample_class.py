#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
from mpi4py import MPI
from qens_kde_resampled import qens_kde_resampled as qkr


class Sqens_balloon_resamples(qkr):
    def __init__(self, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1,
                 ishist=False, num=6400, rsmodifier="b", orgmodifier="org",
                 prefix="./", variables=[0.655, 0.0129, 0.200, 0.00208],
                 quiet=False):
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
        self.DefineFiles()

    def getrsspectra(self, rsfile, inb=0):
        qkr.__init__(self, pklfile=rsfile)
        return self.spectrab[inb, 0, :], self.spectrab[inb, 1, :]

    def CalcBandW(self, orgfile, inb=0):
        if self.ishist:
            if self.rank == 0 and not self.quiet:
                print("skipping KDE because ishist", self.ishist)
            self.y = "dummy"
        else:
            qkr.__init__(self, pklfile=orgfile)
            self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :],
                     num=self.num)
        return self.y

    def balloon(self, ky, sy):
        bw = np.interp(sy[0], ky[1], ky[2])
        idx = sy[1].nonzero()
        sy_nz = sy[1][idx]
        t_nz = sy[0][idx]
        yv = np.zeros_like(sy[0])
        dt = min(np.diff(sy[0]))
        for k in range(yv.shape[0]):
            yv[k] = np.sum(sy_nz * dt * self.Gauss(sy[0][k]-t_nz, bw[k]))
        ##yv = yv * np.max(ky[0]) / np.max(yv)
        return yv

    def Gauss(self, x, w):
        return 1.0/((2.0*np.pi)**0.5*w)*np.exp(-(x/w)**2/2.0)

    def DefineFiles(self):
        self.rsfiles = []
        self.orgfiles = []
        for runno in self.runNos:
            self.rsfiles.append(self.prefix + "run" + str(runno) + "spectra" +
                                self.rsmodifier + ".pkl")
            self.orgfiles.append(self.prefix + "run" + str(runno) + "spectra" +
                                 self.orgmodifier + ".pkl")
        if self.rank == 0:
            if not self.quiet:
                print(self.prefix)
            print("")
            print("rsfiles:", [rsf.split(self.prefix)[1] for rsf in
                               self.rsfiles])
            print("orgfiles:", [orgf.split(self.prefix)[1] for orgf in
                                self.orgfiles])

    def eachrunno(self, fidx, inb):
        sy = self.getrsspectra(self.rsfiles[fidx], inb)
        if self.ishist:
            return sy[0], sy[1], sy[1]
        else:
            syb = self.balloon(self.kys[fidx], sy)
            return sy[0], syb, sy[1]

    def DoQf(self, inb):
        import matplotlib.pyplot as plt
        xt, yt, yth = self.eachrunno(0, inb)
        xd, yd, ydh = self.eachrunno(1, inb)
        xt, yt = self.rebin(xt, yt)
        xd, yd = self.rebin(xd, yd)
        self.icorr()
        print("CHK elim:", self.elim)
        xtl, ytl = self.limit2(xt, yt, self.elim)
        xdl, ydl = self.limit2(xd, yd, self.elim)
        if inb == 0 and self.rank == 0:
            if np.sum(xtl - xdl) > 0.000001:
                print('WARNING, check x_tf - x_df')
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        self.bg = 0.
        self.check_out(inb, self.optimize(xdl, ydl, ytl,
                                          variables=self.variables))

        [alpha, gamma, delta, base] = self.outall[-1][0:4]
        yqens = alpha*self.convloreorg(ydl, gamma, xdl)
        y = yqens + delta*ydl + base
        plt.plot(xdl*1000, y, c='k')
        plt.scatter(xdl*1000, ytl, marker='o', s=18, fc='None', ec='k')
        plt.plot(xdl*1000, yqens, ls='dotted', c='k')
        plt.ylabel('Intensity (Arb. Units)')
        plt.xlabel(r'$Energy\ (\mu eV)$')
        plt.show()

    def getbins(self):
        bins = []
        for line in open("/home/kazu/Ebin20150709.txt"):
            bins.append(float(line[:-1].split()[0]))
        self.bins = np.array(bins)

    def rebin(self, x, y):
        nbins = self.bins.shape[0]
        xr = np.zeros((nbins-1))
        yr = np.zeros((nbins-1))
        for ibdx in range(nbins-1):
            xr[ibdx] = (self.bins[ibdx] + self.bins[ibdx+1])/2.
        for _x, _y in zip(x, y):
            for ibdx in range(nbins-1):
                if _x + 0.0000125 >= self.bins[ibdx] and _x < self.bins[ibdx+1]:
                    yr[ibdx] += _y
        for ibdx in range(nbins-1):
            yr[ibdx] /= ((self.bins[ibdx+1] - self.bins[ibdx])/0.000025)
        return xr, yr

    def run(self):
        self.kys = [self.CalcBandW(orgfile, inb=0) for orgfile in self.orgfiles
                    ]
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

    def run_eachkde(self):
        self.outall = []
        for inb in range(self.Nb):
            self.kys = [self.CalcBandW(orgfile, inb=inb) for orgfile in
                        self.orgfiles]
            self.DoQf(inb)
        self.outall = np.array(self.outall)
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
    #prefix = preprefix + "100/"
    #prefix = preprefix + "0125/back/test5/6208/nonmpi_test/"
    #prefix = preprefix + "0125/back/test5/6208/whole/"
    prefix = preprefix + "test/0125/back/test5/6208/nonmpi_test/dnapca03/"
    #prj = qens_balloon_resamples(runNos=[6202, 6204], elim=elim, Nb=Nb,
    prj = qens_balloon_resamples(runNos=[6207, 6204], elim=elim, Nb=Nb,
                                 ishist=ishist, num=num, rsmodifier=rsmodifier,
                                 prefix=prefix, variables=variables)
    #print(qens_balloon_resamples.__mro__)
    prj.run()
    prj.ishist = False
    prj.run()


#testrun()
