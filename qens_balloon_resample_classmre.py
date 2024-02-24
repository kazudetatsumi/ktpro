#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
# mre stands for multiple qs, ebinned data and channel specific errors
# for least square fitting.
# CI_of_intensities is overriding the original in the 
# qens_balloon_resample_class, to generate rebinned resampled histogram data
# sets.
from mpi4py import MPI
import numpy as np
import scipy.optimize as so
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as sqkr


class qens_balloon_resamples(sqkr):
    def __init__(self, qidx, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1, 
                 ishist=False, num=6400, rsmodifier="b", orgmodifier="org",
                 prefixes=["./", "./"], variables=[0.655, 0.0129, 0.200, 0.00208],
                 quiet=False):
        self.qidx = qidx
        self.runNos = runNos
        self.Nb = Nb
        self.gammas = np.zeros((Nb, 2))
        self.elim = elim
        self.ishist = ishist
        self.num = num
        self.rsmodifier = rsmodifier
        self.orgmodifier = orgmodifier
        self.prefixes = prefixes
        self.variables = variables
        self.quiet = quiet
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.leastsq = False
        self.DefineFiles()

    def eachrunno(self, fidx, inb):
        #print("CHK", self.rsfiles)
        sy = self.getrsspectra(self.rsfiles[fidx], inb)
        if self.ishist:
            return sy[0], sy[1], sy[2]
        else:
            syb = self.balloon(self.kys[fidx], sy)
            return sy[0], syb, sy[1]

    def getrsspectra(self, rsfile, inb=0):
        super(sqkr, self).__init__(pklfile=rsfile)
        print("getrsspectra: chkm slicing spectrab at qidx 0, 1, 3")
        return self.spectrab[inb, 0, self.qidx],\
            self.spectrab[inb, 1, self.qidx],\
            self.spectrab[inb, 3, self.qidx]
            #self.spectrab[inb, 2, self.qidx]

    def CalcBandW(self, orgfile, inb=0):
        if self.ishist:
            if self.rank == 0 and not self.quiet:
                print("skipping KDE because ishist", self.ishist)
            self.y = "dummy"
        else:
            super(sqkr, self).__init__(pklfile=orgfile)
            print("CalcBandW: chkm slicing spectrab at qidx")
            self.spectrab = self.spectrab[:, :, self.qidx, :]
            self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :],
                     num=self.num)
        return self.y

    def DoQf(self, inb):
        import matplotlib.pyplot as plt
        xt, yt, et = self.eachrunno(0, inb)
        xd, yd, ed = self.eachrunno(1, inb)
        xtr, ytr = self.rebin(xt, yt)
        xdr, ydr = self.rebin(xd, yd)
        etr = self.rebine(xt, et)
        self.icorr()
        print("CHK elim:", self.elim)
        xtl, ytl = self.limit2(xtr, ytr, self.elim)
        xdl, ydl = self.limit2(xdr, ydr, self.elim)
        xtl, etl = self.limit2(xtr, etr, self.elim)
        if inb == 0 and self.rank == 0:
            if np.sum(xtl - xdl) > 0.000001:
                print('WARNING, check x_tf - x_df')
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        etlc, dummy = self.correction(xtl, etl, ytl)
        self.bg = 0.
        self.check_out(inb, self.optimize(xdl, ydlc, ytlc, etlc,
        #self.check_out(inb, self.optimize(xdl, ydl, ytl, etl,
                                          variables=self.variables))

        [alpha, gamma, delta, base] = self.outall[-1][0:4]
        yqens = alpha*self.convloreorg(ydlc, gamma, xdl)
        y = yqens + delta*ydl + base
        plt.plot(xdl*1000, y, c='k')
        #plt.errorbar(xdl*1000, ytlc, yerr=etl, marker='o', ms=2., lw=0, elinewidth=1 )
        plt.errorbar(xdl*1000, ytl, yerr=etl, marker='o', ms=3., lw=0,
                     elinewidth=1, mfc='None')
        plt.plot(xdl*1000, yqens, ls='dotted', c='k')
        plt.ylabel('Intensity (Arb. Units)')
        plt.xlabel(r'$Energy\ (\mu eV)$')
        plt.show()

    def getbins(self):
        #bins = []
        #for line in open("/home/kazu/Ebin20150709.txt"):
        #    bins.append(float(line[:-1].split()[0]))
        #self.bins = np.array(bins)
        self.bins = np.arange(-0.03, 0.121, 0.001)

    def rebin(self, x, y):
        nbins = self.bins.shape[0]
        xr = np.zeros((nbins-1))
        yr = np.zeros((nbins-1))
        de = x[1] - x[0]
        for ibdx in range(nbins-1):
            xr[ibdx] = (self.bins[ibdx] + self.bins[ibdx+1])/2.
        for _x, _y in zip(x, y):
            for ibdx in range(nbins-1):
                if _x >= self.bins[ibdx] and _x < self.bins[ibdx+1]:
                    yr[ibdx] += _y
        for ibdx in range(nbins-1):
            yr[ibdx] /= ((self.bins[ibdx+1] - self.bins[ibdx])/de)
        return xr, yr

    def rebine(self, x, e):
        nbins = self.bins.shape[0]
        xr = np.zeros((nbins-1))
        er = np.zeros((nbins-1))
        de = x[1] - x[0]
        for ibdx in range(nbins-1):
            xr[ibdx] = (self.bins[ibdx] + self.bins[ibdx+1])/2.
        for _x, _e in zip(x, e):
            if _e < 100000.:
                for ibdx in range(nbins-1):
                    if _x >= self.bins[ibdx] and _x < self.bins[ibdx+1]:
                        er[ibdx] += _e**2
        er = er**0.5
        for ibdx in range(nbins-1):
            er[ibdx] /= ((self.bins[ibdx+1] - self.bins[ibdx])/de)
        return er

    def res(self, coeffs, x, d, t, e):
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
            y = alpha*self.convloreorg(d, gamma, x)\
                + delta*d + base
        if len(coeffs) == 3:
            [alpha, gamma, delta] = coeffs
            y = alpha*self.convlore(d, gamma, x)\
                + delta*d + self.bg
        xl, dif = self.limit2(x[e > 0.], (t-y)[e > 0.]/(e[e > 0.]), self.elim)
        return dif

    def optimize(self, x, yd, yt, et, variables=[6.e-6, 2.e-2, 1.e-6, 4.e-3, 7.e-3, 3.e-1]):
        # leastsq
        if self.leastsq:
            out = so.leastsq(self.res, variables, args=(x, yd, yt, et),
                             full_output=1, epsfcn=0.0001)
            return out
        # least_squares
        else:
            bounds = (0, np.inf)
            out = so.least_squares(self.res, variables, bounds=bounds,
                                   args=(x, yd, yt, et))
            #if self.rank == 0:
            print(out.active_mask, out.success, out.x)
            #out = so.least_squares(self.res, variables, args=(x, yd, yt))
            _out = [out.x, np.linalg.inv(np.dot(out.jac.T, out.jac))]
            s_sq = (self.res(_out[0], x, yd, yt, et)**2).sum() / (len(yt)-len(_out[0]))
            cov = np.absolute(_out[1]*s_sq)
            print("cov**0.5")
            print(cov**0.5)
            return [out.x, out.success, out.active_mask, cov]


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
