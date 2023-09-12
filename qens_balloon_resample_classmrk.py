#!/usr/bin/env python
# This script balloon-estimates the re-sampled qens 1D intensity profiles with
# the kernel band widths optimized on the original data, and fits the estimated
# target density with the estimated device density.
# Kazuyoshi TATSUMI 2023/02/23
# The histogram bin widths are locally varied according to VKDE.
# Kazuyoshi TATSUMI 2023/09/09
from mpi4py import MPI
import numpy as np
import scipy.optimize as so
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as sqkr


class qens_balloon_resamples(sqkr):
    def __init__(self, qidx, runNos=[6202, 6204], elim=[-0.03, 0.07], Nb=1,
                 ishist=False, num=6400, rsmodifier="b", orgmodifier="org",
                 prefix="./", variables=[0.655, 0.0129, 0.200, 0.00208],
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
        self.prefix = prefix
        self.variables = variables
        self.quiet = quiet
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.leastsq = False
        self.DefineFiles()

    def getrsspectra(self, rsfile, inb=0):
        super(sqkr, self).__init__(pklfile=rsfile)
        print("getrsspectra: chkm slicing spectrab at qidx")
        return self.spectrab[inb, 0, self.qidx],\
            self.spectrab[inb, 1, self.qidx]

    def CalcBandW(self, orgfile, inb=0):
        #if self.ishist:
        #    if self.rank == 0 and not self.quiet:
        #        print("skipping KDE because ishist", self.ishist)
        #    self.y = "dummy"
        #else:
        #    super(sqkr, self).__init__(pklfile=orgfile)
        #    print("CalcBandW: chkm slicing spectrab at qidx")
        #    self.spectrab = self.spectrab[:, :, self.qidx, :]
        #    self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :],
        #             num=self.num)
        super(sqkr, self).__init__(pklfile=orgfile)
        print("CalcBandW: chkm slicing spectrab at qidx")
        self.spectrab = self.spectrab[:, :, self.qidx, :]
        self.kde(self.spectrab[inb, 0, :], self.spectrab[inb, 2, :],
                 num=self.num)
        return self.y

    def DoQf(self, inb):
        #import matplotlib.pyplot as plt
        xt, yt, yth = self.eachrunno(0, inb)
        xd, yd, ydh = self.eachrunno(1, inb)
        xt, yt = self.rebin(xt, yt, self.kys[0])
        xd, yd = self.rebin(xd, yd, self.kys[1])
        #xt, yt = self.rebin(xt, yt)
        #xd, yd = self.rebin(xd, yd)
        self.icorr()
        print("CHK elim:", self.elim)
        xtl, ytl = self.limit2(xt, yt, self.elim)
        xdl, ydl = self.limit2(xd, yd, self.elim)
        #if inb == 0 and self.rank == 0:
        #    if np.sum(xtl - xdl) > 0.000001:
        #        print('WARNING, check x_tf - x_df')
        ydlc, ytlc = self.correction(xdl, ydl, xtl, ytl)
        self.bg = 0.
        #self.check_out(inb, self.optimize(xdl, ydlc, ytlc,
        self.check_out(inb, self.optimize(xdl, ydlc, xtl, ytlc,
                                          variables=self.variables))

        #[alpha, gamma, delta, base] = self.outall[-1][0:4]
        #yqens = alpha*self.convloreorg(ydlc, gamma, xdl, xtl)
        #y = yqens + delta*self.convloreorg(ydlc, 0.0000001, xdl, xtl) + base
        #plt.plot(xtl*1000, y, c='k')
        #plt.scatter(xtl*1000, ytlc, marker='o', s=18, fc='None', ec='k')
        #plt.plot(xtl*1000, yqens, ls='dotted', c='k')
        #plt.ylabel('Intensity (Arb. Units)')
        #plt.xlabel(r'$Energy\ (\mu eV)$')
        #plt.show()

    def correction(self, x, yd, xt, yt):
        x = x + 2.085
        xt = xt + 2.085
        yd = yd / (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        yt = yt / (self.k[0] + self.k[1]*xt + self.k[2]*xt**2 + self.k[3]*xt**3)
        return yd, yt

    def getbins(self, x, kys):
        #import matplotlib.pyplot as plt
        print("Entering getbins")
        h = kys[2]
        tie = kys[1]
        #plt.plot(tie, h)
        #plt.show()
        horg = tie[1] - tie[0]
        print('CHKk!!!', h[np.argmin(np.abs(tie - 0.1))])
        h0 = h[np.argmin(np.abs(tie))]
        ux = [horg/2.]
        while ux[-1] < x[-1]:
            idx = np.min(np.argsort(np.abs(tie - ux[-1]))[0:2])
            if (ux[-1] - tie[idx])*(ux[-1] - tie[idx+1]) >= 0.:
                break
            else:
                _h = (h[idx + 1] - h[idx])/(tie[idx+1] - tie[idx])*(
                        ux[-1] - tie[idx]) + h[idx]
                __h = _h/h0*horg
                __h = _h
                ux.append(ux[-1] + __h)
        dx = [-horg/2.]
        while dx[-1] > x[0]:
            idx = np.min(np.argsort(np.abs(tie - dx[-1]))[0:2])
            if (dx[-1] - tie[idx])*(dx[-1] - tie[idx+1]) >= 0.:
                break
            else:
                _h = (h[idx + 1] - h[idx])/(tie[idx+1] - tie[idx])*(
                        dx[-1] - tie[idx]) + h[idx]
                __h = _h/h0*horg
                __h = _h
                dx.append(dx[-1] - __h)
        dx.reverse()
        self.bins = np.array(dx + ux)
        #import matplotlib.pyplot as plt
        #plt.scatter(self.bins, np.zeros(self.bins.shape[0]), marker='|')
        #plt.show()

#    def getbins(self):
#        bins = []
#        #for line in open("/home/kazu/Ebin20150709.txt"):
#        #    bins.append(float(line[:-1].split()[0]))
#        #self.bins = np.array(bins)
#        #self.bins = np.arange(-0.03, 0.12025, 0.00025)
#        self.bins = np.arange(-0.03, 0.121, 0.001)
#        #self.bins = np.arange(-0.03, 0.121, 0.0005)
#        #self.bins = np.arange(-0.03, 0.122, 0.002)
#        #self.bins = np.arange(-0.03, 0.123, 0.003)
#        #self.bins = np.arange(-0.03, 0.122, 0.004)
#        #self.bins = np.arange(-0.03, 0.125, 0.005)
#        #self.bins = np.arange(-0.03, 0.13, 0.01)
#        #self.bins = np.arange(-0.03, 0.12125, 0.00275)
#
#    def rebin(self, x, y):
#        nbins = self.bins.shape[0]
#        xr = np.zeros((nbins-1))
#        yr = np.zeros((nbins-1))
#        for ibdx in range(nbins-1):
#            xr[ibdx] = (self.bins[ibdx] + self.bins[ibdx+1])/2.
#        for _x, _y in zip(x, y):
#            for ibdx in range(nbins-1):
#                if _x + 0.0000125 >= self.bins[ibdx] and _x + 0.0000125 < self.bins[ibdx+1]:
#                    yr[ibdx] += _y
#        for ibdx in range(nbins-1):
#            yr[ibdx] /= ((self.bins[ibdx+1] - self.bins[ibdx])/0.000025)
#        return xr, yr

    def rebin(self, x, y, kys):
        self.getbins(x, kys)
        nbins = self.bins.shape[0]
        xr = np.zeros((nbins-1))
        yr = np.zeros((nbins-1))
        for ibdx in range(nbins-1):
            xr[ibdx] = (self.bins[ibdx] + self.bins[ibdx+1])/2.
        for _x, _y in zip(x, y):
            for ibdx in range(nbins-1):
                if _x + 0.0000125 >= self.bins[ibdx] and _x + 0.0000125 < self.bins[ibdx+1]:
                    yr[ibdx] += _y
        for ibdx in range(nbins-1):
            yr[ibdx] /= ((self.bins[ibdx+1] - self.bins[ibdx])/0.000025)
        return xr, yr

    def optimize(self, x, yd, xt, yt,
                 variables=[6.e-6, 2.e-2, 1.e-6, 4.e-3, 7.e-3, 3.e-1]):
        # leastsq
        if self.leastsq:
            out = so.leastsq(self.res, variables, args=(x, yd, xt, yt),
                             full_output=1, epsfcn=0.0001)
            return out 
        # least_squares
        else:
            bounds = (0, np.inf)
            out = so.least_squares(self.res, variables, bounds=bounds,
                                   args=(x, yd, xt, yt))
            #if self.rank == 0:
            print(out.active_mask, out.success, out.x)
            #out = so.least_squares(self.res, variables, args=(x, yd, yt))
            _out = [out.x, np.linalg.pinv(np.dot(out.jac.T, out.jac))]
            s_sq = (self.res(_out[0], x, yd, xt, yt)**2).sum() / (
                    len(yt)-len(_out[0]))
            cov = np.absolute(_out[1]*s_sq)
            print("cov**0.5")
            print(cov**0.5)
            return [out.x, out.success, out.active_mask, cov]

    def res(self, coeffs, x, d, xt, t):
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
        # We use convloreorg for histograms whose binwidths varied with x.
        if len(coeffs) == 4:
            [alpha, gamma, delta, base] = coeffs
            y = alpha*self.convloreorg(d, gamma, x, xt)\
                + delta*self.convloreorg(d, 0.0000001, x, xt) + base
                #+ delta*d + base
        if len(coeffs) == 3:
            [alpha, gamma, delta] = coeffs
            y = alpha*self.convlore(d, gamma, x)\
                + delta*d + self.bg
        xl, dif = self.limit2(xt, t-y, self.elim)
        return dif

    def convloreorg(self, f, gamma, x, xt):
        _convd = np.zeros_like(xt)
        for _x, _f, in zip(x, f):
            _convd += self.fun_lore(xt - _x, gamma)*_f
        return _convd


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
