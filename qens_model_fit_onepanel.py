#!/usr/bin/env python
import numpy as np
import scipy.optimize as so
from matplotlib import pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as qbr


class qens_model_fit(qbr):
    def __init__(self, runNos, temps, qsize):
        self.runNos = runNos
        self.temps = np.array(temps)
        self.qsize = qsize
        self.q2 = (np.arange(self.qsize)*0.1 + 0.2 + 0.05)**2.
        self.D = np.zeros((3, len(runNos)))
        self.stdD = np.zeros((3, len(runNos)))

    def getdata(self, runNo):
        preprefix = "/home/kazu/desktop/210108/Tatsumi/"
        # preprefix = "/Users/kazu/Desktop/210108/Tatsumi/"
        self.orgprefix = preprefix + "from_pca03/wcorr/test/100/orgs/run4174/bin0001/"
        self.outfile = self.orgprefix + "outhist" + str(runNo) + "m.pkl"
        self.loadfile()
        orghout = self.outall
        self.outfile = self.orgprefix + "outkde" + str(runNo) + "m.pkl"
        self.loadfile()
        orgkout = self.outall
        return orghout, orgkout

    def stats(self):
        outall = self.outall[np.sum(self.outall[:, 4:], axis=1) > -0.5]
        orderidx1 = np.argsort(outall[:, 1])
        lbs = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1]
        ubs = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.84))], 1]
        return (ubs - lbs) / 2., np.std(outall[:, 1]), np.average(outall[:, 1])

    def readerror(self, hork, runNo):
        error = np.zeros((2, self.qsize))
        ave = np.zeros((self.qsize))
        for qidx in range(0, self.qsize):
            self.outfile = self.stdprefix + "runno" + str(runNo) +\
                           "/out" + hork + ".pkl." + str(qidx)
            self.loadfile()
            error[0, qidx], error[1, qidx], ave[qidx] = self.stats()
        return error, ave

    def readorgout(self, orgout):
        return np.sum(orgout[:, 4:], axis=1) < -0.5,  orgout[:, 2]

    def readlog(self):
        mask = []
        gamma = []
        for line in open(self.kdefile):
            mask.append("-1" in line[:-1].split(']')[0])
            gamma.append(float(line[:-1].split('[')[2].split()[1]))
        return np.array(mask), np.array(gamma)

    def readstdhessfromlog(self, runNo):
        std = []
        linecount = 999999
        for line in open(self.orgprefix + "std-hist-" + str(runNo)):
            if "cov**0.5" in line:
                linecount = 1
            elif linecount > 0:
                linecount -= 1
            elif linecount == 0:
                std.append(float(line[:-1].split()[1]))
                linecount = 999999
        return np.array(std)

    def res(self, coeffs, x, t, error, mask):
        [D, tau] = coeffs
        y = D*x[mask]/(1. + D*tau*x[mask])
        return ((t[mask] - y)/error[mask])
        #return t[mask] - y

    def res_arrhenius(self, coeffs, x, t, error):
        [a, b] = coeffs
        y = a*x + b
        return ((t - y)/error)
        #return t - y

    def optimize(self, variables, gamma, error, mask):
        return so.least_squares(self.res, variables, bounds=(0, np.inf),
                                args=(self.q2, gamma, error, mask))

    def optimize_arrhenius(self, variables, x, t, error):
        out = so.least_squares(self.res_arrhenius, variables,
                               args=(x, t, error))
        s_sq = (self.res_arrhenius(out.x, x, t, error)**2).sum() /\
               (len(t)-len(out.x))
        cov = np.absolute(np.linalg.inv(np.dot(out.jac.T, out.jac))*s_sq)
        return out.x, cov

    def run(self):
        #fig = plt.figure(figsize=(16, 4))
        fig = plt.figure(figsize=(3, 4))
        for cidx, runNo in enumerate(self.runNos):
            self.eachrun(cidx, runNo, fig)

        plt.subplots_adjust(wspace=0.6, hspace=0.0)
        plt.show()

    def eachrun(self, cidx, runNo, fig):
        orghout, orgkout = self.getdata(runNo)
        self.plotter(fig, 2,  0, cidx, runNo, orghout, 'hist')
        self.plotter(fig, 2,  1, cidx, runNo, orgkout, 'kdeb')

    def plotter(self, fig, nr, sidx, cidx, runNo, orgout, title):
        #xerr = 2.*(np.arange(self.qsize)*0.1 + 0.05)*0.05
        #mask *= e > 0.0001
        #mask *= t > 0.00001
        mask = np.sum(orgout[:, 4:], axis=1) < -0.5
        pnr = 1+len(self.runNos)*sidx+cidx
        ax = fig.add_subplot(nr, len(self.runNos), pnr)
        ax.scatter(self.q2[~mask], np.log10(orgout[~mask, 0]), marker='o', s=10., ec='k', lw=1., fc='k')
        ax.scatter(self.q2[~mask], np.log10(orgout[~mask, 2]), marker='x', s=10., ec='k', lw=1., fc='k')
        ax.scatter(self.q2[mask], np.log10(orgout[mask, 0]), marker='o', s=10., ec='r', lw=1., fc='r')
        ax.scatter(self.q2[mask], np.log10(orgout[mask, 2]), marker='x', s=10., ec='r', lw=1., fc='r')
        #ax.scatter(self.q2[mask], np.log10(orgout[mask, 0]), marker='o', s=10., ec='gray', lw=1., fc='gray')
        #ax.scatter(self.q2[mask], np.log10(orgout[mask, 2]), marker='x', s=10., ec='gray', lw=1., fc='gray')
        #ax.scatter(self.q2[~mask], np.log10(orgout[~mask, 3]), marker='s', s=8., ec='k', lw=0.5, fc='w')
        ax.set_ylim(-4., 0.)
        ax2 = ax.twinx()
        ax2.scatter(self.q2[~mask], orgout[~mask, 1]*1000, marker='.', s=10., ec='k', lw=1., fc='k')
        ax2.scatter(self.q2[mask], orgout[mask, 1]*1000, marker='.', s=10., ec='r', lw=1., fc='r')
        ax2.set_ylim(0., 25)

        ##ax.plot(x, y*1000.)
        ##ax.errorbar(x[mask], t[mask]*1000., yerr=e[mask]*1000., marker="x",
        ##            ms=2, elinewidth=1, lw=0, capsize=3)
        ##ax.errorbar(x[~mask], t[~mask]*1000., yerr=e[~mask]*1000., marker="x",
        ##            ms=2, elinewidth=1, lw=0, capsize=3, c='gray')
        #ax2.text(0.1, 0.017*1000., title+"_"+str(self.runNos[cidx]))
        ##ax.text(0.1, 0.015*1000., '{:.1f} +/- {:.1f}'.format(self.D[sidx, cidx]*1000., self.stdD[sidx, cidx]*1000.))
        ##ax.set_ylim(-1, 0.022*1000.)
        #ax.set_xlim(0., 1.6)
        #ax.set_yticks([0.000, 0.005, 0.010, 0.015, 0.020])
        ##ax.set_yticks([0., 5, 10, 15, 20])
        ##ax.set_xticks([0., 0.4, 0.8, 1.2, 1.6])
        if pnr >= (nr-1)*len(self.runNos)+1:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=True)
        else:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=False)
        ##if pnr == 13:
        if pnr == nr*len(self.runNos):
        #if pnr == nr*len(self.runNos)-2:
            ax.set_xlabel(r'$Q^2 \ (\AA^{-2})$ ')
        #if pnr == 6:
        if pnr == 2:
            ax.set_ylabel(r'$\log_{10}A_{QENS}\ (Arb. Units)$')
            ax2.set_ylabel(r'$\Gamma\ (\mu eV)$')

    #def eachsolution(self, fig, sidx, cidx, runNo, orgout, label):
        #out = self.optimize([0.05, 52], gamma, error, mask)
        #D = out.x[0]
        #self.D[sidx, cidx] = D
        #tau = out.x[1]
        #y = D*self.q2/(1. + D*tau*self.q2)
        #s_sq = (self.res(out.x, self.q2, gamma, error, mask)**2).sum() /\
        #       (len(gamma)-len(out.x))
        #cov = np.absolute(np.linalg.inv(np.dot(out.jac.T, out.jac))*s_sq)
        #self.stdD[sidx, cidx] = (cov**0.5)[0, 0]

        #self.plotter(fig, 3, self.q2, y, gamma,
        #             error, mask, label, sidx, cidx)


def testrun():
    runNos = [6202, 6205, 6203, 6206, 6207]
    temps = [303., 288., 275., 263., 253.]
    runNos = [4174]
    temps = [298.]
    #qsize = 11
    #runNos = [6202, 6205, 6203, 6206,  6207]
    #runNos = [6205, 6203, 6206,  6207]
    #temps = [303., 288., 275., 263., 253.]
    #temps = [288., 275., 263., 253.]
    qsize = 17
    prj = qens_model_fit(runNos, temps, qsize)
    prj.run()


testrun()