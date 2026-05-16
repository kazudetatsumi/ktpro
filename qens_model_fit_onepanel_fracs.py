#!/usr/bin/env python
import numpy as np
import scipy.optimize as so
from matplotlib import pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as qbr


class qens_model_fit(qbr):
    def __init__(self, runNo, fracs, qsize, resful, io):
        self.runNo = runNo
        self.fracs = np.array(fracs)
        self.qsize = qsize
        self.q2 = (np.arange(self.qsize)*0.1 + 0.2 + 0.05)**2.
        self.resful = resful
        self.io = io

    def getdata(self, frac):
        preprefix = "/home/kazu/desktop/210108/Tatsumi/"
        if self.resful:
            self.orgprefix = preprefix + "from_pca03/wcorr/run" + self.runNo +\
                "/" + frac + "/3540_full/"
        else:
            self.orgprefix = preprefix + "from_pca03/wcorr/run" + self.runNo +\
                "/" + frac + "/"
        self.outfile = self.orgprefix + "outhist" + str(self.runNo) + "m.pkl"
        self.loadfile()
        orghout = self.outall
        if self.io:
            self.outfile = self.orgprefix + "outkdeio" + str(self.runNo) +\
                           "m.pkl"
        else:
            self.outfile = self.orgprefix + "outkde" + str(self.runNo) +\
                           "m.pkl"
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
        fig = plt.figure(figsize=(6.0, 4))
        for fidx, frac in enumerate(self.fracs):
            self.eachrun(fidx, frac, fig)

        plt.subplots_adjust(wspace=0.5, hspace=0.0)
        plt.show()

    def eachrun(self, cidx, frac, fig):
        orghout, orgkout = self.getdata(frac)
        self.plotter(fig, 2,  0, cidx, frac, orghout, 'hist')
        if self.io:
            self.plotter(fig, 2,  1, cidx, frac, orgkout, 'kde')
        else:
            self.plotter(fig, 2,  1, cidx, frac, orgkout, 'kdeb')

    def plotter(self, fig, nr, sidx, cidx, frac, orgout, title):
        mask = np.sum(orgout[:, 4:], axis=1) < -0.5
        pnr = 1+len(self.fracs)*sidx+cidx
        ax = fig.add_subplot(nr, len(self.fracs), pnr)
        ax.scatter(self.q2[~mask], np.log10(orgout[~mask, 0]), marker='o',
                   s=10., ec='k', lw=1., fc='k')
        ax.scatter(self.q2[~mask], np.log10(orgout[~mask, 2]), marker='x',
                   s=10., ec='k', lw=1., fc='k')
        ax.scatter(self.q2[mask], np.log10(orgout[mask, 0]), marker='o',
                   s=10., ec='r', lw=1., fc='r')
        ax.scatter(self.q2[mask], np.log10(orgout[mask, 2]), marker='x',
                   s=10., ec='r', lw=1., fc='r')
        ax.set_ylim(-4., 0.)
        ax.set_yticks([-4, -3, -2, -1, 0])
        ax2 = ax.twinx()
        ax2.scatter(self.q2[~mask], orgout[~mask, 1]*1000, marker='s', s=10.,
                    ec='k', lw=1., fc='k')
        ax2.scatter(self.q2[mask], orgout[mask, 1]*1000, marker='s', s=10.,
                    ec='r', lw=1., fc='r')
        ax2.set_ylim(0., 35)
        ax2.set_yticks([0, 10, 20, 30])

        ax2.text(2.1, 1., title+"_"+frac)
        if pnr >= (nr-1)*len(self.fracs)+1:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=True)
        else:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=False)
        if pnr == nr*len(self.fracs):
            ax.set_xlabel(r'$Q^2 \ (\AA^{-2})$ ')
        if pnr == 2:
            ax.set_ylabel(r'$\log_{10}A_{QENS}\ (Arb. Units)$')
            ax2.set_ylabel(r'$\Gamma\ (\mu eV)$')



def testrun():
    frac = "0125"
    resful = False
    io = False
    if len(sys.argv) >= 2:
        frac = sys.argv[1]
    if len(sys.argv) >= 3:
        if sys.argv[2] == "resful":
            resful = True
        elif sys.argv[2] == "io":
            io = True
    runNo = "4174"
    fracs = ["100", frac]
    #qsize = 11
    #runNos = [6202, 6205, 6203, 6206,  6207]
    #runNos = [6205, 6203, 6206,  6207]
    #temps = [303., 288., 275., 263., 253.]
    #temps = [288., 275., 263., 253.]
    qsize = 17
    prj = qens_model_fit(runNo, fracs, qsize, resful, io)
    prj.run()


testrun()
