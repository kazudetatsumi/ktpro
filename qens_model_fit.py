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

    def getdata(self, runNo):
        # preprefix = "/home/kazu/desktop/210108/Tatsumi/"
        preprefix = "/Users/kazu/Desktop/210108/Tatsumi/"
        self.orgprefix = preprefix + "from_pca03/wcorr/test/100/orgs/"
        self.stdprefix = preprefix + "from_pca03/wcorr/test/100/qs/"
        self.kdeprefix = preprefix + "winparam_exam_" + str(runNo) + \
            "/160_1_0000025io_Boxcar/n8000/"
        self.outfile = self.orgprefix + "outhist" + str(runNo) + "m.pkl"
        self.loadfile()
        orghout = self.outall
        self.outfile = self.orgprefix + "outkde" + str(runNo) + "m.pkl"
        self.loadfile()
        orgkout = self.outall
        maskh, gammah = self.readorgout(orghout)
        maskkb, gammakb = self.readorgout(orgkout)
        errorh = self.readerror('hist', runNo)
        errork = self.readerror('kde', runNo)
        self.kdefile = self.kdeprefix + "kde2.log"
        maskk, gammak = self.readlog()
        stdhes = self.readstdhessfromlog(runNo)
        return maskh, gammah, maskkb, gammakb, maskk, gammak, errorh, errork,\
            stdhes

    def stats(self):
        orderidx1 = np.argsort(self.outall[:, 1])
        lbs = self.outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1]
        ubs = self.outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.84))], 1]
        return (ubs - lbs) / 2., np.std(self.outall[:, 1])

    def readerror(self, hork, runNo):
        error = np.zeros((2, self.qsize))
        for qidx in range(0, self.qsize):
            self.outfile = self.stdprefix + "runno" + str(runNo) +\
                           "/out" + hork + ".pkl." + str(qidx)
            self.loadfile()
            error[0, qidx], error[1, qidx] = self.stats()
        return error

    def readorgout(self, orgout):
        return np.sum(orgout[:, 4:], axis=1) < -0.5,  orgout[:, 1]

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

    def optimize(self, variables, gamma, error, mask):
        return so.least_squares(self.res, variables, bounds=(0, np.inf),
                                args=(self.q2, gamma, error, mask))

    def plotter(self, fig, nr, pnr,  x, y, t, e, mask, title, runNo):
        #heavymask = ~mask
        #xerr = 2.*(np.arange(self.qsize)*0.1 + 0.05)*0.05
        ax = fig.add_subplot(nr, len(self.runNos), pnr)
        ax.errorbar(x[mask], t[mask], yerr=e[mask], marker="x", ms=2,
                    elinewidth=1, lw=0, capsize=3)
        #ax.errorbar(x[heavymask], t[heavymask], yerr=e[heavymask], marker="x",
        #            ms=2, elinewidth=1, lw=0, capsize=3, c='gray')
        ax.text(0.1, 0.017, title+"_"+str(runNo))
        ax.set_ylim(0., 0.022)
        ax.set_xlim(0., 1.6)
        ax.set_yticks([0.000, 0.005, 0.010, 0.015, 0.020])
        ax.plot(x, y)
        if pnr >= (nr-1)*len(self.runNos)+1:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=True)
        else:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=False)

    def run(self):
        fig = plt.figure(figsize=(10, 10))
        for cidx, runNo in enumerate(self.runNos):
            if runNo == 6202:
                mask = np.array([False, True, True, True, True, True, True,
                                 True, False, False, False])
            else:
                mask = np.array([False, True, True, True, True, True, True,
                                 True, True, True, True])
            self.eachrun(cidx, runNo, mask, fig)

        plt.subplots_adjust(wspace=0.2, hspace=0.0)
        plt.show()
        plt.scatter(1/self.temps, np.log(self.D[0]), label='hist')
        plt.scatter(1/self.temps, np.log(self.D[1]), label='kb')
        plt.scatter(1/self.temps, np.log(self.D[2]), label='k')
        plt.legend()
        plt.show()

    def eachsolution(self, fig, sidx, cidx, runNo, gamma, error, mask, label):
        out = self.optimize([0.05, 52], gamma, error, mask)
        D = out.x[0]
        self.D[sidx, cidx] = D
        tau = out.x[1]
        y = D*self.q2/(1. + D*tau*self.q2)
        s_sq = (self.res(out.x, self.q2, gamma, error, mask)**2).sum() /\
               (len(gamma)-len(out.x))
        cov = np.absolute(np.linalg.inv(np.dot(out.jac.T, out.jac))*s_sq)
        stdD = (cov**0.5)[0, 0]
        print(stdD, label, runNo)

        self.plotter(fig, 3, 1+len(self.runNos)*sidx+cidx, self.q2, y, gamma,
                     error, mask, label, runNo)

    def eachrun(self, cidx, runNo, mask, fig):
        maskh, gammah, maskkb, gammakb, maskk, gammak, errorh, errork, stdhes\
            = self.getdata(runNo)

        if runNo < 6206:
            maskhh = np.array([False, True, True, True, True, True, True,
                               True, True, True, True])
        else:
            maskhh = np.array([False, False, True, True, True, True, True,
                               True, True, True, True])

        self.eachsolution(fig, 0, cidx, runNo, gammah, errorh[0],
        # self.eachsolution(fig, 0, cidx, runNo, gammah, stdhes,
                          mask*~maskh*maskhh, 'hist')
        self.eachsolution(fig, 1, cidx, runNo, gammakb, errork[0],
                          mask*~maskkb, 'kdeb')
        self.eachsolution(fig, 2, cidx, runNo, gammak, errork[0],
                          mask*~maskk, 'kde')


def testrun():
    runNos = [6202, 6205, 6203, 6206, 6207]
    temps = [303., 288., 275., 263., 253.]
    qsize = 11
    prj = qens_model_fit(runNos, temps, qsize)
    prj.run()


testrun()