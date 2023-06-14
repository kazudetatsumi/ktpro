#!/usr/bin/env python
import numpy as np
import scipy.optimize as so
from matplotlib import pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as qbr


class qens_model_fit(qbr):
    def __init__(self, runNos, qsize):
        self.runNos = runNos
        self.qsize = qsize
        self.q2 = (np.arange(self.qsize)*0.1 + 0.2 + 0.05)**2.

    def getdata(self, runNo):
        preprefix = "/home/kazu/desktop/210108/Tatsumi/"
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
        return maskh, gammah, maskkb, gammakb, maskk, gammak, errorh, errork

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

    def res(self, coeffs, x, t, error, mask):
        [D, tau] = coeffs
        y = D*x[mask]/(1. + D*tau*x[mask])
        return ((t[mask] - y)/error[mask])
        #return t[mask] - y

    def optimize(self, variables, gamma, error, mask):
        return so.least_squares(self.res, variables, bounds=(0, np.inf),  args=(self.q2, gamma, error, mask))

    def plotter(self, fig, nr, pnr,  x, y, t, e, mask, title, runNo):
        heavymask = ~mask
        ax = fig.add_subplot(nr, len(self.runNos), pnr)
        ax.errorbar(x[mask], t[mask], yerr=e[mask], marker="x", ms=2, elinewidth=1, lw=0, capsize=3)
        ax.errorbar(x[heavymask], t[heavymask], yerr=e[heavymask], marker="x", ms=2, elinewidth=1, lw=0, capsize=3, c='gray')
        ax.text(0.1, 0.017, title+"_"+str(runNo))
        ax.set_ylim(0., 0.022)
        ax.set_xlim(0., 1.6)
        ax.set_yticks([0.000, 0.005, 0.010, 0.015, 0.020])
        ax.plot(x, y)
        if pnr >= (nr-1)*len(self.runNos)+1:
            ax.tick_params(direction='in', top=True, right=True, labelbottom=True)
        else:
            ax.tick_params(direction='in', top=True, right=True, labelbottom=False)

    def run(self):
        fig = plt.figure(figsize=(10, 10))
        for cidx, runNo in enumerate(self.runNos):
            if runNo == 6202:
                mask = np.array([False, True, True, True, True, True, True, True, False, False, False])
            else:
                mask = np.array([False, True, True, True, True, True, True, True, True, True, True])
            self.eachrun(cidx, runNo, mask, fig)

        plt.subplots_adjust(wspace=0.2, hspace=0.0)
        plt.show()

    def eachrun(self, cidx, runNo, mask, fig):
        maskh, gammah, maskkb, gammakb, maskk, gammak, errorh, errork =\
            self.getdata(runNo)
        print(runNo, errorh[0])

        out = self.optimize([0.05, 52], gammah, errorh[0], mask*~maskh)
        D = out.x[0]
        tau = out.x[1]
        y = D*self.q2/(1. + D*tau*self.q2)
        self.plotter(fig, 3, 1+cidx, self.q2, y, gammah, errorh[0], mask*~maskh, 'hist', runNo)

        out = self.optimize([0.05, 40], gammakb, errork[0], mask*~maskkb)
        D = out.x[0]
        tau = out.x[1]
        y = D*self.q2/(1. + D*tau*self.q2)
        self.plotter(fig, 3, 1+len(self.runNos)+cidx, self.q2, y, gammakb, errork[0], mask*~maskkb, 'kde_balloon', runNo)

        out = self.optimize([0.05, 40], gammak, errork[0], mask*~maskk)
        D = out.x[0]
        tau = out.x[1]
        y = D*self.q2/(1. + D*tau*self.q2)
        self.plotter(fig, 3, 1+len(self.runNos)*2+cidx, self.q2, y, gammak, errork[0], mask*~maskk, 'kde', runNo)


def testrun():
    #runNos = [6202, 6205, 6203, 6206, 6207]
    runNos = [6202, 6205]
    qsize = 11
    prj = qens_model_fit(runNos, qsize)
    prj.run()



testrun()
