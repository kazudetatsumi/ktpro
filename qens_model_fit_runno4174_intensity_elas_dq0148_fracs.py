#!/usr/bin/env python
import numpy as np
import scipy.optimize as so
from matplotlib import pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as qbr


class qens_model_fit(qbr):
    def __init__(self, runNo, fracs, qsize):
        self.runNo = runNo
        self.fracs = fracs
        self.qsize = qsize
        self.q2 = ((np.arange(self.qsize)+1)*0.148)**2.
        self.u2 = np.zeros((3, len(fracs)))
        self.stdu2 = np.zeros((3, len(fracs)))

    def getdata(self, frac):
        preprefix = "/home/kazu/desktop/210108/Tatsumi/"
        # preprefix = "/Users/kazu/Desktop/210108/Tatsumi/"
        self.orgprefix = preprefix + "from_pca03/wcorr/run" + self.runNo + "/"\
            + frac + "/dq0148/"
        self.stdprefix = preprefix + "from_pca03/wcorr/run" + self.runNo + "/"\
            + frac + "/dq0148/resamples/"
        self.outfile = self.orgprefix + "outhist" + self.runNo + "m.pkl"
        self.loadfile()
        orghout = self.outall
        self.outfile = self.orgprefix + "outkde" + self.runNo + "m.pkl"
        self.loadfile()
        orgkout = self.outall
        self.outfile = self.orgprefix + "outkdeio" + self.runNo + "m.pkl"
        self.loadfile()
        orgkioout = self.outall
        maskh, gammah = self.readorgout(orghout)
        maskkb, gammakb = self.readorgout(orgkout)
        maskk, gammak = self.readorgout(orgkioout)
        errorh, aveh = self.readerror('hist', self.runNo)
        errorkb, avekb = self.readerror('kde', self.runNo)
        errork, avek = self.readerror('kdeio', self.runNo)
        #self.kdefile = self.kdeprefix + "kde3.log"
        #maskk, gammak = self.readlog()
        stdhes = self.readstdhessfromlog(self.runNo)
        return maskh, gammah, maskkb, gammakb, maskk, gammak, errorh, errorkb,\
            errork, stdhes, aveh, avekb

    def stats(self):
        outall = self.outall[np.sum(self.outall[:, 4:], axis=1) > -0.5]
        orderidx1 = np.argsort(outall[:, 2])
        lbs = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1]
        ubs = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.84))], 1]
        return (ubs - lbs) / 2., np.std(outall[:, 2]), np.average(outall[:, 2])

    def readerror(self, hork, runNo):
        error = np.zeros((2, self.qsize))
        ave = np.zeros((self.qsize))
        for qidx in range(0, self.qsize):
            self.outfile = self.stdprefix + \
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
        [u2, b] = coeffs
        y = -u2/3.0*x[mask] + b
        #return ((np.log(t[mask]) - y)/(error[mask]/t[mask]))
        return np.log(t[mask]) - y

    def res_arrhenius(self, coeffs, x, t, error):
        [a, b] = coeffs
        y = a*x + b
        return ((t - y)/error)
        #return t - y

    def optimize(self, variables, gamma, error, mask):
        return so.least_squares(self.res, variables,
                                args=(self.q2, gamma, error, mask))

    def optimize_arrhenius(self, variables, x, t, error):
        out = so.least_squares(self.res_arrhenius, variables,
                               args=(x, t, error))
        s_sq = (self.res_arrhenius(out.x, x, t, error)**2).sum() /\
               (len(t)-len(out.x))
        cov = np.absolute(np.linalg.inv(np.dot(out.jac.T, out.jac))*s_sq)
        return out.x, cov

    def plotter(self, fig, nr, x, y, t, e, mask, title, sidx, fidx):
        #xerr = 2.*(np.arange(self.qsize)*0.1 + 0.05)*0.05
        #mask *= e > 0.0001
        #mask *= t > 0.00001
        pnr = 1+len(self.fracs)*sidx+fidx
        ax = fig.add_subplot(nr, len(self.fracs), pnr)
        ax.plot(x, y, ls='dashed', lw=0.5, c='k')
        ax.errorbar(x[mask], np.log(t[mask]), yerr=np.abs(e[mask]/t[mask]), marker="x",
                    ms=2, elinewidth=1, lw=0, capsize=3, c='k')
        ax.errorbar(x[~mask], np.log(t[~mask]), yerr=np.abs(e[~mask]/t[~mask]), marker="x",
                    ms=2, elinewidth=1, lw=0, capsize=3, c='gray')
        ax.text(1.4, -0.5, title+"_"+str(self.fracs[fidx]), fontsize=7)
        ax.text(1.4, -1.0, '$u^2$={:.2f} +/- {:.2f} '.format(self.u2[sidx, fidx], self.stdu2[sidx, fidx]) + '$\AA^2$', fontsize=7)
        ax.text(1.4, -1.5, 'rel. error +/- {:.1f}% '.format(self.stdu2[sidx, fidx]/np.abs(self.u2[sidx, fidx])*100.), fontsize=7)
        ax.set_ylim(-5., -0.)
        #ax.set_ylim(0, 400.)
        #ax.set_xlim(0., 1.6)
        #ax.set_yticks([0.000, 0.005, 0.010, 0.015, 0.020])
        #ax.set_yticks([0., 5, 10, 15, 20])
        #ax.set_xticks([0., 0.4, 0.8, 1.2, 1.6])
        if pnr >= (nr-1)*len(self.fracs)+1:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=True)
        else:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=False)
        if sidx == 2:
            ax.set_xlabel(r'$Q^2 \ (\AA^{-2})$ ')
        if sidx == 1:
            ax.set_ylabel(r'$lnA_{ELAS}$')

    def plotters(self, x,  ys, es, titles):
        for sidx, title in enumerate(titles):
            plt.errorbar(x, ys[sidx], yerr=es[sidx], marker="x", ms=1,
                         elinewidth=1, lw=0, capsize=3, label=title)
        plt.legend()
        plt.xlim(0.0032, 0.004)
        plt.xlabel('1/T (K)')
        plt.ylabel('lnD')

    def run(self):
        fig = plt.figure(figsize=(2.7*len(self.fracs), 6))
        for fidx, frac in enumerate(self.fracs):
            mask = np.zeros((self.qsize), dtype=bool)
            mask[-2:] = True
            mask[3:5] = True
            self.eachrun(fidx, frac, mask, fig)

        plt.subplots_adjust(wspace=0.3, hspace=0.0)
        plt.legend()
        plt.show()
        #fig2 = plt.figure(figsize=(10, 10))
        #self.plotters(1./self.temps, np.log(self.D), self.stdD/self.D,
        #              ['hist', 'kdeb', 'kde'])
        #for sidx, color in enumerate(['blue', 'orange', 'green']):
        #    out_arrhenius, cov = self.optimize_arrhenius([-1., 1.],
        #                                                 1./self.temps,
        #                                                 np.log(self.D[sidx]),
        #                                                 self.stdD[sidx] /
        #                                                 self.D[sidx])
        #    print(out_arrhenius[0], "+-", (cov**0.5)[0, 0])
        #    plt.plot(1./self.temps, out_arrhenius[0]/self.temps +
        #             out_arrhenius[1], c=color)
        #plt.show()

    def eachsolution(self, fig, sidx, fidx, frac, gamma, error, mask, label):
        #mask *= error > 0.00005
        mask *= np.log(gamma) < -1.0
        out = self.optimize([2.35, 0.52], gamma, error, mask)
        u2 = out.x[0]
        self.u2[sidx, fidx] = u2  # [Angs^2]
        b = out.x[1]
        #y = D*self.q2/(1. + D*tau*self.q2)
        y = -u2/3.0*self.q2 + b
        s_sq = (self.res(out.x, self.q2, gamma, error, mask)**2).sum() /\
               (len(gamma)-len(out.x))
        cov = np.absolute(np.linalg.inv(np.dot(out.jac.T, out.jac))*s_sq)
        self.stdu2[sidx, fidx] = (cov**0.5)[0, 0]

        self.plotter(fig, 3, self.q2, y, gamma,
                     error, mask, label, sidx, fidx)

    def eachrun(self, fidx, frac, mask, fig):
        maskh, gammah, maskkb, gammakb, maskk, gammak, errorh, errorkb,\
         errork, stdhes, aveh, avekb = self.getdata(frac)
        print(frac, len(maskh), len(gammah), len(maskkb), len(gammakb), len(maskk), len(gammak), len(errork[0]), len(stdhes))
        self.eachsolution(fig, 0, fidx, frac, gammah, errorh[0],
        #self.eachsolution(fig, 0, fidx, runNo, gammah, stdhes,
                          mask*~maskh, 'hist')
        self.eachsolution(fig, 1, fidx, frac, gammakb, errorkb[0],
                          mask*~maskkb, 'kdeb')
        self.eachsolution(fig, 2, fidx, frac, gammak, errork[0],
                          mask*~maskk, 'kde')


def testrun():
    runNo = "4174"
    qsize = 12
    fracs = ["100", "050", "025", "0125"]
    if len(sys.argv) >= 2:
        fracs = ["100", sys.argv[1]]

    prj = qens_model_fit(runNo, fracs, qsize)
    prj.run()


testrun()
