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
        self.D = np.zeros((3, len(fracs)))
        self.stdD = np.zeros((3, len(fracs)))
        self.tau = np.zeros((3, len(fracs)))
        self.stdtau = np.zeros((3, len(fracs)))

    def getdata(self, frac):
        prepreprefix = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/"
        # preprefix = "/Users/kazu/Desktop/210108/Tatsumi/"
        rebin = "rebin_hist_Ebin20150709"
        rebin = "rebin_0002"
        #rebin = "rebin_0010"
        #rebin = "rebin_0007"
        nM = "n200"
        self.prefix = prepreprefix + "run" + self.runNo + "/"\
            + frac + "/dq0148/"
        self.outfile = self.prefix + "/resamples/" + rebin +\
            "/outhist" + self.runNo + "mr.pkl"
        self.loadfile()
        orghout = self.outall
        #self.outfile = self.orgprefix + "outkde" + self.runNo + "m.pkl"
        self.outfile = self.prefix + "/" + nM + "/outkde" + self.runNo +\
            "m.pkl"
        self.loadfile()
        orgkout = self.outall
        #self.outfile = self.prefix + "/moni/outkdeio" + self.runNo + "m.pkl"
        #self.loadfile()
        #orgkioout = self.outall
        maskh, gammah = self.readorgout(orghout)
        maskkb, gammakb = self.readorgout(orgkout)
        #maskk, gammak = self.readorgout(orgkioout)
        self.outfile = self.prefix + "/resamples/" + rebin
        maskh2, errorh, aveh = self.readerror('histr', self.runNo)
        self.outfile = self.prefix + "/resamples/" + nM
        maskkb2, errorkb, avekb = self.readerror('kde', self.runNo)
        #self.outfile = self.prefix + "/resamples/"
        #maskk2, errork, avek = self.readerror('kdeio', self.runNo)
        #self.kdefile = self.kdeprefix + "kde3.log"
        #maskk, gammak = self.readlog()
        #stdhes = self.readstdhessfromlog(self.runNo)
        stdhes = self.readstdhessfromorgout(orghout)
        #return maskh, maskh2, gammah, maskkb, maskkb2, gammakb, maskk, maskk2,\
        #    gammak, errorh, errorkb, errork, stdhes, aveh, avekb
        return maskh, maskh2, gammah, maskkb, maskkb2, gammakb, \
            errorh, errorkb, stdhes, aveh, avekb

    def stats(self):
        outall = self.outall[np.sum(self.outall[:, 4:], axis=1) > -0.5]
        #print('CHK', outall.shape)
        orderidx1 = np.argsort(outall[:, 1])
        #print("CHK", orderidx1.shape)
        print(outall.shape[0], int(self.outall.shape[0]/2))
        lbs = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1]
        ubs = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.84))], 1]
        return outall.shape[0] < int(self.outall.shape[0]/2), (ubs - lbs) / 2.,\
            np.std(outall[:, 1]), np.average(outall[:, 1])

    def readerror(self, hork, runNo):
        error = np.zeros((2, self.qsize))
        ave = np.zeros((self.qsize))
        mask = np.zeros((self.qsize), dtype=bool)
        print(hork)
        prestr = str(self.outfile)
        for qidx in range(0, self.qsize):
            self.outfile = prestr + "/out" + hork + ".pkl." + str(qidx)
            #print(self.outfile)
            self.loadfile()
            mask[qidx], error[0, qidx], error[1, qidx], ave[qidx] =\
                self.stats()
        return mask, error, ave

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

    def readstdhessfromorgout(self, orgout):
        return (orgout[:, 8:].reshape(-1, 4, 4)**0.5)[:, 1, 1]

    def res(self, coeffs, x, t, error):
        [D, tau] = coeffs
        y = D*x/(1. + D*tau*x)
        return ((t - y)/error)
        #return t[mask] - y

    def res_arrhenius(self, coeffs, x, t, error):
        [a, b] = coeffs
        y = a*x + b
        return ((t - y)/error)
        #return t - y

    def optimize(self, variables, q2, gamma, error):
        return so.least_squares(self.res, variables, bounds=(0, np.inf),
                                args=(q2, gamma, error))

    def optimize_arrhenius(self, variables, x, t, error):
        out = so.least_squares(self.res_arrhenius, variables,
                               args=(x, t, error))
        s_sq = (self.res_arrhenius(out.x, x, t, error)**2).sum() /\
               (len(t)-len(out.x))
        cov = np.absolute(np.linalg.inv(np.dot(out.jac.T, out.jac))*s_sq)
        return out.x, cov

    def plotter(self, fig, nr, x, y, t, R2, e, mask, mask2, title, sidx,
                fidx):
        #xerr = 2.*(np.arange(self.qsize)*0.1 + 0.05)*0.05
        #mask *= e > 0.0001
        #mask *= t > 0.00001
        print(mask, title, mask2, title)
        pnr = 1+len(self.fracs)*sidx+fidx
        ax = fig.add_subplot(nr, len(self.fracs), pnr)
        ax.plot(x, y*1000., ls='dashed', c='k')
        ax.errorbar(x[~mask*~mask2], t[~mask*~mask2]*1000.,
                    yerr=e[~mask*~mask2]*1000., marker="x", ms=2, elinewidth=1,
                    lw=0, capsize=3, c='k')
        ax.errorbar(x[mask*~mask2], t[mask*~mask2]*1000.,
                    yerr=e[mask*~mask2]*1000., marker="x",
                    ms=2, elinewidth=1, lw=0, capsize=3, c='blue')
        ax.errorbar(x[mask2*~mask], t[mask2*~mask]*1000.,
                    yerr=e[mask2*~mask]*1000., marker="x",
                    ms=2, elinewidth=1, lw=0, capsize=3, c='red')
        ax.errorbar(x[mask2*mask], t[mask2*mask]*1000.,
                    yerr=e[mask2*mask]*1000., marker="x",
                    ms=2, elinewidth=1, lw=0, capsize=3, c='pink')
        ax.text(0.1, 0.020*1000., title+"_"+str(self.fracs[fidx]))
        ax.text(1.6, 0.020*1000., 'R2={:.0f}%'.format(R2*100.))
        ax.text(0.1, 0.018*1000., 'D={:.0f} +/- {:.0f} '.format(self.D[sidx, fidx], self.stdD[sidx, fidx]) + '$10^9\AA^2s^{-1}$')
        ax.text(0.1, 0.016*1000., 'rel. diff. {:.0f}% '.format(self.D[sidx, fidx]/self.D[sidx, 0]*100.-100.))
        ax.text(0.1, 0.014*1000., 'tau={:.0f} +/- {:.0f} '.format(self.tau[sidx, fidx], self.stdtau[sidx, fidx]) + 'ps')
        ax.text(0.1, 0.012*1000., 'rel. diff. {:.0f}% '.format(self.tau[sidx, fidx]/self.tau[sidx, 0]*100.-100.))
        ax.set_ylim(-1, 0.022*1000.)
        #ax.set_xlim(0., 1.6)
        #ax.set_yticks([0.000, 0.005, 0.010, 0.015, 0.020])
        ax.set_yticks([0., 5, 10, 15, 20])
        #ax.set_xticks([0., 0.4, 0.8, 1.2, 1.6])
        if pnr >= (nr-2)*len(self.fracs)+1:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=True)
        else:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=False)
        #if pnr == 5 or pnr == 6:
        if sidx == 1:
            ax.set_xlabel(r'$Q^2 \ (\AA^{-2})$ ')
        if sidx == 1:
            ax.set_ylabel(r'$\Gamma\ (\mu eV)$')

    def plotters(self, x,  ys, es, titles):
        for sidx, title in enumerate(titles):
            plt.errorbar(x, ys[sidx], yerr=es[sidx], marker="x", ms=1,
                         elinewidth=1, lw=0, capsize=3, label=title)
        plt.legend()
        plt.xlim(0.0032, 0.004)
        plt.xlabel('1/T (K)')
        plt.ylabel('lnD')

    def run(self):
        fig = plt.figure(figsize=(3*len(self.fracs), 8))
        for fidx, frac in enumerate(self.fracs):
            #mask = np.ones((self.qsize), dtype=bool)
            #mask[0:2] = False
            self.eachrun(fidx, frac, fig)

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

    def eachsolution(self, fig, sidx, fidx, frac, gamma, error, mask, mask2,
                     label):
        #mask *= error > 0.00005
        #mask *= gamma > 0.001
        #print(mask, mask2, label)
        mask[0:2] = True
        _mask = ~mask*~mask2
        _q2 = self.q2[_mask]
        _gamma = gamma[_mask]
        _error = error[_mask]
        out = self.optimize([0.05, 52], _q2, _gamma, _error)
        D = out.x[0]
        tau = out.x[1]
        # We should mulitiply hbar to convert angular freq. to energy.
        self.D[sidx, fidx] = D/4.1355667*1000.*2.0*3.1415926  # [Angs^2s-1 * 10^9]
        self.tau[sidx, fidx] = tau*4.1355667/2.0/3.1415926    # [ps]
        y = D*self.q2/(1. + D*tau*self.q2)
        s_sq = (self.res(out.x, _q2, _gamma, _error)**2).sum() /\
               (len(_gamma)-len(out.x))
        cov = np.absolute(np.linalg.inv(np.dot(out.jac.T, out.jac))*s_sq)
        self.stdD[sidx, fidx] = (cov**0.5)[0, 0] / 4.1355667 * 1000.*2.0*3.1415926
        self.stdtau[sidx, fidx] = (cov**0.5)[1, 1] * 4.1355667/2.0/3.1415926
        R2 = 1. - np.sum((_gamma - y[_mask])**2.)/np.sum((_gamma - np.mean(_gamma))**2.)

        self.plotter(fig, 3, self.q2, y, gamma, R2,
                     error, mask, mask2, label, sidx, fidx)

    def eachrun(self, fidx, frac, fig):
        #maskh, maskh2, gammah, maskkb, maskkb2, gammakb, maskk, maskk2,\
        #     gammak, errorh, errorkb,  errork, stdhes, aveh, avekb =\
        #     self.getdata(frac)
        maskh, maskh2, gammah, maskkb, maskkb2, gammakb,\
             errorh, errorkb, stdhes, aveh, avekb =\
             self.getdata(frac)
        print(frac, len(maskh), len(gammah), len(maskkb), len(gammakb), len(stdhes))
        #maskhh = np.ones((self.qsize), dtype=bool)
        #if fidx >= 2:
        #    maskhh[1] = False
        self.eachsolution(fig, 0, fidx, frac, gammah, errorh[1],
        #self.eachsolution(fig, 0, fidx, frac, gammah, stdhes,
                          maskh, maskh2, 'hist')
        self.eachsolution(fig, 1, fidx, frac, gammakb, errorkb[1],
                          maskkb, maskkb2, 'kdeb')
        #self.eachsolution(fig, 2, fidx, frac, gammak, errork[1],
        #                  maskk, maskk2, 'kde')


def testrun():
    runNo = "6207"
    qsize = 12
    fracs = ["100", "050", "025", "0125", "0100", "0050"]
    fracs = ["100", "0125", "0050"]
    fracs = ["100"]
    if len(sys.argv) >= 2:
        fracs = ["100", sys.argv[1]]

    prj = qens_model_fit(runNo, fracs, qsize)
    prj.run()


testrun()
