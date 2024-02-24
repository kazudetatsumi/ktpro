#!/usr/bin/env python
import numpy as np
import scipy.optimize as so
from matplotlib import pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from qens_balloon_resample_class import Sqens_balloon_resamples as qbr


class qens_model_fit(qbr):
    def __init__(self, runNo, fracs, qsize, rebin, nM):
        self.runNo = runNo
        self.fracs = fracs
        self.qsize = qsize
        self.rebin = rebin
        self.nM = nM
        self.q2 = ((np.arange(self.qsize)+1)*0.148)**2.
        self.D = np.zeros((3, len(fracs)))
        self.stdD = np.zeros((3, len(fracs)))
        self.tau = np.zeros((3, len(fracs)))
        self.stdtau = np.zeros((3, len(fracs)))

    def getdata(self, frac):
        prepreprefix = "/home/kazu/desktop/210108/Tatsumi/from_pca03/wcorr/"
        # preprefix = "/Users/kazu/Desktop/210108/Tatsumi/"
        #rebin = "rebin_hist_Ebin20150709"
        #rebin = "rebin_000025"
        #rebin = "norebin"
        #nM = "n800"
        self.prefix = prepreprefix + "run" + self.runNo + "/"\
            + frac + "/dq0148/"
        #self.outfile = self.prefix + "/resamples/" + self.rebin +\
        #    "/outhist" + self.runNo + "mr.pkl"
            #"/outhist" + self.runNo + "mrk.pkl"
        #print(self.outfile)
        self.outfile = self.prefix + "/teste/withcorrection/outhist" + self.runNo + "m.pkl"
        self.loadfile()
        orghout = self.outall
        #self.outfile = self.orgprefix + "outkde" + self.runNo + "m.pkl"
        self.outfile = self.prefix + "/" + self.nM + "/outkde" + self.runNo +\
            "me.pkl"
        #self.outfile = self.prefix + "/" + self.nM + "/outkdenovbw" + self.runNo +\
        #    "me.pkl"
        #print(self.outfile)
        self.loadfile()
        orgkout = self.outall
        self.outfile = self.prefix + "/moni/" + self.nM +"/outkdeionovbw" + self.runNo + "m.pkl"
        #print(self.outfile)
        self.loadfile()
        orgkioout = self.outall
        maskh, gammah = self.readorgout(orghout)
        maskkb, gammakb = self.readorgout(orgkout)
        maskk, gammak = self.readorgout(orgkioout)
        self.outfile = self.prefix + "/resamples/" + self.rebin
        #maskh2, errorh, aveh = self.readerror('histrk', self.runNo)
        maskh2, errorh, aveh = self.readerror('histr', self.runNo)
        self.outfile = self.prefix + "/resamples/" + self.nM
        maskkb2, errorkb, avekb = self.readerror('kdee', self.runNo)
        self.outfile = self.prefix + "/resamples/" + self.nM
        maskk2, errork, avek = self.readerror('kdeionovbw', self.runNo)
        #self.kdefile = self.kdeprefix + "kde3.log"
        #maskk, gammak = self.readlog()
        #stdhes = self.readstdhessfromlog(self.runNo)
        stdhes = self.readstdhessfromorgout(orghout)
        print("CHECK stdhes:", stdhes)
        return maskh, maskh2, gammah, maskkb, maskkb2, gammakb, maskk, maskk2,\
            gammak, errorh, errorkb, errork, stdhes, aveh, avekb

    def stats(self):
        outall = self.outall[np.sum(self.outall[:, 4:], axis=1) > -0.5]
        #print('CHK', outall.shape)
        orderidx1 = np.argsort(outall[:, 1])
        if orderidx1.shape[0] > 0:
            #print("CHK", orderidx1.shape[0])
            #print(outall.shape[0], int(self.outall.shape[0]/2))
            lbs = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1]
            if orderidx1.shape[0] <= int(np.ceil(orderidx1.shape[0]*.84)):
                ubs = outall[orderidx1[-1], 1]
            else:
                ubs = outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.84))], 1]
            return outall.shape[0] < int(self.outall.shape[0]/2), (ubs - lbs) / 2.,\
                np.std(outall[:, 1]), np.average(outall[:, 1])
        else:
            return outall.shape[0] < int(self.outall.shape[0]/2), 0.07, 0.07, 0.07
         

    def readerror(self, hork, runNo):
        error = np.zeros((2, self.qsize))
        ave = np.zeros((self.qsize))
        mask = np.zeros((self.qsize), dtype=bool)
        #print(hork)
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

    def plotter(self, fig, nr, x, y, t, R2, e, mask, mask2, mask3, title, sidx,
                fidx):
        #xerr = 2.*(np.arange(self.qsize)*0.1 + 0.05)*0.05
        #mask *= e > 0.0001
        #mask *= t > 0.00001
        #print(mask, title, mask2, title)
        pnr = 1+len(self.fracs)*sidx+fidx
        ax = fig.add_subplot(nr, len(self.fracs), pnr)
        ax.plot(x, y*1000., ls='dashed', c='k')
        ax.errorbar(x[~mask*~mask2*~mask3], t[~mask*~mask2*~mask3]*1000.,
                    yerr=e[~mask*~mask2*~mask3]*1000., marker="x", ms=2, elinewidth=1,
                    lw=0, capsize=3, c='k')
        #ax.errorbar(x[mask*~mask2], t[mask*~mask2]*1000.,
        #            yerr=e[mask*~mask2]*1000., marker="x",
        #            ms=2, elinewidth=1, lw=0, capsize=3, c='blue')
        #ax.errorbar(x[mask2*~mask], t[mask2*~mask]*1000.,
        #            yerr=e[mask2*~mask]*1000., marker="x",
        #            ms=2, elinewidth=1, lw=0, capsize=3, c='red')
        #ax.errorbar(x[mask2*mask], t[mask2*mask]*1000.,
        #            yerr=e[mask2*mask]*1000., marker="x",
        #            ms=2, elinewidth=1, lw=0, capsize=3, c='pink')
        ax.errorbar(x[~(~mask3*~mask2*~mask)], t[~(~mask3*~mask2*~mask)]*1000.,
                    yerr=e[~(~mask3*~mask2*~mask)]*1000., marker="x",
                    ms=2, elinewidth=1, lw=0, capsize=3, c='lightgray')
        ax.text(0.1, 0.020*1000., title+"_"+str(self.fracs[fidx]))
        ax.text(1.6, 0.020*1000., '$R^2$={:.0f}%'.format(R2*100.))
        ax.text(0.1, 0.018*1000., 'D={:.0f} +/- {:.0f} '.format(self.D[sidx, fidx], self.stdD[sidx, fidx]) + '$10^9\AA^2s^{-1}$')
        ax.text(0.1, 0.016*1000., 'rel. diff. {:.0f}% '.format(self.D[sidx, fidx]/self.D[sidx, 0]*100.-100.))
        ax.text(0.1, 0.014*1000., 'τ={:.0f} +/- {:.0f} '.format(self.tau[sidx, fidx], self.stdtau[sidx, fidx]) + 'ps')
        ax.text(0.1, 0.012*1000., 'rel. diff. {:.0f}% '.format(self.tau[sidx, fidx]/self.tau[sidx, 0]*100.-100.))
        ax.set_ylim(-1, 0.022*1000.)
        #ax.set_xlim(0., 1.6)
        #ax.set_yticks([0.000, 0.005, 0.010, 0.015, 0.020])
        ax.set_yticks([0., 5, 10, 15, 20])
        #ax.set_xticks([0., 0.4, 0.8, 1.2, 1.6])
        if pnr >= (nr-1)*len(self.fracs)+1:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=True)
        else:
            ax.tick_params(direction='in', top=True, right=True,
                           labelbottom=False)
        #if pnr == 5 or pnr == 6:
        if sidx == 2:
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
        fig.suptitle("run#" +  self.runNo + ", " + self.rebin + ", " + self.nM)
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
        #print(gamma)
        #mask[0:2] = True
        #mask[0:2] = True
        mask[-2:] = True
        if self.nM.startswith("n8000"):
            mask3 = gamma < 0.0007
        elif self.nM == "n800":
            mask3 = gamma < 0.00015
        elif self.nM == "n200":
            mask3 = gamma < 0.00030
        #mask3 = gamma < 0.00015#for 0.00025
        #mask3 = gamma < 0.0007 #for norebin
        #mask3 = gamma < 0.00030 #for 0.001
        _mask = ~mask*~mask2*~mask3
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
                     error, mask, mask2, mask3, label, sidx, fidx)

    def eachrun(self, fidx, frac, fig):
        maskh, maskh2, gammah, maskkb, maskkb2, gammakb, maskk, maskk2,\
             gammak, errorh, errorkb,  errork, stdhes, aveh, avekb =\
             self.getdata(frac)
        #print(frac, len(maskh), len(gammah), len(maskkb), len(gammakb), len(maskk), len(gammak), len(errork[0]), len(stdhes))
        #maskhh = np.ones((self.qsize), dtype=bool)
        #if fidx >= 2:
        #    maskhh[1] = False
        #self.eachsolution(fig, 0, fidx, frac, gammah, errorh[0],
        self.eachsolution(fig, 0, fidx, frac, gammah, stdhes,
                          maskh, maskh2, 'hist')
        self.eachsolution(fig, 1, fidx, frac, gammakb, errorkb[0],
        #self.eachsolution(fig, 1, fidx, frac, gammakb, stdhes,
                          maskkb, maskkb2, 'kdeb')
        self.eachsolution(fig, 2, fidx, frac, gammak, errork[0],
                          maskk, maskk2, 'kde')


def testrun():
    runNo = "4174"
    qsize = 12
    fracs = ["100", "050", "025", "0125", "0100", "0050"]
    fracs = ["100", "050", "025", "0125"]
    fracs = ["100", "025",  "0125"]
    fracs = ["100", "050", "025",  "0125", "0100", "0050"]
    fracs = ["100", "0125"]
    fracsval = [1., 0.5, 0.25, 0.125, 0.1, 0.05]
    if len(sys.argv) >= 2:
        fracs = ["100", sys.argv[1]]

    #rebin = "rebin_hist_Ebin20150709"
    #rebin = "rebin_000025"
    #rebin = "norebin"
    rebin = "rebin_0001"
    #rebin = "rebin_0001_histe"
    #nM = "n200"
    nM = "n8000_M80_w5"
    #nM = "n200"
    prj = qens_model_fit(runNo, fracs, qsize, rebin, nM)
    prj.run()
    #plt.errorbar(fracsval, prj.D[0, :], yerr=prj.stdD[0, :], marker="x", ms=2, elinewidth=1, lw=0, capsize=3, label='hist')
    #plt.errorbar(fracsval, prj.D[1, :], yerr=prj.stdD[1, :], marker="x", ms=2, elinewidth=1, lw=0, capsize=3, label='kdeb')
    #plt.errorbar(fracsval, prj.D[2, :], yerr=prj.stdD[2, :], marker="x", ms=2, elinewidth=1, lw=0, capsize=3, label='kdeio')
    #plt.ylim([0, 120])
    #plt.legend()
    #plt.show()
    #plt.errorbar(fracsval, prj.tau[0, :], yerr=prj.stdtau[0, :], marker="x", ms=2, elinewidth=1, lw=0, capsize=3, label='hist')
    #plt.errorbar(fracsval, prj.tau[1, :], yerr=prj.stdtau[1, :], marker="x", ms=2, elinewidth=1, lw=0, capsize=3, label='kdeb')
    #plt.errorbar(fracsval, prj.tau[2, :], yerr=prj.stdtau[2, :], marker="x", ms=2, elinewidth=1, lw=0, capsize=3, label='kdeio')
    #plt.ylim([0, 120])
    #plt.legend()
    #plt.show()


testrun()
