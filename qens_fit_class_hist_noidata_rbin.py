#!/usr/bin/env python
import os
import matplotlib.pyplot as plt
import numpy as np
import pickle
import scipy.optimize as so
import copy
import sys
sys.path.append("/home/kazu/ktpro")
from qens_fit_class import qens_fit as qf


class runhistnoidata(qf):
    def __init__(self, devf, tf, outfile, alpha, elim, elimw, numcycle=100,
                 leastsq=True):
        self.elim = elim
        self.elimw = elimw
        self.devf = devf
        self.tf = tf
        self.outfile = outfile
        self.alpha = alpha
        self.elim = elim
        self.numcycle = numcycle
        self.leastsq = leastsq

    def get_xmlyd(self):
        x, yd, yt = self.preprocess()
        # normalization
        dx = x[1] - x[0]
        yd = yd / np.sum(yd) / dx
        yt = yt / np.sum(yt) / dx
        _out = self.optimize(x, yd, yt,
#                            variables=[2.18704786e-04, 1.67980295e-02,
#                                       4.92405238e-05, 1.88866588e-03,
#                                       1.21127501e-01, 5.02759930e-02])
                            #variables=[0.59704786e-00, 2.67980295e-02,
                            #           3.82405238e-01, 7.88866588e-03,
                            #           0.21127501e+00, 1.82759930e-02])
                            variables=[0.54e-00, 2.7e-02,
                                       3.5e-01, 7.0e-03,
                                       0.18e+00, 1.9e-02])
        print(_out[0])
        self.ml = self.reconstruct(x, yd, _out[0])
        self.yd = yd
        self.x = x

    def preprocess(self):
        self.icorr()
        xd, yd = self.get_data(self.devf)
        xt, yt = self.get_data(self.tf)
        xdl, ydl = self.limit2(xd, yd, self.elimw)
        xtl, ytl = self.limit2(xt, yt, self.elimw)
        ydlc, ytlc = self.correction(xtl, ydl, ytl)
        return xtl, ydlc, ytlc

    def cycle(self):
        self.outall = []
        for cyidx in range(0, self.numcycle):
            simd, simt = self.generate_data()
            # normalization
            # dx = self.x[1] - self.x[0]
            # simd = simd / np.sum(simd) / dx
            # simt = simt / np.sum(simt) / dx * 100.
            # rebinning
            tin_real, simdr = self.rebin_generated_samples(self.x, simd, num=480)
            tin_real, simtr = self.rebin_generated_samples(self.x, simt, num=480)
            dx = tin_real[1] - tin_real[0]
            simdr = simdr / np.sum(simdr) / dx
            simtr = simtr / np.sum(simtr) / dx * 100.
            _out = self.optimize(tin_real, simdr, simtr,
            # out = self.optimize(self.x, simd, simt,
#                                 variables=[1.73704786e-05, 2.66580295e-02,
#                                           9.96405238e-06, 7.00766588e-03,
#                                           2.00077501e-01, 1.78759930e-01])
                                #variables=[0.18704786e-00, 2.67980295e-02,
                                #           1.02405238e-01, 6.48866588e-03,
                                #           0.06127501e+00, 0.13759930e-00])
                                #variables=[0.7e-00, 2.6e-02,
                                #           2.3e-01, 7.0e-03,
                                #           0.17e+00, 0.38e-00])
                                variables=[6.2e+01, 2.7e-02,
                                           2.5e+01, 7.0e-03,
                                           2.0e+01, 4.0e+01])
            self.check_out(cyidx, _out)

    def check_out(self, cyidx, _out):
        if self.leastsq:
            if _out[1] is None:
                print(cyidx, 'curveture is flat. omitting..')
            else:
                #self.modify_out(cyidx, _out[0])
                self.modify_out(cyidx, _out)
        else:
            if _out[1]:
                #self.modify_out(cyidx, _out[0])
                self.modify_out(cyidx, _out)
            else:
                print(cyidx, 'optimization is not converged..')

    def modify_out(self, cyidx, out):
        #if len(out) >= 5:
        if len(out[0]) >= 5:
            #if out[0] < 0 and out[1] < 0:
            if out[0][0] < 0 and out[0][1] < 0:

                #print("negative-negative")
                #out[0] = out[0]*(-1.)
                #out[1] = out[1]*(-1.)
                out[0][0] = out[0][0]*(-1.)
                out[0][1] = out[0][1]*(-1.)
            #if out[2] < 0 and out[3] < 0:
            if out[0][2] < 0 and out[0][3] < 0:
                #print("negative-negative")
                #out[2] = out[2]*(-1.)
                #out[3] = out[3]*(-1.)
                out[0][2] = out[0][2]*(-1.)
                out[0][3] = out[0][3]*(-1.)
            #if out[1] < out[3]:
            if out[0][1] < out[0][3]:
                #print("exchange")
                #tmpout = out[1]
                #tmpout2 = out[0]
                #out[1] = out[3]
                #out[3] = tmpout
                #out[0] = out[2]
                #out[2] = tmpout2
                tmpout = out[0][1]
                tmpout2 = out[0][0]
                out[0][1] = out[0][3]
                out[0][3] = tmpout
                out[0][0] = out[0][2]
                out[0][2] = tmpout2
                tmpout = out[2][1]
                tmpout2 = out[2][0]
                out[2][1] = out[2][3]
                out[2][3] = tmpout
                out[2][0] = out[2][2]
                out[2][2] = tmpout2
        # if self.rank:
        #     if self.rank == 0:
        #         print(cyidx, out)
        # else:
        #     print(cyidx, out)
        #self.outall.append(out)
        #conbine the lists of out.x, out.active_mask, out.success
        print('chk2', out[0].tolist() + out[2].tolist())
        self.outall.append(out[0].tolist() + out[2].tolist())

    def output(self):
        self.outall = np.array(self.outall)
        if self.outall.shape[0] == 1:
            print(self.outall)
        elif self.outall.shape[1] == 6 or self.outall.shape[1] == 12:
            orderidx1 = np.argsort(self.outall[:, 1])
            print("median of gamma1:", np.median(self.outall[:, 1]))
            print("average of gamma1:", np.average(self.outall[:, 1]))
            print("68% CI of gamma1")
            print(self.outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1])
            print(self.outall[orderidx1[int(np.floor(orderidx1.shape[0]*.84))], 1])
            orderidx2 = np.argsort(self.outall[:, 3])
            print("median of gamma2:", np.median(self.outall[:, 3]))
            print("average of gamma2:", np.average(self.outall[:, 3]))
            print("68% CI of gamma2")
            print(self.outall[orderidx2[int(np.ceil(orderidx2.shape[0]*.16))], 3])
            print(self.outall[orderidx2[int(np.floor(orderidx2.shape[0]*.84))], 3])
            ave1 = np.average(self.outall[:, 1])
            std1 = np.std(self.outall[:, 1])
            ave2 = np.average(self.outall[:, 3])
            std2 = np.std(self.outall[:, 3])
            mask = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                            & (self.outall[:, 2] > 0) & (self.outall[:, 3] > 0)
                            & (self.outall[:, 4] > 0) & (self.outall[:, 5] > 0))
            maskwobg = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                                & (self.outall[:, 2] > 0) & (self.outall[:, 3] > 0)
                                & (self.outall[:, 4] > 0))
            outnonneg = self.outall[mask]
            avenonneg1 = np.average(outnonneg[:, 1])
            stdnonneg1 = np.std(outnonneg[:, 1])
            outnonnegwobg = self.outall[maskwobg]
            avenonnegwobg1 = np.average(outnonnegwobg[:, 1])
            stdnonnegwobg1 = np.std(outnonnegwobg[:, 1])
            avenonneg2 = np.average(outnonneg[:, 3])
            stdnonneg2 = np.std(outnonneg[:, 3])
            avenonnegwobg2 = np.average(outnonnegwobg[:, 3])
            stdnonnegwobg2 = np.std(outnonnegwobg[:, 3])
            print('ave_gamma1 std_gamma1 ave_gamma2 std_gamma2 # '
                  'ave_gamma1 std_gamma1 ave_gamma2 std_gamma2 # '
                  'ave_gamma1 std_gamma1 ave_gamma2 std_gamma2 # ')
            print('{0:.8e} {1:.8e} {2:.8e} {3:.8e} {4} '
                  '{5:.8e} {6:.8e} {7:.8e} {8:.8e} {9} '
                  '{10:.8e} {11:.8e} {12:.8e} {13:.8e} {14}'
                  .format(avenonneg1, stdnonneg1, avenonneg2, stdnonneg2,
                          outnonneg.shape[0],
                          avenonnegwobg1, stdnonnegwobg1, avenonnegwobg2,
                          stdnonnegwobg2, outnonnegwobg.shape[0],
                          ave1, std1, ave2, std2, self.outall.shape[0]))
        elif self.outall.shape[1] == 4 or self.outall.shape[1] == 8:
            orderidx1 = np.argsort(self.outall[:, 1])
            #print(self.outall)
            print("median of gamma1:", np.median(self.outall[:, 1]))
            print("average of gamma1:", np.average(self.outall[:, 1]))
            print("68% CI of gamma1")
            print(self.outall[orderidx1[int(np.ceil(orderidx1.shape[0]*.16))], 1])
            print(self.outall[orderidx1[int(np.floor(orderidx1.shape[0]*.84))], 1])
            #orderidx2 = np.argsort(self.outall[:, 3])
            #print("median of gamma2:", np.median(self.outall[:, 3]))
            #print("average of gamma2:", np.average(self.outall[:, 3]))
            #print("68% CI of gamma2")
            #print(self.outall[orderidx2[int(np.ceil(orderidx2.shape[0]*.16))], 3])
            #print(self.outall[orderidx2[int(np.floor(orderidx2.shape[0]*.84))], 3])
            ave1 = np.average(self.outall[:, 1])
            std1 = np.std(self.outall[:, 1])
            #ave2 = np.average(self.outall[:, 3])
            #std2 = np.std(self.outall[:, 3])
            mask = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                            & (self.outall[:, 2] > 0) & (self.outall[:, 3] > 0))
            maskwobg = np.where((self.outall[:, 0] > 0) & (self.outall[:, 1] > 0)
                                & (self.outall[:, 2] > 0))
            outnonneg = self.outall[mask]
            avenonneg1 = np.average(outnonneg[:, 1])
            stdnonneg1 = np.std(outnonneg[:, 1])
            outnonnegwobg = self.outall[maskwobg]
            avenonnegwobg1 = np.average(outnonnegwobg[:, 1])
            stdnonnegwobg1 = np.std(outnonnegwobg[:, 1])
            #avenonneg2 = np.average(outnonneg[:, 3])
            #stdnonneg2 = np.std(outnonneg[:, 3])
            #avenonnegwobg2 = np.average(outnonnegwobg[:, 3])
            #stdnonnegwobg2 = np.std(outnonnegwobg[:, 3])
            print('ave_gamma1 std_gamma1 # '
                  'ave_gamma1 std_gamma1 # '
                  'ave_gamma1 std_gamma1 # ')
            print('{0:.8e} {1:.8e} {2} '
                  '{3:.8e} {4:.8e} {5} '
                  '{6:.8e} {7:.8e} {8}'
                  .format(avenonneg1, stdnonneg1, outnonneg.shape[0],
                          avenonnegwobg1, stdnonnegwobg1, outnonnegwobg.shape[0],
                          ave1, std1, self.outall.shape[0]))

    def savefile(self):
        dataset = {}
        dataset['out'] = self.outall
        with open(self.outfile, 'wb') as f:
            #pickle.dump(dataset, f, -1)
            pickle.dump(dataset, f, 4)

    def loadfile(self):
        with open(self.outfile, 'rb') as f:
            dataset = pickle.load(f, encoding='latin1')
        self.outall = dataset['out']

    def plot_distribution(self, binwidth1, binwidth2):
        numsumple = self.outall[:, 1].shape[0]
        numbins1 = int(np.ceil((np.max(self.outall[:, 1]) -
                               np.min(self.outall[:, 1]))/binwidth1))
        numbins2 = int(np.ceil((np.max(self.outall[:, 3]) -
                               np.min(self.outall[:, 3]))/binwidth2))
        heights1, bins1 = np.histogram(self.outall[:, 1], bins=numbins1)
        heights2, bins2 = np.histogram(self.outall[:, 3], bins=numbins2)
        ## plot only the solutions not on the boundaries.
        #heights1, bins1 = np.histogram(self.outall[np.sum(self.outall[:, 6:], axis=1) > -0.1, 1], bins=numbins1)
        #heights2, bins2 = np.histogram(self.outall[np.sum(self.outall[:, 6:], axis=1) > -0.1, 3], bins=numbins2)
        plt.bar(bins1[:-1]+binwidth1/2., heights1/binwidth1/numsumple,
                width=binwidth1, label='$\Gamma_1$')
        plt.bar(bins2[:-1]+binwidth2/2., heights2/binwidth2/numsumple,
                width=binwidth2, label='$\Gamma_2$')
        plt.xlabel('HWHM (meV)')
        plt.ylabel('Distribution (1/meV)')
        plt.xlim(0, 0.08)
        plt.legend()
        plt.show()

    def plot_distribution_single(self, binwidth1, show=True):
        numsumple = self.outall[:, 1].shape[0]
        numbins1 = int(np.ceil((np.max(self.outall[:, 1]) -
                               np.min(self.outall[:, 1]))/binwidth1))
        heights1, bins1 = np.histogram(self.outall[:, 1], bins=numbins1)
        plt.bar(bins1[:-1]+binwidth1/2., heights1/binwidth1/numsumple,
                width=binwidth1, label='$\Gamma_1$', color='r')
                #width=binwidth1, label='$\Gamma_1$', alpha=0.5, ec=None)
        plt.xlabel('HWHM (meV)')
        plt.ylabel('Distribution (1/meV)')
        plt.xlim(0.000, 0.016)
        plt.legend()
        plt.tick_params(top=True, right=True, direction='in', which='both',
                        labelbottom=True, width=1.5)
        if show:
            plt.show()

    def correction(self, x, yd, yt):
        x = x + 2.085
        yd = yd / (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        yt = yt / (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        return yd, yt

    def decorrection(self, x, yd, yt):
        x = x + 2.085
        yd = yd * (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        yt = yt * (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        return yd, yt

    def optimize(self, x, yd, yt,
                 variables=[6.e-6, 2.e-2, 1.e-6, 4.e-3, 7.e-3, 3.e-1]):
        # leastsq
        if self.leastsq:
            out = so.leastsq(self.res, variables, args=(x, yd, yt),
                             full_output=1, epsfcn=0.0001)
            return out
        # least_squares
        else:
            bounds = (0, np.inf)
            out = so.least_squares(self.res, variables, bounds=bounds,
                                   args=(x, yd, yt))
            #if self.rank == 0:
            print(out.active_mask, out.success, out.x)
            #out = so.least_squares(self.res, variables, args=(x, yd, yt))
            _out = [out.x, np.linalg.inv(np.dot(out.jac.T, out.jac))]
            s_sq = (self.res(_out[0], x, yd, yt)**2).sum() / (len(yt)-len(_out[0]))
            print("cov**0.5")
            print(np.absolute(_out[1]*s_sq)**0.5)
            return [out.x, out.success, out.active_mask]

    def res(self, coeffs, x, d, t):
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
        xl, dif = self.limit2(x, (t-y)/t, self.elim)
        return dif

    def reconstruct(self, x, yd, out):
        _alpha, _gamma, _alpha2, _gamma2,  _delta, _base = out
        #_base = _base * 9.
        return _alpha*self.convlore(yd, _gamma, x)\
            + _alpha2*self.convlore(yd, _gamma2, x)\
            + _delta*yd + _base

#    def limit(self, x, y, elim):
#        mask = np.where((x > elim[0]) & (x < elim[1]))
#        return x[mask], y[mask]

    def generate_data(self, idata=False):
        if idata:
            return np.random.poisson(self.iyd/np.sum(self.iyd)
                                     * 94612.*self.alpha)*1.,\
                   np.random.poisson(self.iyt/np.sum(self.iyt)
                                     * 50299.*self.alpha)*1.
        else:
            return np.random.poisson(self.yd/np.sum(self.yd)
                                     * 59146.*self.alpha)*1.,\
                   np.random.poisson(self.ml/np.sum(self.ml)
                                     * 18944.*self.alpha)*1.


def testrun():
    outfile = "./outhistnoidata.pkl"
    alpha = 0.5
    devf = "./qens_kde_o_divided_by_i_6204.pkl"
    tf = "./qens_kde_o_divided_by_i_6202.pkl"
    elim = [-0.03, 0.07]
    elimw = [-0.04, 0.08]
    numcycle = 3000
    binwidth1=0.0016
    binwdith2=0.00074
    binwdith2=0.0016
    np.set_printoptions(linewidth=120)
    if os.path.isfile(outfile):
        proj = runhistnoidata(devf, tf, outfile, alpha, elim, elimw,
                              leastsq=True, numcycle=numcycle)
        proj.loadfile()
        proj.output()
        proj.plot_distribution(binwidth1, binwdith2)
    else:
        np.random.seed(314)
        proj = runhistnoidata(devf, tf, outfile, alpha, elim, elimw,
                              leastsq=True, numcycle=numcycle)
        proj.get_xmlyd()
        proj.cycle()
        proj.output()
        proj.savefile()


#testrun()
