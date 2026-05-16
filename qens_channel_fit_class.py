#!/usr/bin/env python
import numpy as np
import pickle
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.signal as ss
import re
import sys
sys.path.append("/home/kazu/ktpro")
import qens_fit_class as qfc


params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


class qens_fit(qfc.qens_fit):
    def __init__(self, devf, tf, elim, showplot=True):
        super().__init__(devf, tf, elim, showplot=showplot)

    def preprocess(self, doicorr=False):  ## ORG
        self.xo_devf, self.yr_ssvk_devf = self.get_data(self.devf)
        self.xo_tf, self.yr_ssvk_tf = self.get_data(self.tf)
        self.interpolate()
        if doicorr:
            x = self.x_tf + 2.085
            self.y_tf = self.y_tf / (self.k[0] + self.k[1]*x
                                     + self.k[2]*x**2
                                     + self.k[3]*x**3)
            self.y_df = self.y_df / (self.k[0] + self.k[1]*x
                                     + self.k[2]*x**2
                                     + self.k[3]*x**3)

    def preprocessh(self, doicorr=False):  ## ORG
        x_devf, ys_ssvk_devf = self.get_hdata(self.devf)
        x_tf, ys_ssvk_tf = self.get_hdata(self.tf)
        #print(np.sum(np.abs(x_tf - x_devf)))
        x_df, self.y_df = self.limit(x_devf, ys_ssvk_devf, mergin=0.00)
        self.x_tf, self.y_tf = self.limit(x_tf, ys_ssvk_tf, mergin=0.00)
        if doicorr:
            x = self.x_tf + 2.085
            self.y_tf = self.y_tf / (self.k[0] + self.k[1]*x
                                     + self.k[2]*x**2
                                     + self.k[3]*x**3)
            self.y_df = self.y_df / (self.k[0] + self.k[1]*x
                                     + self.k[2]*x**2
                                     + self.k[3]*x**3)
        self.bg = self.optbgpeakratio *\
            np.sum(self.y_tf[np.argmax(self.y_tf)-10:
                             np.argmax(self.y_tf)+10])

    def optimize(self, variables=[1.46103037e-04, 1.23754329e-02,
                 5.20429443e-01, 9.30889687e-06], figname=None):  ## ORG
        # initial guess on the parameters are hard-coded here.
        #variables = [0.00001, 0.0015, 0.01]
        #variables = [1.58837344e-04, 1.00454636e-02, 4.57573203e-01, 0.009]
        #variables = [1.38746043e-04, 8.27288080e-03, 4.47976536e-01, 1.75691683e-02]
        #variables = [1.46103037e-04, 1.23754329e-02, 5.20429443e-01, 9.30889687e-03]
        #variables = [1.38876225e-04, 8.09183272e-03, 4.42217308e-01]
        out = so.leastsq(self.res, variables,
                         args=(self.x_tf, self.y_df, self.y_tf), full_output=1,
                         epsfcn=0.0001)
        s_sq = (self.res(out[0], self.x_tf, self.y_df, self.y_tf)**2).sum() /\
               (len(self.y_tf)-len(out[0]))
        if not self.quiet:
            print("estimated constants alpha, gamma, delta, base")
            print(out[0])
        self.gamma = out[0][1]
        if not self.quiet:
            print("cov**0.5")
            print(np.absolute(out[1]*s_sq)**0.5)
        self.gammaerror = np.absolute(out[1][1][1]*s_sq)**0.5
        if len(variables) == 4:
            _alpha, _gamma, _delta, _base = out[0]
        elif len(variables) == 3:
            _alpha, _gamma, _delta = out[0]
            #_base = self.y_tf[-1]
            _base = self.bg
        #print("optbg/peak:", _base/np.max(self.y_tf))
        #self.bgpeakr = _base/np.max(self.y_tf)
        if not self.quiet:
            print("optbg/peak:", _base/np.sum(
                self.y_tf[np.argmax(self.y_tf)-10:np.argmax(self.y_tf)+10]),
              _base)
        self.optbgpeakratio = _base/np.sum(self.y_tf[np.argmax(self.y_tf)-10:np
                                           .argmax(self.y_tf)+10])
        if figname:
            _y = self.convlore(_alpha*self.y_df, _gamma, self.x_tf)
            _y += _delta*self.y_df + _base
            plt.plot(self.x_tf, self.y_tf, label='target')
            plt.plot(self.x_tf, np.zeros_like(self.x_tf) + _base,
                     label='constant')
            plt.plot(self.x_tf, _delta*self.y_df, label='delta_conv')
            plt.plot(self.x_tf, self.convlore(_alpha*self.y_df, _gamma,
                     self.x_tf), label='lore_conv')
            plt.plot(self.x_tf, _y, label='ML')
            plt.xlabel('energy (meV)')
            plt.tick_params(top=True, right=True, direction='in', which='both',
                            labelbottom=True, width=1.5)
            plt.yscale('log')
            plt.legend()
            plt.savefig(figname)
            if self.showplot:
                plt.show()
            plt.close()


def samplerun():
    head = "/home/kazu/desktop/210108/Tatsumi/pickles/"
# old pkls of spectra which were gathered over many detectors
    # devf = head + "qens_kde_o_divided_by_i_6204.pkl"
    # tf = head + "qens_kde_o_divided_by_i_6202.pkl"
# pkls of spectra over a specific q range, intensities were divided by incident
# neutron spectral profiles
    devf = head + "qens_kde_o_divided_by_i_6204_qsel.pkl"
    tf = head + "qens_kde_o_divided_by_i_6202_qsel.pkl"
    elim = np.array([-0.03, 0.10])
    proj = qens_fit(devf, tf, elim)
    proj.preprocess()
    proj.optimize()


def ssamplerun():
    head = "/home/kazu/desktop/210108/Tatsumi/pickles/"
# pkls of spectra over a specific q range, intensities were corrected by a
# standard method.
    devf = head + "qens_run6204_kde_results_on_sdata_qsel.pkl"
    tf = head + "qens_run6202_kde_results_on_sdata_qsel.pkl"
    elim = np.array([-0.03, 0.10])
    proj = qens_fit(devf, tf, elim)
    proj.preprocesss()
    proj.optimize()

#samplerun()
