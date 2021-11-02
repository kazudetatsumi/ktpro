#!/usr/bin/env python
import numpy as np
import pickle
import matplotlib.pyplot as plt
import scipy.optimize as so
import re

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


class qens_fit:
    def __init__(self, devf, tf, elim, plot=True):
        self.devf = devf
        self.tf = tf
        self.elim = elim
        self.plot = plot
        self.quiet = False

    def preprocess(self, doicorr=False):
        x_devf, y_ssvk_devf = self.get_data(self.devf)
        x_tf, y_ssvk_tf = self.get_data(self.tf)
        #self.interpolate()
        x_df, self.y_df = self.limit(x_devf, y_ssvk_devf, mergin=0.00)
        self.x_tf, self.y_tf = self.limit(x_tf, y_ssvk_tf, mergin=0.00)
        #self.x_df, self.y_df = self.get_data(self.devf)
        #self.x_tf, self.y_tf = self.get_data(self.tf)
        print(x_df[0:10])
        print(self.x_tf[0:10])
        print(self.x_tf.shape)
        if doicorr:
            x = self.x_tf + 2.085
            self.y_tf = self.y_tf / (self.k[0] + self.k[1]*x
                                     + self.k[2]*x**2
                                     + self.k[3]*x**3)
            self.y_df = self.y_df / (self.k[0] + self.k[1]*x
                                     + self.k[2]*x**2
                                     + self.k[3]*x**3)


    def preprocesss(self, doicorr=False):
        x_devf, ys_ssvk_devf = self.get_sdata(self.devf)
        x_tf, ys_ssvk_tf = self.get_sdata(self.tf)
        self.bg = 0.000000
        #print(np.sum(np.abs(x_tf - x_devf)))
        x_df, self.y_df = self.limit(x_devf, ys_ssvk_devf[0], mergin=0.00)
        self.x_tf, self.y_tf = self.limit(x_tf, ys_ssvk_tf[0], mergin=0.00)
        #print("y_tf[-1]", self.y_tf[-1])
        if doicorr:
            x = self.x_tf + 2.085
            self.y_tf = self.y_tf / (self.k[0] + self.k[1]*x
                                     + self.k[2]*x**2
                                     + self.k[3]*x**3)
            self.y_df = self.y_df / (self.k[0] + self.k[1]*x
                                     + self.k[2]*x**2
                                     + self.k[3]*x**3)

    def preprocessh(self, doicorr=False):
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
        self.bg = self.optbgpeakratio *np.sum(self.y_tf)
            #np.sum(self.y_tf[np.argmax(self.y_tf)-100:np.argmax(self.y_tf)+100])

    def get_data(self, infile):
        with open(infile, 'rb') as f:
            data = pickle.load(f, encoding='latin1')
        return data['xo'], data['yr_ssvk']

    def get_sdata(self, infile):
        with open(infile, 'rb') as f:
            data = pickle.load(f, encoding='latin1')
        return data['tin_real'], data['ys_ssvk']

    def get_hdata(self, infile):
        with open(infile, 'rb') as f:
            data = pickle.load(f, encoding='latin1')
        return data['energy'], data['spectra']

    def limit(self, x, y, mergin=0.0):
        mask = np.where((x > self.elim[0] - mergin) &
                        (x < self.elim[1] + mergin))
        return x[mask], y[mask]

    def interpolate(self):
        x1, y1 = self.limit(self.xo_devf, self.yr_ssvk_devf, mergin=0.01)
        x2, y2 = self.limit(self.xo_tf, self.yr_ssvk_tf, mergin=0.01)

        y2_ip = np.zeros_like(y1)
        for idx1, _x1 in enumerate(x1):
            for idx2, _x2 in enumerate(x2[:-1]):
                if _x1 > _x2 and _x1 <= x2[idx2+1]:
                    y2_ip[idx1] = y2[idx2] + (y2[idx2+1] - y2[idx2]) /\
                                             (x2[idx2+1] - x2[idx2]) *\
                                             (_x1 - _x2)
        self.x_df, self.y_df = self.limit(x1, y1)
        self.x_tf, self.y_tf = self.limit(x1, y2_ip)

    def checkdata(self):
        plt.plot(self.x_df, self.y_df)
        plt.plot(self.x_tf, self.y_tf)
        if not self.quiet:
            print(self.x_df.shape)
            print(self.x_tf.shape)
        plt.show()

    def fun_lore(self, x, gamma):
        return 1/np.pi*gamma/(x**2 + gamma**2)

    def convlore(self, f, gamma, x):
        _convd = np.zeros_like(x)
        for _x, _f, in zip(x, f):
            _convd += self.fun_lore(x - _x, gamma)*_f
        return _convd

    def res(self, coeffs, x, d, t):
        if len(coeffs) == 4:
            [alpha, gamma, delta, base] = coeffs
            y = self.convlore(alpha*d, gamma, x)
            y += delta*d + base
        if len(coeffs) == 3:
            #print(self.y_tf[-1])
            #print(self.bg)
            [alpha, gamma, delta] = coeffs
            y = self.convlore(alpha*d, gamma, x)
            #y += delta*d + self.y_tf[-1]
            y += delta*d + self.bg
        return t - y

    def res_icorr(self, coeffs, x, t):
        [k0, k1, k2, k3] = coeffs
        y = k0 + k1*x + k2*x**2+ k3*x**3
        return t - y

    def optimize(self, variables=[1.46103037e-04, 1.23754329e-02, 5.20429443e-01, 9.30889687e-06]):
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
        #self.optbgpeakratio = _base/np.sum(self.y_tf[np.argmax(self.y_tf)-100:np
        #                                   .argmax(self.y_tf)+100])
        self.optbgpeakratio = _base/np.sum(self.y_tf)
        if self.plot:
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
            plt.savefig('qens_fit.pdf')
            plt.show()

    def get_icorrdata(self, icorrfile):
        x = []
        y = []
        for line in open(icorrfile):
            if not re.compile('[A-z[]').search(line):
                x.append(float(line.split()[0]))
                y.append(float(line.split()[1]))
        x = np.array(x)
        y = np.array(y)
        return x, y

    def icorr(self):
        x, y = self.get_icorrdata("/home/kazu/desktop/210108/Tatsumi/" +
                                   "Ilambda_correction.txt")
        variables = [10, 10, 10, 10]
        out = so.leastsq(self.res_icorr, variables,
                         args=(x, y), full_output=1,
                         epsfcn=0.0001)
        
        #s_sq = (self.res_icorr(out[0], x, y)**2).sum()/(len(y)-len(out[0]))
        #print(s_sq)
        #print(out[0])
        #print(out[1])
        #print("cov**0.5")
        #print(np.absolute(out[1]*s_sq)**0.5)
        [_k0, _k1, _k2, _k3] = out[0]
        #_x = np.linspace(np.min(x), np.max(x), 200)
        #_y = _k0 + _k1*_x + _k2*_x**2 + _k3*_x**3
        #plt.plot(x, y, marker='o', lw=0)
        #plt.plot(_x, _y)
        #plt.show()
        self.k = out[0]

    def check_spectra(self):
        plt.plot(self.x_tf, self.y_tf, label=self.tf)
        plt.plot(self.x_tf, self.y_df, label=self.devf)
        plt.yscale('log')
        plt.legend()
        plt.show()

    def kde_hist_sub(self, tf, devf, kde=True, variables=None):
        self.devf = devf
        self.tf = tf
        if not self.quiet:
            print(self.devf)
            print(self.tf)
        if kde:
            self.preprocesss(doicorr=True)
        else:
            self.preprocessh(doicorr=True)
        if variables:
            self.optimize(variables=variables)
        else:
            self.optimize()

    def kde_hist(self, kvariables=None,
                 hvariables=[1.46103037e-04, 1.23754329e-02, 5.20429443e-01]):
        self.icorr()
        if not self.quiet:
            print("entering kde part")
        self.kde_hist_sub(self.ktf, self.kdevf, kde=True, variables=kvariables)
        kgamma = self.gamma
        kgammaerror = self.gammaerror
        if not self.quiet:
            print("entering histogram part")
        self.kde_hist_sub(self.htf, self.hdevf, kde=False, variables=hvariables)
        hgamma = self.gamma
        hgammaerror = self.gammaerror
        return(kgamma, kgammaerror, hgamma, hgammaerror)


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
