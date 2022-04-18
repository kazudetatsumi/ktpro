#!/usr/bin/env python
import numpy as np
import pickle
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy.signal as ss
import re

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


class qens_fit:
    def __init__(self, devf, tf, elim, showplot=True):
        self.devf = devf
        self.tf = tf
        self.elim = elim
        self.showplot = showplot
        self.quiet = False

    def preprocessnoi(self, doicorr=False):
        x_devf, y_ssvk_devf = self.get_idata(self.devf)
        x_tf, y_ssvk_tf = self.get_idata(self.tf)
        x_df, self.y_df = self.limit(x_devf, y_ssvk_devf, mergin=0.00)
        self.x_tf, self.y_tf = self.limit(x_tf, y_ssvk_tf, mergin=0.00)

    def preprocess(self, doicorr=False):
        x_devf, y_ssvk_devf = self.get_data(self.devf)
        x_tf, y_ssvk_tf = self.get_data(self.tf)
        #self.interpolate()
        x_df, self.y_df = self.limit(x_devf, y_ssvk_devf, mergin=0.00)
        self.x_tf, self.y_tf = self.limit(x_tf, y_ssvk_tf, mergin=0.00)
        #self.x_df, self.y_df = self.get_data(self.devf)
        #self.x_tf, self.y_tf = self.get_data(self.tf)
 #temp       print("check x_df and x_tf")
 #temp       print(x_df[0:5], x_df[-5:])
 #temp       print(self.x_tf[0:5], self.x_tf[-5:])
        if doicorr:
            self.correction()
            #x = self.x_tf + 2.085
            #self.y_tf = self.y_tf / (self.k[0] + self.k[1]*x
            #                         + self.k[2]*x**2
            #                         + self.k[3]*x**3)
            #self.y_df = self.y_df / (self.k[0] + self.k[1]*x
            #                         + self.k[2]*x**2
            #                         + self.k[3]*x**3)
    
    def correction(self):
        x = self.x_tf + 2.085
        self.y_tf /= (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3
                      )
        self.y_df /= (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3
                      )

    def decorrection(self):
        x = self.x_tf + 2.085
        self.ml *= (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3)
        self.y_df *= (self.k[0] + self.k[1]*x + self.k[2]*x**2 + self.k[3]*x**3
                      )

    def preprocesss(self, doicorr=False):
        x_devf, ys_ssvk_devf = self.get_sdata(self.devf)
        x_tf, ys_ssvk_tf = self.get_sdata(self.tf)
        self.bg = 0.000000
        #print(np.sum(np.abs(x_tf - x_devf)))
        x_df, self.y_df = self.limit(x_devf, ys_ssvk_devf[0], mergin=0.00)
        self.x_tf, self.y_tf = self.limit(x_tf, ys_ssvk_tf[0], mergin=0.00)
        #print("y_tf[-1]", self.y_tf[-1])
        if doicorr:
            self.correction()
            #x = self.x_tf + 2.085
            #self.y_tf = self.y_tf / (self.k[0] + self.k[1]*x
            #                         + self.k[2]*x**2
            #                         + self.k[3]*x**3)
            #self.y_df = self.y_df / (self.k[0] + self.k[1]*x
            #                         + self.k[2]*x**2
            #                         + self.k[3]*x**3)

    def preprocessh(self, doicorr=False):
        x_devf, ys_ssvk_devf = self.get_hdata(self.devf)
        x_tf, ys_ssvk_tf = self.get_hdata(self.tf)
        #print(np.sum(np.abs(x_tf - x_devf)))
        x_df, self.y_df = self.limit(x_devf, ys_ssvk_devf, mergin=0.00)
        self.x_tf, self.y_tf = self.limit(x_tf, ys_ssvk_tf, mergin=0.00)
        if doicorr:
            self.correction()
            #x = self.x_tf + 2.085
            #self.y_tf = self.y_tf / (self.k[0] + self.k[1]*x
            #                         + self.k[2]*x**2
            #                         + self.k[3]*x**3)
            #self.y_df = self.y_df / (self.k[0] + self.k[1]*x
            #                         + self.k[2]*x**2
            #                         + self.k[3]*x**3)
        if 'optbgpeakratio' in dir(self):
            self.bg = self.optbgpeakratio*np.sum(self.y_tf)*(x_tf[1]-x_tf[0])
            #np.sum(self.y_tf[np.argmax(self.y_tf)-100:np.argmax(self.y_tf)+100])

    def get_idata(self, infile):
        with open(infile, 'rb') as f:
            data = pickle.load(f, encoding='latin1')
        return data['tin_real'], data['y_ssvk'][0]

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

    def convloreorg(self, f, gamma, x):
        _convd = np.zeros_like(x)
        for _x, _f, in zip(x, f):
            _convd += self.fun_lore(x - _x, gamma)*_f
        return _convd

    def convlore(self, f, gamma, x):
        ex = np.linspace(np.min(x)-(np.max(x)-np.min(x))*0.5,
                         np.max(x)+(np.max(x)-np.min(x))*0.5, x.shape[0]*2+1)
        win = self.fun_lore(ex - ex[int(x.shape[0])], gamma)
        _convd = ss.convolve(f, win, mode='same', method='fft')
        return _convd

    def testconv(self):
        plt.plot(self.convloreorg(self.y_df, 0.001, self.x_tf), lw=10)
        plt.plot(self.convlore(self.y_df, 0.001, self.x_tf))
        plt.yscale('log')
        plt.show()

    def res(self, coeffs, x, d, t):
        if len(coeffs) == 6:
            [alpha1, gamma1, alpha2, gamma2,  delta, base] = coeffs
            #y = self.convlore(alpha1*d, gamma1, x)
            #y += self.convlore(alpha2*d, gamma2, x)
            #y += delta*d + base
            y = alpha1*self.convlore(d, gamma1, x)
            y += alpha2*self.convlore(d, gamma2, x)
            y += delta*d + base
        if len(coeffs) == 5:
            [alpha1, gamma1, alpha2, gamma2,  delta] = coeffs
            y = self.convlore(alpha1*d, gamma1, x)
            y += self.convlore(alpha2*d, gamma2, x)
            y += delta*d + self.bg
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
        y = k0 + k1*x + k2*x**2 + k3*x**3
        return t - y

    def optimize(self, variables=[1.46103037e-04, 1.23754329e-02,
                 5.20429443e-01, 9.30889687e-06], figname=None):
        # initial guess on the parameters are hard-coded here.
        #variables = [0.00001, 0.0015, 0.01]
        #variables = [1.58837344e-04, 1.00454636e-02, 4.57573203e-01, 0.009]
        #variables = [1.38746043e-04, 8.27288080e-03, 4.47976536e-01, 1.75691683e-02]
        #variables = [1.46103037e-04, 1.23754329e-02, 5.20429443e-01, 9.30889687e-03]
        #variables = [1.38876225e-04, 8.09183272e-03, 4.42217308e-01]
        #print(self.x_tf.shape)
        out = so.leastsq(self.res, variables,
                         args=(self.x_tf, self.y_df, self.y_tf), full_output=1,
                         epsfcn=0.0001)
        s_sq = (self.res(out[0], self.x_tf, self.y_df, self.y_tf)**2).sum() /\
               (len(self.y_tf)-len(out[0]))
        if not self.quiet:
            if len(variables) == 4:
                print("estimated constants alpha, gamma, delta, base")
            if len(variables) == 3:
                print("estimated constants alpha, gamma, delta")
            if len(variables) == 5:
                print("estimated constants alpha1, gamma1, alpha2, gamma2, delta")
            #if len(variables) == 6:
                #print("estimated constants alpha1, gamma1, alpha2, gamma2, delta, base")
            print(out[0])
        self.out = out[0]
        self.gamma = out[0][1]
   #temp     if not self.quiet:
   #temp         print("cov**0.5")
   #temp         #print(out[1])
   #temp         print(np.absolute(out[1]*s_sq)**0.5)
        #self.gammaerror = np.absolute(out[1][1][1]*s_sq)**0.5
        if len(variables) == 4:
            _alpha, _gamma, _delta, _base = out[0]
        elif len(variables) == 3:
            _alpha, _gamma, _delta = out[0]
            #_base = self.y_tf[-1]
            _base = self.bg
        elif len(variables) == 6:
            _alpha, _gamma, _alpha2, _gamma2,  _delta, _base = out[0]
        elif len(variables) == 5:
            _alpha, _gamma, _alpha2, _gamma2,  _delta = out[0]
            _base = self.bg
        #print("optbg/peak:", _base/np.max(self.y_tf))
        #self.bgpeakr = _base/np.max(self.y_tf)
        #if not self.quiet:
        #    print("optbg/peak:", _base/np.sum(
        #        self.y_tf[np.argmax(self.y_tf)-10:np.argmax(self.y_tf)+10]),
        #      _base)
        #self.optbgpeakratio = _base/np.sum(self.y_tf[np.argmax(self.y_tf)-100:np
        #                                   .argmax(self.y_tf)+100])
        self.optbgpeakratio = _base/np.sum(self.y_tf)\
            / (self.x_tf[1]-self.x_tf[0])
        #if not self.quiet:
        #    print("optbg/peak:", self.optbgpeakratio)
        if figname:
            #_y = self.convlore(_alpha*self.y_df, _gamma, self.x_tf)
            _y = _alpha*self.convlore(self.y_df, _gamma, self.x_tf)
            _y += _delta*self.y_df + _base
            if len(variables) >= 5:
                #_y += self.convlore(_alpha2*self.y_df, _gamma2, self.x_tf)
                _y += _alpha2*self.convlore(self.y_df, _gamma2, self.x_tf)
            plt.plot(self.x_tf, self.y_tf, label='target')
            plt.plot(self.x_tf, np.zeros_like(self.x_tf) + _base,
                     label='constant')
            plt.plot(self.x_tf, _delta*self.y_df, label='delta_conv')
            #plt.plot(self.x_tf, self.convlore(_alpha*self.y_df, _gamma,
            #         self.x_tf), label='lore_conv')
            plt.plot(self.x_tf, _alpha*self.convlore(self.y_df, _gamma,
                     self.x_tf), label='lore_conv')
            plt.plot(self.x_tf, _y, label='ML')
            if len(variables) >= 5:
                #plt.plot(self.x_tf, self.convlore(_alpha2*self.y_df, _gamma2,
                #         self.x_tf), label='lore_conv')
                plt.plot(self.x_tf, _alpha2*self.convlore(self.y_df, _gamma2,
                         self.x_tf), label='lore_conv')
            plt.xlabel('energy (meV)')
            plt.tick_params(top=True, right=True, direction='in', which='both',
                            labelbottom=True, width=1.5)
            plt.yscale('log')
            plt.legend()
            plt.savefig(figname)
            if self.showplot:
                plt.show()
            plt.close()

    def reconstruct(self, elim=None, check=True, idevf=None, itf=None):
        if elim:
            self.elim = elim
        x_devf, y_ssvk_devf = self.get_data(self.devf)
        x_tf, y_ssvk_tf = self.get_data(self.tf)
        x_df, self.y_df = self.limit(x_devf, y_ssvk_devf, mergin=0.00)
        self.x_tf, self.y_tf = self.limit(x_tf, y_ssvk_tf, mergin=0.00)
        self.correction()
        _alpha, _gamma, _alpha2, _gamma2,  _delta, _base = self.out
        _y = self.convlore(_alpha*self.y_df, _gamma, self.x_tf)
        _y += _delta*self.y_df + _base
        _y += self.convlore(_alpha2*self.y_df, _gamma2, self.x_tf)
        self.ml = _y
        if check:
            plt.plot(self.x_tf, self.ml)
            plt.plot(self.x_tf, self.y_tf)
            plt.yscale('log')
            plt.show()
        if idevf and itf:
            self.decorrection()
            self.multii(idevf, itf)

    def multii(self, idevf, itf):
        xid, yid = self.get_idata(idevf)
        xit, yit = self.get_idata(itf)
        xid, yid = self.limit(xid, yid, mergin=0.00)
        xit, yit = self.limit(xit, yit, mergin=0.00)
        yid_ip = np.interp(self.x_tf, xid, yid)
        yit_ip = np.interp(self.x_tf, xit, yit)
        self.y_df *= yid_ip
        self.ml *= yit_ip

    def generate_data(self, idevf, itf, check=True, rebin=False):
        self.y_df = self.y_df/np.sum(self.y_df)*59146.*10
        self.ml = self.ml/np.sum(self.ml)*18944.*10
        ddata = np.random.poisson(self.y_df)*1.
        tdata = np.random.poisson(self.ml)*1.
        #xid, yid = self.get_idata(idevf)
        #xit, yit = self.get_idata(itf)
        #yid = yid/np.sum(yid)*103203
        #yit = yit/np.sum(yit)*54678
        #iddata = np.random.poisson(yid)
        #itdata = np.random.poisson(yit)
        if rebin:
            tin_real, ddata = self.rebin_generated_samples(self.x_tf, ddata)
            tin_real, tdata = self.rebin_generated_samples(self.x_tf, tdata)
        else:
            tin_real = self.x_tf
        if check:
            self.check_generated_samples(self.x_tf, ddata)
            self.check_generated_samples(self.x_tf, tdata)
            #self.check_generated_samples(xid, iddata)
            #self.check_generated_samples(xit, itdata)
        self.save_generated_data(tin_real, ddata, 'qens_sim_6204.pkl')
        self.save_generated_data(tin_real, tdata, 'qens_sim_6202.pkl')
        #self.save_generated_data(xid, iddata, 'qens_sim_moni_6204.pkl')
        #self.save_generated_data(xit, itdata, 'qens_sim_moni_6202.pkl')

    def save_generated_data(self, x, data, savefile):
        dataset = {}
        dataset['energy'] = x
        dataset['spectra'] = data
        with open(savefile, 'wb') as f:
            pickle.dump(dataset, f, -1)

    def rebin_generated_samples(self, x, data,  num=600, shift=False):
        dx = x[1]-x[0]
        xvec_real = np.array([x[idx] for idx in range(0, data.shape[0]) for
                              num_repeat in range(0, int(data[idx]))],
                             dtype=float)
        if shift:
            xvec_real += np.random.uniform(0., 1.0, size=xvec_real.shape[0])*dx
        tin_real = np.linspace(x[0], x[-1], num=num)
        dt = tin_real[1] - tin_real[0]
        print(dt)
        thist = np.concatenate((tin_real, (tin_real[-1]+dt)[np.newaxis]))
        return tin_real, np.histogram(xvec_real, thist-dt/2)[0]*1.

    def check_generated_samples(self, x, data):
        dx = x[1]-x[0]
        xvec_real = np.array([x[idx] for idx in range(0, data.shape[0]) for
                              num_repeat in range(0, int(data[idx]))],
                             dtype=float)
        xvec_real += np.random.uniform(0., 1.0, size=xvec_real.shape[0])*dx
        tin_real = np.linspace(x[0], x[-1], num=600)
        dt = tin_real[1] - tin_real[0]
        thist = np.concatenate((tin_real, (tin_real[-1]+dt)[np.newaxis]))
        y_hist = np.histogram(xvec_real, thist-dt/2)[0]
        plt.bar(tin_real, y_hist, width=dt)
        plt.yscale('log')
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

    def save_result(self):
        dataset = {}
        dataset['fitparam'] = self.out
        dataset['x'] = self.x_tf

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

    def save_generated_data(self, x, data, savefile):
        dataset = {}
        dataset['energy'] = x
        dataset['spectra'] = data
        with open(savefile, 'wb') as f:
            pickle.dump(dataset, f, -1)

    def check_generated_samples(self, x, data):
        #dx = x[1]-x[0]
        xvec_real = np.array([x[idx] for idx in range(0, data.shape[0]) for
                              num_repeat in range(0, int(data[idx]))],
                             dtype=float)
        #xvec_real += np.random.uniform(0., 1.0, size=xvec_real.shape[0])*dx
        tin_real = np.linspace(x[0], x[-1], num=800)
        dt = tin_real[1] - tin_real[0]
        thist = np.concatenate((tin_real, (tin_real[-1]+dt)[np.newaxis]))
        y_hist = np.histogram(xvec_real, thist-dt/2)[0]
        plt.bar(tin_real, y_hist, width=dt)
        plt.yscale('log')
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

    def save_result(self):
        dataset = {}
        dataset['fitparam'] = self.out
        dataset['x'] = self.x_tf

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
