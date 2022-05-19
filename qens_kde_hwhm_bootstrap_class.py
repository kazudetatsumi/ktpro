#!/usr/bin/env python
# hwhms calculated with boot-strap sampled kde profiles
# Kazuyoshi Tatsumi 2022/04/07
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
home=os.path.expanduser("~")
from qens_kde_results_odata_divided_by_idata_class \
 import odata_divided_by_idata as div
from qens_fit_class import qens_fit as fit


class bootstrap(div, fit):
    def __init__(self, odfile, idfile, ofile, ifile, elim):
        self.bootodfile = odfile
        self.bootidfile = idfile
        self.bootofile = ofile
        self.bootifile = ifile
        self.elim = elim

    def divide(self):
        div.__init__(self, self.bootodfile, self.bootidfile, bootstrap=True)
        div.get_data(self)
        self.y_df = self.y_ssvk
        self.x_df = self.xo
        div.__init__(self, self.bootofile, self.bootifile, bootstrap=True)
        div.get_data(self)
        #for i in range(0,10):
        #   plt.plot(self.y_ssvk[i], lw=0.4)
        #   plt.yscale('log')
        #plt.show()
        self.y_tf = self.y_ssvk
        self.x_tf = self.xo

    # override the method in fit class.
    def preprocess(self, doiocorr=True):
        self.x_df, self.y_df = self.limit(self.x_df, np.transpose(self.y_df))
        self.y_df = np.transpose(self.y_df)
        self.x_tf, self.y_tf = self.limit(self.x_tf, np.transpose(self.y_tf))
        self.y_tf = np.transpose(self.y_tf)
        if doiocorr:
            fit.icorr(self)
            fit.correction(self)

    def allfit(self):
        self.y_dfs = self.y_df
        self.y_tfs = self.y_tf
        self.allout = np.zeros((self.y_dfs.shape[0], 7, 6))
        allout = []
        for sidx, (self.y_df, self.y_tf) in enumerate(zip(self.y_dfs, self.y_tfs)):
            fit.optimize(self, variables=[1.62856562e-05, 3.05453407e-02,
                                          1.11413842e-05, 7.41258597e-03,
                                          2.07570182e-01, 2.27120576e-02])
            if np.all(self.out > 0.):
                if self.out[1] < self.out[3]:
                    tmp = self.out[1]
                    self.out[1] = self.out[3]
                    self.out[3] = tmp
                allout.append(self.out)
        allout = np.array(allout)
        print(allout.shape)
        print(np.average(allout[:, 1]), np.std(allout[:, 1]))
        print(np.average(allout[:, 3]), np.std(allout[:, 3]))
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(1, 2, 1)
        ax.hist(allout[:, 1], bins='auto')
        ax = fig.add_subplot(1, 2, 2)
        ax.hist(allout[:, 3], bins='auto')
        plt.show()
        return np.array([np.average(allout[:, 1]), np.average(allout[:, 3])]),\
            np.array([np.std(allout[:, 1]), np.std(allout[:, 3])])


def samplerun():
    fracs = ["", "075", "0625", "05", "0375", "025", "0125"]
    fracs = [""]
    elim = [-0.03, 0.07]
    hwhms = np.zeros((2, len(fracs)))
    shwhms = np.zeros((2, len(fracs)))
    tcount = np.array([18970., 16098., 13380., 10621., 7794., 5318., 2826.])
    x = tcount/tcount[0]
    ylims = [[0., 0.11], [0., 0.012]]
    savefile = "./results_bootstrap.pkl"
    if not os.path.exists(savefile):
        for fidx, frac in enumerate(fracs):
            print(frac)
            odfile = home+"/desktop/210108/Tatsumi/winparam_exam/" +\
                "test_python_de_0.000025/160_1_0000025io_Boxcar_/" +\
                "qens_run6204united_kde_results_on_data_qsel.pkl"
            idfile = home+"/desktop/210108/Tatsumi/pickles/" +\
                "0000025io/qens_kde_results_on_idata_6204.pkl"
            ofile = home+"/desktop/210108/Tatsumi/winparam_exam/" +\
                "test_python_de_0.000025/160_1_0000025io_Boxcar_"+frac+"/" +\
                "qens_run6202"+frac+"united_kde_results_on_data_qsel.pkl"
            ifile = home+"/desktop/210108/Tatsumi/pickles/" +\
                "0000025io/qens_kde_results_on_idata_6202"+frac+".pkl"
            prj = bootstrap(odfile, idfile, ofile, ifile, elim)
            prj.quiet = True
            prj.divide()
            prj.preprocess()
            hwhms[:, fidx], shwhms[:, fidx] = prj.allfit()
        dataset = {}
        dataset['odfile'] = odfile
        dataset['idfile'] = idfile
        dataset['ofile'] = ofile
        dataset['ifile'] = ifile
        dataset['hwhms'] = hwhms
        dataset['shwhms'] = shwhms
        dataset['x'] = x
        with open("./results_bootstrap.pkl", 'wb') as f:
            #pickle.dump(dataset, f, -1)
            pickle.dump(dataset, f, 4)
    else:
        with open(savefile, 'rb') as f:
            #print('reading ', savefile)
            dataset = pickle.load(f)
            hwhms = dataset['hwhms']
            shwhms = dataset['shwhms']
            x = dataset['x']
    print(dataset['odfile'])
    print(dataset['idfile'])
    print(dataset['ofile'])
    print(dataset['ifile'])
    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_subplot(2, 3, 1)
    ax.errorbar(x, hwhms[0, :], yerr=shwhms[0, :], marker="o", label='kde1',
                ms=4, elinewidth=1, lw=0, capsize=3, c='lightblue', mfc="none")
    ax.set_ylim(ylims[0][0], ylims[0][1])

    ax = fig.add_subplot(2, 3, 4)
    ax.errorbar(x, hwhms[1, :], yerr=shwhms[1, :], marker="o", label='kde1',
                ms=4, elinewidth=1, lw=0, capsize=3, c='blue', mfc="none")
    ax.set_ylim(ylims[1][0], ylims[1][1])
    plt.subplots_adjust(hspace=0.0)
    plt.show()


#samplerun()
