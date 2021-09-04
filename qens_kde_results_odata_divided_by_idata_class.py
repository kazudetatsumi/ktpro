#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import pickle


class odata_divided_by_idata:
    def __init__(self, ofile, ifile):
        self.ofile = ofile
        self.ifile = ifile

    def read_pkl(self, pklfile):
        with open(pklfile, 'rb') as f:
            dataset = pickle.load(f)
        return(dataset)

    def save_data(self, rfile):
        dataset = {}
        dataset['xo'] = self.xo
        dataset['yr_ssvk'] = self.y_ssvk
        dataset['yr_ssk'] = self.y_ssvk
        with open(rfile, 'wb') as f:
            pickle.dump(dataset, f, -1)

    def interpolate(self, _yo, _yi):
        self.xo = np.linspace(self.oxlim[0], self.oxlim[1], _yo[1].shape[0])
        self.xi = np.linspace(self.ixlim[0], self.ixlim[1], _yi[1].shape[0])
        yo = _yo[0]
        yi = _yi[0]
        yi_ip = np.zeros_like(yo)
        for oidx, _xo in enumerate(self.xo):
            for iidx, _xi in enumerate(self.xi[:-1]):
                if _xo > _xi and _xo <= self.xi[iidx+1]:
                    yi_ip[oidx] = yi[iidx] + (yi[iidx+1] - yi[iidx]) /\
                                             (self.xi[iidx+1] - self.xi[iidx]) *\
                                             (_xo - _xi)
        return(yi_ip)

    def get_data(self):
        self.odataset = self.read_pkl(self.ofile)
        self.idataset = self.read_pkl(self.ifile)
        self.yo_ssvk = self.odataset['y_ssvk']
        self.yi_ssvk = self.idataset['y_ssvk']
        self.yo_ssk = self.odataset['y_ssk']
        self.yi_ssk = self.idataset['y_ssk']
        self.oxlim = self.odataset['xlim']
        self.ixlim = self.idataset['xlim']
        self.yi_ssvk_ip = self.interpolate(self.yo_ssvk, self.yi_ssvk)
        self.yi_ssk_ip = self.interpolate(self.yo_ssk, self.yi_ssk)
        self.y_ssvk = self.yo_ssvk[0] / self.yi_ssvk_ip
        self.y_ssk = self.yo_ssk[0] / self.yi_ssk_ip

    def plot_data(self):
        fig = plt.figure(figsize=(12, 12))
        fig.suptitle('Results of ssk and ssvk methods on a QENS data set')

        ax = self.ax_set(1, fig)
        ax.plot(self.xo, self.yo_ssk[0], label='outgoing_ssk')
        ax.plot(self.xo, self.yo_ssvk[0], label='outgoing_ssvk')
        plt.legend()

        ax = self.ax_set(3, fig)
        ax.plot(self.xi, self.yi_ssk[0], label='incoming_ssk')
        ax.plot(self.xi, self.yi_ssvk[0], label='incoming_ssvk')
        plt.legend()

        ax = self.ax_set(5, fig)
        plt.plot(self.xo, self.y_ssk, label='outgoing/incoming_ssk')
        plt.plot(self.xo, self.y_ssvk, label='outgoing/incoming_ssvk')
        ax.set_ylim(0.000, 5.00)
        plt.legend()

        ax = self.ax_set(2, fig)
        ax.plot(self.xo, self.yo_ssk[0], label='outgoing_ssk')
        ax.plot(self.xo, self.yo_ssvk[0], label='outgoing_ssvk')
        plt.legend()

        ax = self.ax_set(4, fig)
        ax.plot(self.xi, self.yi_ssk[0], label='incoming_ssk')
        ax.plot(self.xi, self.yi_ssvk[0], label='incoming_ssvk')
        plt.legend()

        ax = self.ax_set(6, fig)
        #ax.set_xlim(-0.06, 0.14)
        ax.set_ylim(0.0001, 7.00)
        plt.plot(self.xo, self.y_ssk, label='outgoing/incoming_ssk')
        plt.plot(self.xo, self.y_ssvk, label='outgoing/incoming_ssvk')
        plt.legend()
        plt.show()

    def ax_set(self, no, fig):
        ax = fig.add_subplot(3, 2, no)
        ax.set_xlim(self.oxlim[0], self.oxlim[1])
        ax.set_xlabel(r'Energy ($\mu$eV)')
        ax.set_ylabel(r'Intensity (arb. untis)')
        if no % 2 == 0:
            ax.set_yscale('log')
        return(ax)


def samplerun():
    #ofile = "/home/kazu/desktop/210108/Tatsumi/qens_kde_results_on_odata.pkl"
    #ifile = "/home/kazu/desktop/210108/Tatsumi/qens_kde_results_on_idata.pkl"
    #ofile = "/home/kazu/desktop/210108/Tatsumi/qens_wide_e_kde_results_on_odata.pkl"
    #ifile = "/home/kazu/desktop/210108/Tatsumi/qens_kde_results_on_idata_wide_E.pkl"
    ofile = "/home/kazu/desktop/210108/Tatsumi/qens_run6204_kde_results_on_odata.pkl"
    #ofile = "/home/kazu/desktop/210108/Tatsumi/qens_fine_run6204_kde_results_on_odata.pkl"
    ifile = "/home/kazu/desktop/210108/Tatsumi/qens_kde_results_on_idata_6204.pkl"
    #ofile = "/home/kazu/desktop/210108/Tatsumi/qens_run6207_kde_results_on_odata.pkl"
    #ifile = "/home/kazu/desktop/210108/Tatsumi/qens_kde_results_on_idata_6207.pkl"
    rfile = "/home/kazu/desktop/210108/Tatsumi/qens_kde_o_divided_by_i_6204.pkl"
    proj = odata_divided_by_idata(ofile, ifile)
    proj.get_data()
    proj.plot_data()
    proj.save_data(rfile)


#samplerun()
