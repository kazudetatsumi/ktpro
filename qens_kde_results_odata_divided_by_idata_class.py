#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import pickle


class odata_divided_by_idata:
    def __init__(self, ofile, ifile, iskde=True):
        self.ofile = ofile
        self.ifile = ifile
        self.iskde = iskde

    def read_pkl(self, pklfile):
        print(pklfile)
        with open(pklfile, 'rb') as f:
            if self.iskde:
                dataset = pickle.load(f)
            else:
                dataset = pickle.load(f, encoding='latin1')
        return(dataset)

    def save_data(self, rfile):
        dataset = {}
        if self.iskde:
            dataset['xo'] = self.xo
            dataset['yr_ssvk'] = self.y_ssvk
            dataset['yr_ssk'] = self.y_ssvk
        else:
            dataset['energy'] = self.xo
            dataset['spectra'] = self.y
        with open(rfile, 'wb') as f:
            pickle.dump(dataset, f, -1)

    #def interpolate(self, yi):
        #if self.iskde:
            #self.xo = np.linspace(self.oxlim[0], self.oxlim[1], _yo[1].shape[0])
            #self.xi = np.linspace(self.ixlim[0], self.ixlim[1], _yi[1].shape[0])
        #    yi = _yi[0]
        #else:
        #    yi = _yi
        #return(np.interp(self.xo, self.xi, yi))

    def get_data(self):
        self.odataset = self.read_pkl(self.ofile)
        self.idataset = self.read_pkl(self.ifile)
        if self.iskde:
            self.yo_ssvk = self.odataset['y_ssvk']
            self.yi_ssvk = self.idataset['y_ssvk']
            self.yo_ssk = self.odataset['y_ssk']
            self.yi_ssk = self.idataset['y_ssk']
            self.xo = self.odataset['tin_real']
            self.xi = self.idataset['tin_real']
            #self.oxlim = self.odataset['xlim']
            #self.ixlim = self.idataset['xlim']
            #self.yi_ssvk_ip = self.interpolate(self.yo_ssvk, self.yi_ssvk)
            #self.yi_ssk_ip = self.interpolate(self.yo_ssk, self.yi_ssk)
            self.yi_ssvk_ip = np.interp(self.xo, self.xi, self.yi_ssvk[0])
            self.yi_ssk_ip = np.interp(self.xo, self.xi, self.yi_ssk[0])
            self.y_ssvk = self.yo_ssvk[0] / self.yi_ssvk_ip
            self.y_ssk = self.yo_ssk[0] / self.yi_ssk_ip
        else:
            self.yo = self.odataset['spectra']
            self.xo = self.odataset['energy']
            self.yi = self.idataset['spectra'][0, 1, :]
            self.xi = self.idataset['spectra'][0, 0, :] - 2.085
            #self.yi_ip = self.interpolate(self.yo, self.yi)
            self.yi_ip = np.interp(self.xo, self.xi, self.yi)
            self.y = self.yo / self.yi_ip


    def plot_data(self):
        fig = plt.figure(figsize=(12, 12))
        if self.iskde:
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
        else:
            fig.suptitle('A QENS histogram data set')
            ax = self.ax_set(1, fig)
            ax.plot(self.xo, self.yo, label='outgoing')
            plt.legend()
            ax = self.ax_set(3, fig)
            ax.plot(self.xi, self.yi, label='incoming')
            plt.legend()
            ax = self.ax_set(5, fig)
            plt.plot(self.xo, self.y, label='outgoing/incoming')
            ax.set_ylim(0.000, 5.00)
            plt.legend()
            ax = self.ax_set(2, fig)
            ax.plot(self.xo, self.yo, label='outgoing')
            plt.legend()
            ax = self.ax_set(4, fig)
            ax.plot(self.xi, self.yi, label='incoming')
            plt.legend()
            ax = self.ax_set(6, fig)
            ax.set_ylim(0.0001, 7.00)
            plt.plot(self.xo, self.y, label='outgoing/incoming')
        plt.legend()
        plt.show()

    def ax_set(self, no, fig):
        ax = fig.add_subplot(3, 2, no)
        #if self.iskde:
        #    ax.set_xlim(self.oxlim[0], self.oxlim[1])
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

