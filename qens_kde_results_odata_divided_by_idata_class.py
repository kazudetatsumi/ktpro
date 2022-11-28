#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
sys.path.append("/home/kazu/ktpro")
import qens_class as qc


class odata_divided_by_idata(qc.qens):
    def __init__(self, ofile, ifile, iskde=True, denew=None, bootstrap=False):
        self.ofile = ofile
        self.ifile = ifile
        self.iskde = iskde
        self.denew = denew
        self.bootstrap = bootstrap

    def read_pkl(self, pklfile):
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
            #pickle.dump(dataset, f, -1)
            pickle.dump(dataset, f, 4)

    def get_data(self, norm=False):
        print(self.ofile, self.ifile)
        self.odataset = self.read_pkl(self.ofile)
        self.idataset = self.read_pkl(self.ifile)
        if self.iskde:
            self.yo_ssvk = self.odataset['y_ssvk']
            self.yi_ssvk = self.idataset['y_ssvk']
            self.yo_ssk = self.odataset['y_ssk']
            self.yi_ssk = self.idataset['y_ssk']
            self.xo = self.odataset['tin_real']
            self.xi = self.idataset['tin_real']
            if self.bootstrap:
                self.yi_ssvk_ip = np.zeros((self.yi_ssvk[6].shape[0],
                                            self.xo.shape[0]))
                for sidx, yi_ssvk_samp in enumerate(self.yi_ssvk[6]):
                    self.yi_ssvk_ip[sidx, :] = np.interp(self.xo, self.xi,
                                                         yi_ssvk_samp)
                self.y_ssvk = self.yo_ssvk[6] / self.yi_ssvk_ip
            else:
                self.yi_ssvk_ip = np.interp(self.xo, self.xi, self.yi_ssvk[0])
                self.yi_ssk_ip = np.interp(self.xo, self.xi, self.yi_ssk[0])
                self.y_ssvk = self.yo_ssvk[0] / self.yi_ssvk_ip
                self.y_ssk = self.yo_ssk[0] / self.yi_ssk_ip
        else:
            self.yo = self.odataset['spectra']
            self.xo = self.odataset['energy']
            self.yi = self.idataset['spectra'][0, 1, :]
            self.xi = self.idataset['spectra'][0, 0, :] - 2.085
            if norm:
                self.yo = self.yo/np.sum(self.yo)/(self.xo[1] - self.xo[0])
                mergin = 0.001
                xlim = np.array([-0.10025 - mergin, 0.14975 + mergin])
                mask = np.where((self.xi >= xlim[0]) & (self.xi <= xlim[1]))[0]
                self.selected_spectra = self.yi[mask]
                self.selected_energy = self.xi[mask]
                self.yi = self.selected_spectra /\
                    np.sum(self.selected_spectra) /\
                    (self.selected_energy[1] - self.selected_energy[0])
                self.xi = self.selected_energy
            if self.denew:
                self.reconstruct_hist()
            self.yi_ip = np.interp(self.xo, self.xi, self.yi)
            self.y = self.yo / self.yi_ip

    def reconstruct_hist(self):
        self.odata = True
        self.qsel = True
        self.dataset = self.odataset
        self.select_spectra()
        self.xo = np.linspace(self.selected_energy[0],
                              self.selected_energy[-1],
                              int(self.de/self.denew *
                              self.selected_energy.shape[0]))
        thist = np.concatenate((self.xo, (self.xo[-1]+self.denew)[np.newaxis]))
        self.yo = np.histogram(self.xvec_real, thist-self.denew/2)[0]
        
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

