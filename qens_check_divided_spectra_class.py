#!/usr/bin/env python
import matplotlib.pylab as plt
import sys
sys.path.append("/home/kazu/ktpro")
import qens_kde_results_odata_divided_by_idata_class as qkrodic
import qens_fit_class as qfc


class check_divided_spectra(qkrodic.odata_divided_by_idata, qfc.qens_fit):
    def __init__(self, histofile, histifile, dividedkdefile):
        self.ofile = histofile
        self.ifile = histifile
        self.dividedkdefile = dividedkdefile

        # get divided hist data
        self.iskde = False
        qkrodic.odata_divided_by_idata.get_data(self, norm=True)
        # get divided kde data
        self.x_kde, self.y_kde = qfc.qens_fit.get_data(self, self.dividedkdefile)

    def plotter(self):
        fig = plt.figure(figsize=(6, 8))
        ax = fig.add_subplot(3, 1, 1)
        ax.bar(self.xo, self.y, width=self.xo[1]-self.xo[0], label='expt data')
        ax.plot(self.x_kde, self.y_kde, c='r', label='ssvkernel')
        ax.set_ylim(0, 6.0)
        ax.set_xlim(-0.05, 0.125)
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False, width=1.5)
        ax.set_ylabel('density')
        plt.legend()
        ax = fig.add_subplot(3, 1, 2)
        ax.set_ylabel('density')
        ax.set_xlabel('energy (meV)')
        ax.bar(self.xo, self.y, width=self.xo[1]-self.xo[0], label='expt data')
        ax.plot(self.x_kde, self.y_kde, c='r', label='ssvkernel')
        ax.set_yscale('log')
        ax.set_ylim(0.03, 6.0)
        ax.set_xlim(-0.05, 0.125)
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=True, width=1.5)
        plt.subplots_adjust(hspace=0.0)
        plt.legend()
        plt.savefig('divided_spectra.png')
        plt.show()


def samplerun():
    histofile = "/home/kazu/desktop/210108/Tatsumi/srlz/0000025io/run7784united_spectra.pkl"
    histifile = "/home/kazu/desktop/210108/Tatsumi/srlz/0000025io/run7784united_monispectra.pkl"
    dividedkdefile = "/home/kazu/desktop/210108/Tatsumi/winparam_exam_7784/160_0.5_0000025io_Boxcar/qens_kde_o_divided_by_i_7784.pkl"
    proj = check_divided_spectra(histofile, histifile, dividedkdefile)
    proj.plotter()


#samplerun()
