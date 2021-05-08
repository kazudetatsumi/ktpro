# Adaptive kernel density estimation on qens data:
# Presently, this script can prepare the input for the adaptive kernel density
# estimation from txt files containing a histogram per a detector element.
# These txt files are to be deduced from the raw data by a script written by
# BL02 staffs and stored in a directory.
# This script can apply the adaptive kernel density estimation on the prepared
# data and show the plots of the results, i.e., the estimated density and the
# estimated destribution of the kernel widths.
# This script uses the python methods of the adaptive kenerl denssity
# estimation distributed at the site linked by Dr. Shimazaki's  HP.
# The adaptive kenel density estimation was developped by Dr. Shimazaki.
# Kazuyoshi TATSUMI 2021/04/03
import numpy as np
import os
import re
import pickle
import matplotlib.pylab as plt
import sys
sys.path.append("/home/kazu/desktop/210108/AdaptiveKDE/adaptivekde")
import ssvkernel, sskernel, sshist

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


class qens:
    def __init__(self, datadir, save_file, odata=True):
        self.datadir = datadir
        self.save_file = save_file
        self.odata = odata
        if not os.path.exists(self.save_file):
            self.init_qens()
        with open(self.save_file, 'rb') as f:
            self.dataset = pickle.load(f)

    def get_filenames(self):
        self.filenames = os.popen('/bin/ls ' + self.datadir + self.string)\
                         .read().split()

    def get_data_from_container_txt(self):
        self.string = '*.txt'
        self.get_filenames()
        for fidx, fn in enumerate(self.filenames):
            _psd = fn[len(self.datadir)+11:len(self.datadir)+14]
            _pix = fn[len(self.datadir)+18:len(self.datadir)+21]
            _tof = []
            _intensity = []
            _error = []
            f = open(fn)
            for line in f:
                if not re.compile(r'[A-z]').search(line):
                    values = line.split(',')
                    _tof.append(float(values[0]))
                    _intensity.append(float(values[1]))
                    _error.append(float(values[2]))
            if fidx == 0:
                self.spectra = np.zeros((len(self.filenames), 3, len(_tof)))
                if self.odata:
                    self.detector_position = np.zeros((len(self.filenames), 2))
            self.spectra[fidx, 0, :] = _tof
            self.spectra[fidx, 1, :] = _intensity
            self.spectra[fidx, 2, :] = _error
            if self.odata:
                self.detector_position[fidx, 0] = int(_psd)
                self.detector_position[fidx, 1] = int(_pix)
            f.close()

    def init_qens(self):
        self.get_data_from_container_txt()
        dataset = {}
        dataset['spectra'] = self.spectra
        if self.odata:
            dataset['detector_position'] = self.detector_position
        with open(self.save_file, 'wb') as f:
            pickle.dump(dataset, f, -1)
        print("Done!")

    def check_qens(self):
        total_intensity = np.sum(self.dataset['spectra'][:, 1, :], axis=0)
        hw = self.dataset['spectra'][1010, 0, :]
        print(hw.shape)

        plt.plot(hw, np.log10(total_intensity))
        plt.xlabel('energy(micro eV)')
        plt.ylabel('log10(total count)')
        plt.show()

    def select_spectra(self):
        spectra = self.dataset['spectra']
        if self.odata:
            dp = self.dataset['detector_position']
            mask = np.where((dp[:, 0] >= 10) & (dp[:, 1] <= 65))[0]
            self.selected_spectra = np.sum(spectra[mask, 1, :], axis=0)
            self.selected_energy = spectra[0, 0, :]
        else:
            spectra[0, 0, :] = spectra[0, 0, :] - 2.085
            mergin = 0.001
            xlim = np.array([-0.10025 - mergin, 0.14975 + mergin])
            mask = np.where((spectra[0, 0, :] >= xlim[0]) &
                            (spectra[0, 0, :] <= xlim[1]))[0]
            self.selected_spectra = spectra[0, 1, mask]
            self.selected_energy = spectra[0, 0, mask]

        self.de = spectra[0, 0, 1] - spectra[0, 0, 0]

        # To keep the accuracy, kde is executed on the channel numbers in the
        # histogram data.
        # The actual energies were retrieved by "_real" variables.
        self.xvec = np.array([idx for idx in
                             range(0, self.selected_spectra.shape[0]) for
                             num_repeat in
                             range(0, int(self.selected_spectra[idx]))
                              ], dtype=float)
        self.xvec_real = np.array([self.selected_energy[idx] for idx in
                                  range(0, self.selected_spectra.shape[0]) for
                                  num_repeat in
                                  range(0, int(self.selected_spectra[idx]))
                                   ], dtype=float)

    def add_shift(self):
        self.shift = np.random.uniform(-0.5, 0.5, size=self.xvec.shape[0])
        self.xvec += self.shift
        self.xvec_real += self.shift*self.de

        #print(self.xvec[0:30])

    def run_ssvkernel(self):
        #self.y = ssvkernel.ssvkernel(np.array(testx))
        T = (np.max(self.xvec) - np.min(self.xvec))
        dx = np.sort(np.diff(np.sort(self.xvec)))
        dt_samp = dx[np.nonzero(dx)][0]
        # To relieve the effect of the finite counts near the upper and lower
        # bounds of the measurement range
        # we extend the range by some fraction of the original range
        #mergin = T*0.2
        #tin = np.linspace(np.min(self.xvec) - mergin, np.max(self.xvec) + mergin,
        #                  int(min(np.ceil(T*1.4 / dt_samp), 1e3)))
        # tin is the space where densities are inferred.
        tin = np.linspace(np.min(self.xvec), np.max(self.xvec),
                          int(min(np.ceil(T / dt_samp), 1e3)))
        tin_real = np.linspace(np.min(self.xvec_real), np.max(self.xvec_real),
                               int(min(np.ceil(T / dt_samp), 1e3)))

        self.y = ssvkernel.ssvkernel(self.xvec, tin)
        self.y_ = sskernel.sskernel(self.xvec, tin)
        norms = self.selected_spectra/np.sum(self.selected_spectra)/self.de
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(3, 1, 1)
        ax.bar(self.selected_energy, norms, width=self.de, label='expt data')
        ax.plot(tin_real, self.y[0]/self.de, c='r', label='ssvkernel')
        ax.plot(tin_real, self.y_[0]/self.de, c='k', label='sskernel')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False)
        ax.set_ylabel('density')
        ax.set_xlabel('energy ($\mu eV$)')
        plt.legend()
        ax = fig.add_subplot(3, 1, 2)
        ax.set_ylabel('density')
        ax.bar(self.selected_energy, norms, width=self.de, label='expt data')
        ax.plot(tin_real, self.y[0]/self.de, c='r', label='ssvkernel')
        ax.plot(tin_real, self.y_[0]/self.de, c='k', label='sskernel')
        ax.set_yscale('log')
        ax.set_ylim(0.0001, np.max(self.y_[0])/self.de)
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False)
        plt.legend()
        ax = fig.add_subplot(3, 1, 3)
        ax.plot(tin_real, self.y[2], c='r', label='ssvkernel')
        ax.plot(tin_real, np.ones(self.y_[1].shape[0])*self.y_[2], c='k',
                label='sskernel')
        ax.set_ylabel('band-width')
        ax.set_xlabel('energy ($\mu eV$)')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=True)
        plt.subplots_adjust(hspace=0.0)
        plt.legend()
        plt.show()

    def save_output(self,  output_file):
        dataset = {}
        dataset['y_ssvk'] = self.y
        dataset['y_ssk'] = self.y_
        dataset['xlim'] = np.array([np.min(self.xvec_real),
                                   np.max(self.xvec_real)])
        with open(output_file, 'wb') as f:
            pickle.dump(dataset, f, -1)


def samplerun():
    datadir = "/home/kazu/desktop/210108/Tatsumi/run6202_containers/"
    #dataset_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_dir = datadir + "../"
    save_file = dataset_dir + "/qens.pkl"
    proj = qens(datadir, save_file)
    proj.select_spectra()
    proj.add_shift()
    proj.run_ssvkernel()


#samplerun()
