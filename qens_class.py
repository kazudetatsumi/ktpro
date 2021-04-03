# Adaptive kernel density estimation on qens data:
# Presently, this script can prepare the input for the adaptive kernel density
# estimation from txt files containing a histogram per a detector element.
# These txt files are to be deduced from the raw data by a script written by
# BL02 staffs and stored in a directory.
# This script can apply the adaptive kernel density estimation on the prepared
# data and show the plots of the results, i.e., the estimated density and the
# estimated destribution of the kernel widths.
# This script uses the python methods of the adaptive kenerl denssity estimation
# distributed at the site linked by Dr. Shimazaki's  HP. 
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


class qens:
    def __init__(self, datadir, save_file):
        self.datadir = datadir
        self.save_file = save_file

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
                self.detector_position = np.zeros((len(self.filenames), 2))
            self.spectra[fidx, 0, :] = _tof
            self.spectra[fidx, 1, :] = _intensity
            self.spectra[fidx, 2, :] = _error
            self.detector_position[fidx, 0] = int(_psd)
            self.detector_position[fidx, 1] = int(_pix)
            f.close()

    def init_qens(self):
        self.get_data_from_container_txt()
        dataset = {}
        dataset['spectra'] = self.spectra
        dataset['detector_position'] = self.detector_position
        with open(self.save_file, 'wb') as f:
            pickle.dump(dataset, f, -1)
        print("Done!")

    def check_qens(self):
        if not os.path.exists(self.save_file):
            self.init_qens()
        with open(self.save_file, 'rb') as f:
            dataset = pickle.load(f)

        total_intensity = np.sum(dataset['spectra'][:, 1, :], axis=0)
        hw = dataset['spectra'][1010, 0, :]
        print(hw.shape)

        plt.plot(hw, np.log10(total_intensity))
        plt.xlabel('energy(micro eV)')
        plt.ylabel('log10(total count)')
        plt.show()

    def select_spectra(self):
        if not os.path.exists(self.save_file):
            self.init_qens()
        with open(self.save_file, 'rb') as f:
            dataset = pickle.load(f)

        spectra = dataset['spectra']
        dp = dataset['detector_position']
        mask = np.where((dp[:, 0] >= 10) & (dp[:, 1] <= 65))[0]
        self.selected_spectra = np.sum(spectra[mask, 1, :], axis=0)
        #np.savetxt('spectra.txt', self.selected_spectra, delimiter=',')
        self.xvec = np.array([idx for idx in
                             range(0, self.selected_spectra.shape[0]) for
                             num_repeat in
                             range(0, int(self.selected_spectra[idx]))
                              ], dtype=float)
        #np.savetxt('xvec.txt', self.xvec, delimiter=',')

    def add_shift(self):
        #print(self.xvec[0:30])
        self.xvec = self.xvec + np.random.uniform(-0.5, 0.5, size=self.xvec.shape[0])
        #print(self.xvec[0:30])

    def run_ssvkernel(self):
        testx = [4.37, 3.87, 4.00, 4.03, 3.50, 4.08, 2.25, 4.70, 1.73, 4.93, 1.73, 4.62,
             3.43, 4.25, 1.68, 3.92, 3.68, 3.10, 4.03, 1.77, 4.08, 1.75, 3.20, 1.85,
             4.62, 1.97, 4.50, 3.92, 4.35, 2.33, 3.83, 1.88, 4.60, 1.80, 4.73, 1.77,
             4.57, 1.85, 3.52, 4.00, 3.70, 3.72, 4.25, 3.58, 3.80, 3.77, 3.75, 2.50,
             4.50, 4.10, 3.70, 3.80, 3.43, 4.00, 2.27, 4.40, 4.05, 4.25, 3.33, 2.00,
             4.33, 2.93, 4.58, 1.90, 3.58, 3.73, 3.73, 1.82, 4.63, 3.50, 4.00, 3.67,
             1.67, 4.60, 1.67, 4.00, 1.80, 4.42, 1.90, 4.63, 2.93, 3.50, 1.97, 4.28,
             1.83, 4.13, 1.83, 4.65, 4.20, 3.93, 4.33, 1.83, 4.53, 2.03, 4.18, 4.43,
             4.07, 4.13, 3.95, 4.10, 2.72, 4.58, 1.90, 4.50, 1.95, 4.83, 4.12]

        #self.y = ssvkernel.ssvkernel(np.array(testx))
        self.y = ssvkernel.ssvkernel(self.xvec)
        self.y_ = sskernel.sskernel(self.xvec)
        #print(self.y[0].shape)
        #print(self.y[1].shape)
        #print(self.y[2].shape)
        #print(self.y_[2])
        norms = self.selected_spectra/np.sum(self.selected_spectra)
        ch = np.arange(0, norms.shape[0])
        fig = plt.figure(figsize=(12, 12))
        ax = fig.add_subplot(3, 1, 1)
        ax.bar(ch, norms, width=1.0, label='expt data')
        ax.plot(self.y[1], self.y[0], c='r', label='ssvkernel')
        ax.plot(self.y_[1], self.y_[0], c='k', label='sskernel')
        plt.legend()
        ax.tick_params(top=True, right=True, direction='in', which='both', labelbottom=False)
        ax.set_ylabel('density')
        ax.set_xlabel('# of energy channel')
        ax = fig.add_subplot(3, 1, 2)
        ax.set_ylabel('density')
        ax.set_xlabel('# of energy channel')
        ax.bar(ch, norms, width=1.0, label='expt data')
        ax.plot(self.y[1], self.y[0], c='r', label='ssvkernel')
        ax.plot(self.y_[1], self.y_[0], c='k', label='sskernel')
        ax.set_yscale('log')
        ax.tick_params(top=True, right=True, direction='in', which='both', labelbottom=False)
        ax = fig.add_subplot(3, 1, 3)
        ax.plot(self.y[1], self.y[2], c='r', label='ssvkernel')
        ax.plot(self.y_[1], np.ones(self.y_[1].shape[0])*self.y_[2], c='k', label='sskernel')
        ax.set_ylabel('band-width')
        ax.set_xlabel('# of energy channel')
        ax.tick_params(top=True, right=True, direction='in', which='both', labelbottom=True)
        plt.subplots_adjust(hspace=0.0)
        plt.legend()
        plt.show()


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
