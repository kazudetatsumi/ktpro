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
## if you use fotran library of ssvk or ssk, uncomment the corresponding lines
## and comment out the import lines of the  python library.
#import ssvkernel_fort as ssvkernel
#import sskernel_fort as sskernel
import ssvkernel
import sskernel
import ssvkernel_ktpro as svk
import sskernel_ktpro as sk

params = {'mathtext.default': 'regular', 'axes.linewidth': 1.5}
plt.rcParams.update(params)


class qens:
    def __init__(self, datadir, save_file, odata=True, sfile=None, qsel=False,
                 optsm=False, winparam=1, M=80, WinFunc='Boxcar',
                 figname='qens_out.png', showplot=True):
        self.datadir = datadir
        self.save_file = save_file
        self.odata = odata
        self.qsel = qsel               # in def self.select_energy()
        self.optsm = optsm             # in def run_ssvkernel()
        self.winparam = winparam
        self.M = M
        self.WinFunc = WinFunc
        self.figname = figname
        self.showplot = showplot
        if not os.path.exists(self.save_file):
            print(self.save_file, "is not found. Entering init_qens()")
            self.init_qens()
        with open(self.save_file, 'rb') as f:
            self.dataset = pickle.load(f, encoding='latin1')
        if sfile:
            with open(sfile, 'rb') as f:
                dataset = pickle.load(f, encoding='latin1')
                self.senergy = dataset['energy']
                self.sspectra = dataset['spectra']
                self.sde = self.senergy[1] - self.senergy[0]

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
        if self.odata:  # case outgoing beam
            if self.qsel:   # case spectra are already integrated over a
                            # specific q range.
                self.selected_spectra = self.dataset['spectra']
                self.selected_energy = self.dataset['energy']
            else:       # case dataset are distributed over 2-D PSD elements.
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
        #print(self.selected_energy[np.argmax(self.selected_spectra)])
        self.de = self.selected_energy[1] - self.selected_energy[0]
        print("check!!!", self.selected_energy[0], self.selected_energy[-1],self.de)
        self.get_xvec()

    #def get_xvec(self):
        # To keep the accuracy, kde is executed on the channel numbers in the
        # histogram data.
        # The actual energies were retrieved by "_real" variables.
        #self.xvec = np.array([idx for idx in
        #                     range(0, self.selected_spectra.shape[0]) for
        #                     num_repeat in
        #                     range(0, int(self.selected_spectra[idx]))
        #                      ], dtype=float)
        #self.xvec_real = np.array([self.selected_energy[idx] for idx in
        #                          range(0, self.selected_spectra.shape[0]) for
        #                          num_repeat in
        #                          range(0, int(self.selected_spectra[idx]))
        #                           ], dtype=float)

    def get_xvec(self):
        # To keep the accuracy, kde is executed on the channel numbers in the
        # histogram data.
        # The actual energies were retrieved by "_real" variables.
        print("CHK", self.selected_spectra.shape[0])
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
        #np.random.seed(314)
        self.xvecorg = np.array(self.xvec)
        self.shift = np.random.uniform(-0.5, 0.5, size=self.xvec.shape[0])
        self.xvec += self.shift
        self.xvec_real += self.shift*self.de
        #print(self.xvec[0:30])

    # Since get_qlist_nova_class.get_alldata() takes the bin bottom energy as the neutron energy,
    # the random shift enegy in the following method is modified to be [0, 1]*de.
    # However, the  q values are still inaccurate.
    def add_shift_de(self):
        #np.random.seed(314)
        self.xvecorg = np.array(self.xvec)
        self.shift = np.random.uniform(0., 1.0, size=self.xvec.shape[0])
        self.xvec += self.shift
        self.xvec_real += self.shift*self.de
        #print(self.xvec[0:30])

    def run_ssvkernel(self):
        if self.optsm:
            tinmax = 10.**int(np.log10(self.selected_spectra.shape[0])+1.)
            self.tin = np.linspace(0.0, tinmax, int(tinmax)*2+1)
            print("Check parameters of horizontal axis")
            print("de=", self.de, "selected_energy[0]=",
                  self.selected_energy[0], "num channels=", self.tin.shape[0])
            self.tin_real = self.tin*self.de + self.selected_energy[0]
        else:
            self.tin = np.arange(self.selected_energy.shape[0])
            print("Check parameters of horizontal axis")
            print("de=", self.de, "selected_energy[0]=",
                  self.selected_energy[0], "num channels=", self.tin.shape[0])
            self.tin_real = np.linspace(self.selected_energy[0],
                                        self.selected_energy[-1],
                                        #num=self.selected_spectra.shape[0])
                                        num=800)
                                        #num=8000)
                                        #num=66670)
                                        #num=200000)
            print("num of tin_real element=",self.tin_real.shape)
            print(self.tin_real[0:10])
            print(self.selected_energy[0:20])
            print(self.selected_spectra[0:20])
        print('Total count number:', self.xvec_real.shape)
        self.y = ssvkernel.ssvkernel(self.xvec_real, self.tin_real, M=self.M,
                                     winparam=self.winparam,
                                     WinFunc=self.WinFunc)
        self.y_ = sskernel.sskernel(self.xvec_real, self.tin_real)

    def run_ssvkernel_direct(self):
        self.y = svk.ssvkernel(self.selected_spectra, self.selected_energy+self.de*0.5, M=self.M,
                               winparam=self.winparam, WinFunc=self.WinFunc)
        self.y_ = sk.sskernel(self.selected_spectra, self.selected_energy+self.de*0.5)
        self.tin_real = self.selected_energy+self.de*0.5

    def plotter(self):
        #scf = (np.min(self.xvec_real) - np.max(self.xvec_real)) /\
        #      (np.min(self.xvec) - np.max(self.xvec))

        if self.optsm:
            self.optsmear()
            snorms = self.sspectra/np.sum(self.sspectra)/self.sde

        norms = self.selected_spectra/np.sum(self.selected_spectra)/self.de
        fig = plt.figure(figsize=(6, 8))
        ax = fig.add_subplot(3, 1, 1)
        ax.bar(self.selected_energy, norms, width=self.de, label='expt data')
        if self.optsm:
            #ax.bar(self.senergy, snorms, width=self.sde, label='expt sdata')
            ax.plot(self.tin_real, self.yck/self.de, c='r', label='yck')
            ax.plot(self.tin_real, self.y[0]/self.de, c='k', label='ssvkernel')
        #else:
            #ax.plot(self.tin_real, self.y[0]/self.de, c='r', label='ssvkernel')
            #ax.plot(self.tin_real, self.y_[0]/self.de, c='k', label='sskernel')
            #ax.plot(self.tin_real, self.y[0], c='r', label='ssvkernel')
            #ax.plot(self.tin_real, self.y_[0], c='k', label='sskernel')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False, width=1.5)
        ax.set_ylabel('density')
        ax.set_xlabel('energy ($\mu eV$)')
        #ax.set_yticks([0, 35, 70])
        #ax.set_ylim([0, 70])
        ax.set_ylim(0., np.max(self.y_[0])*1.2)
        plt.legend()
        ax = fig.add_subplot(3, 1, 2)
        ax.set_ylabel('density')
        if self.optsm:
            ax.bar(self.senergy, snorms, width=self.sde, label='expt sdata',
                   bottom=0.0001)
            ax.plot(self.tin_real, self.yck/self.de, c='r', label='yck')
            ax.plot(self.tin_real, self.y[0]/self.de, c='k', label='ssvkernel')
        else:
            ax.bar(self.selected_energy, norms, width=self.de,
                   label='expt data', bottom=0.0001)
            #ax.plot(self.tin_real, self.y_[0]/self.de, c='k', label='sskernel')
            #ax.plot(self.tin_real, self.y[0]/self.de, c='r', label='ssvkernel')
            ax.plot(self.tin_real, self.y_[0], c='k', label='sskernel')
            ax.plot(self.tin_real, self.y[0], c='r', label='ssvkernel')
        ax.set_yscale('log')
        ax.set_ylim(0.0001, np.max(self.y_[0])/self.de)
        ax.set_ylim(0.0001, np.max(self.y_[0]))
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=False, width=1.5)
        tmprange = ax.get_xlim()
        plt.legend()
        ax = fig.add_subplot(3, 1, 3)
        #ax.plot(self.tin_real, self.y[2]*scf, c='r', label='ssvkernel')
        #ax.plot(self.tin_real, np.zeros_like(self.tin_real)+self.y_[2]*scf,
        ax.plot(self.tin_real, self.y[2], c='r', label='ssvkernel')
        ax.plot(self.tin_real, np.zeros_like(self.tin_real)+self.y_[2],
                c='k', label='sskernel')
        ax.set_ylabel('band-width')
        ax.set_xlabel('energy (meV)')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=True, width=1.5)
        ax.set_xlim(tmprange)
        #ax.set_ylim(0, np.max(self.y[2]*scf))
        #ax.set_ylim(0, 5)
        #ax.set_yticks([0, 0.003, 0.006])
        plt.subplots_adjust(hspace=0.0)
        plt.legend()
        plt.savefig(self.figname)
        if self.showplot:
            plt.show()

    def optsmear(self):
        #yinp = np.interp(self.tin_real, self.selected_energy,
        #                 self.selected_spectra)
        yinp = np.interp(self.tin_real, self.senergy, self.sspectra)
        dt = min(np.diff(self.tin))
        #thist = np.concatenate((self.tin, (self.tin[-1]+dt)[np.newaxis]))
        #y_hist = np.histogram(self.xvecorg, thist-dt/2)[0] / dt
        #idx = y_hist.nonzero()
        #t_nz = self.tin[idx]
        #y_hist_nz = y_hist[idx]
        yck = np.zeros_like(self.tin)
        for k in range(yck.shape[0]):
            #yck[k] = np.sum(y_hist_nz *
            #                self.Gauss(self.tin[k] - t_nz, self.y[2][k]))
            yck[k] = np.sum(yinp *
                            self.Gauss(self.tin[k] - self.tin, self.y[2][k]))
        self.yck = yck / np.sum(yck*dt)
        print(np.sum(yck))
        print(np.sum(self.y[0]))
        self.ys = (self.yck, self.y[1], self.y[2])

    def Gauss(self, x, w):
        return 1.0/((2.0*np.pi)**0.5*w)*np.exp(-(x/w)**2/2.0)

    def save_output(self,  output_file):
        dataset = {}
        dataset['y_ssvk'] = self.y
        dataset['y_ssk'] = self.y_
        dataset['tin_real'] = self.tin_real
        dataset['xlim'] = np.array([np.min(self.xvec_real),
                                   np.max(self.xvec_real)])
        with open(output_file, 'wb') as f:
            pickle.dump(dataset, f, -1)

    def save_outputs(self,  output_file):
        dataset = {}
        dataset['ys_ssvk'] = self.ys
        dataset['tin_real'] = self.tin_real
        #dataset['xlim'] = np.array([np.min(self.xvec_real),
        #                           np.max(self.xvec_real)])
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
