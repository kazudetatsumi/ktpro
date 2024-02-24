#!/usr/bin/env python
# This script read the raw data of DNA through specific libraries and generate 3D arrays of q, omega, counts.
# Each of the arrays have axes of eca(element container array), ec(element container), and energy.
# This script also generate a 1D energy spectrum by integrating the counts over a specific q region.
# The spectrum is saved as 1D arrays of energy and counts.
# Note that the Cmm and Manyo library should be loaded on the computer which serves a DNA4 environment.
#
### We just read pklfiles containing numpy arrays. If you use this script on dna, uncomment the following two lines.
try:
    import Cmm
except ModuleNotFoundError:
    pass
try:
    import Manyo
except ModuleNotFoundError:
    pass
#import Cmm
#import Manyo
###
import numpy as np
import os
import pickle
from scipy.interpolate import griddata
import copy

m = 1.674927471*10**(-27)   # [kg]
h = 6.62607015*10**(-34)    # [J. s]
meVtoJ = 1.60218*10**(-22)  # [J/meV]

#print((meVtoJ * 2*m/(h**2))*10**(-20))
meVtoangsm2 = (1./0.81787)*0.01  # [Angs-2/meV]
#print(meVtoangsm2)


class get_qlist:
    def __init__(self, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile

    def loadDAT(self):
            self.DAT = Manyo.ElementContainerMatrix()
            r = Manyo.ReadSerializationFileBinary(self.save_file)
            r.Load(self.DAT)
            print("DAT loaded")

    def loadsDAT(self):
            self.DAT = Manyo.ElementContainerArray()
            r = Manyo.ReadSerializationFileBinary(self.save_file)
            r.Load(self.DAT)
            print("sDAT loaded")

    def saveDAT(self):
        w = Manyo.WriteSerializationFileBinary(self.save_file)
        w.Save(self.DAT)

    def get_data(self):
        self.HwParam = "0.00025/-0.05/0.15"
        EC = Cmm.GetHistogramMon(
                runNo=6204, useEiConv=True, LambdaParam="6.321/4.15",
                t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                #TimeParam="-1.0/-1.0", UseFastChopper=True, isHistogram=False)
                TimeParam="0.0,5843.0", UseFastChopper=True, isHistogram=False)
        Cmm.MutiplyConstant(dat=EC, factor=1e-09)
        self.DAT = Cmm.GetHistogramHW(
                runNo=6204, HwParam="0.00025/-0.05/0.15",
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                #useT0ModCorr=False, TimeParam="-1.0/-1.0", UseFastChopper=True,
                useT0ModCorr=False, TimeParam="0.0, 5843.0", UseFastChopper=True,
                tofOffsetFile="none", isHistogram=False)

    def get_sdata(self):
        self.HwParam = "0.00025/-0.05/0.15"
        EC = Cmm.GetHistogramMon(
                runNo=6204, useEiConv=True, LambdaParam="6.321/4.15",
                t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                #TimeParam="-1.0/-1.0", UseFastChopper=True, isHistogram=False)
                TimeParam="0.0, 5843.0", UseFastChopper=True, isHistogram=False)
        Cmm.MutiplyConstant(dat=EC, factor=1e-09)
        DAT = Cmm.GetHistogramHW(
                runNo=6204, HwParam="0.00025/-0.05/0.15",
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                #useT0ModCorr=False, TimeParam="-1.0/-1.0", UseFastChopper=True,
                useT0ModCorr=False, TimeParam="0.0, 5843.0", UseFastChopper=True,
                tofOffsetFile="none", isHistogram=False)
        Cmm.DoMask(dat=DAT, filename="maskTY.txt")
        ECM = Cmm.ILambdaCorrDNA(dat=DAT, ec=EC, useMonEff=True)
        ECM2 = Cmm.SolidAngleCorrDNA(
                dat=ECM, useDetEff=True, useAbsoCorr=False, useEffCorr=False,
                sampletype="sample", sampleDataPath="test_sample_data.dat",
                DetEffDataPath="none")
        Cmm.MutiplyConstant(dat=ECM2, factor=1e-06)
        self.DAT = Cmm.CreateQEMap(dat=ECM2, startQ=0.0, endQ=2.0, deltaQ=0.05)

    def init_ecm(self):
        print("get_data is starting")
        self.get_data()
        print("get_data is finished")
        self.saveDAT()
        print("Done!")

    def init_eca(self):
        self.get_sdata()
        self.saveDAT()
        print("sDone!")

    def get_all_data(self):
        Ef = np.zeros((self.DAT.PutSize(),
                       self.DAT(0).PutSize(),
                       len(self.DAT(0)(0).PutYList())))
        theta = np.zeros((self.DAT.PutSize(),
                          self.DAT(0).PutSize(),
                          len(self.DAT(0)(0).PutYList())))
        omega = np.zeros((self.DAT.PutSize(),
                          self.DAT(0).PutSize(),
                          len(self.DAT(0)(0).PutYList())))
        intensity = np.zeros((self.DAT.PutSize(),
                              self.DAT(0).PutSize(),
                              len(self.DAT(0)(0).PutYList())))
        ones = np.ones((len(self.DAT(0)(0).PutYList())))
        for ecaidx in range(0, self.DAT.PutSize()):
            for ecidx in range(0, self.DAT(0).PutSize()):
                Ef[ecaidx, ecidx, :] = ones*self.DAT(ecaidx, ecidx)\
                                                   .PutHeader()\
                                                   .PutDouble('Ef')
                theta[ecaidx, ecidx, :] = ones*self.DAT(ecaidx, ecidx)\
                                                   .PutHeader()\
                                                   .PutDouble('PolarAngle')
                # We get omega as the same manner as in the Matsuura's script
                omega[ecaidx, ecidx, :] = np.array(self.DAT(ecaidx, ecidx)
                                                   .PutXList()[:-1])
                intensity[ecaidx, ecidx, :] = np.array(self.DAT(ecaidx, ecidx)
                                                       .PutYList())
        q = ((2.*Ef + omega - 2.*np.cos(theta)
             * (Ef*(Ef+omega))**0.5)*meVtoangsm2)**0.5*2*np.pi
        self.dataset = {}
        self.dataset['omega'] = omega
        self.dataset['q'] = q
        self.dataset['intensity'] = intensity

    def get_intdtype(self, maxnumber):
        if maxnumber <= 255 - 15:
            intdtype = 'uint8'
        elif maxnumber <= 65535 - 510:
            intdtype = 'uint16'
        elif maxnumber <= 4294967295 - 131070:
            intdtype = 'uint32'
        else:
            intdtype = 'uint64'
        return intdtype

    def get_all_data2(self):
        Ef = np.zeros((self.DAT.PutSize(),
                       self.DAT(0).PutSize(),
                       len(self.DAT(0)(0).PutYList())))
        theta = np.zeros((self.DAT.PutSize(),
                          self.DAT(0).PutSize(),
                          len(self.DAT(0)(0).PutYList())))
        omega = np.zeros((self.DAT.PutSize(),
                          self.DAT(0).PutSize(),
                          len(self.DAT(0)(0).PutYList())))
        intensity = np.zeros((self.DAT.PutSize(),
                              self.DAT(0).PutSize(),
                              len(self.DAT(0)(0).PutYList())))
        ones = np.ones((len(self.DAT(0)(0).PutYList())))
        for ecaidx in range(0, self.DAT.PutSize()):
            for ecidx in range(0, self.DAT(0).PutSize()):
                Ef[ecaidx, ecidx, :] = ones*self.DAT(ecaidx, ecidx)\
                                                   .PutHeader()\
                                                   .PutDouble('Ef')
                theta[ecaidx, ecidx, :] = ones*self.DAT(ecaidx, ecidx)\
                                                   .PutHeader()\
                                                   .PutDouble('PolarAngle')
                # We get omega as the same manner as in the Matsuura's script
                omega[ecaidx, ecidx, :] = np.array(self.DAT(ecaidx, ecidx)
                                                   .PutXList()[:-1])
                intensity[ecaidx, ecidx, :] = np.array(self.DAT(ecaidx, ecidx)
                                                       .PutYList())
        q = ((2.*Ef + omega - 2.*np.cos(theta)
             * (Ef*(Ef+omega))**0.5)*meVtoangsm2)**0.5*2*np.pi
        self.dataset = {}
        self.dataset['omega'] = omega
        self.dataset['q'] = q
        if np.max(intensity) <= 255 - 15:
            intdtype = 'uint8'
        elif np.max(intensity) <= 65535 - 510:
            intdtype = 'uint16'
        elif np.max(intensity) <= 4294967295 - 131070:
            intdtype = 'uint32'
        else:
            intdtype = 'uint64'
        self.dataset['intensity'] = np.array(intensity, dtype=intdtype)

    def get_all_data3(self):
        Ef = np.zeros((self.DAT.PutSize(),
                       self.DAT(0).PutSize(),
                       len(self.DAT(0)(0).PutYList())))
        theta = np.zeros((self.DAT.PutSize(),
                          self.DAT(0).PutSize(),
                          len(self.DAT(0)(0).PutYList())))
        omega = np.zeros((self.DAT.PutSize(),
                          self.DAT(0).PutSize(),
                          len(self.DAT(0)(0).PutYList())))
        intensity = np.zeros((self.DAT.PutSize(),
                              self.DAT(0).PutSize(),
                              len(self.DAT(0)(0).PutYList())))
        ones = np.ones((len(self.DAT(0)(0).PutYList())))
        for ecaidx in range(0, self.DAT.PutSize()):
            for ecidx in range(0, self.DAT(0).PutSize()):
                Ef[ecaidx, ecidx, :] = ones*self.DAT(ecaidx, ecidx)\
                                                   .PutHeader()\
                                                   .PutDouble('Ef')
                theta[ecaidx, ecidx, :] = ones*self.DAT(ecaidx, ecidx)\
                                                   .PutHeader()\
                                                   .PutDouble('PolarAngle')
                # We get omega as the same manner as in the Matsuura's script
                omega[ecaidx, ecidx, :] = np.array(self.DAT(ecaidx, ecidx)
                                                   .PutXList()[:-1])
                intensity[ecaidx, ecidx, :] = np.array(self.DAT(ecaidx, ecidx)
                                                       .PutYList())
        q = ((2.*Ef + omega - 2.*np.cos(theta)
             * (Ef*(Ef+omega))**0.5)*meVtoangsm2)**0.5*2*np.pi
        self.dataset = {}
        self.dataset['omega'] = np.array(omega, dtype='float32')
        self.dataset['q'] = np.array(q, dtype='float32')
        if np.max(intensity) <= 255 - 15:
            intdtype = 'uint8'
        elif np.max(intensity) <= 65535 - 510:
            intdtype = 'uint16'
        elif np.max(intensity) <= 4294967295 - 131070:
            intdtype = 'uint32'
        else:
            intdtype = 'uint64'
        self.dataset['intensity'] = np.array(intensity, dtype=intdtype)

    def get_all_sdata(self):
        q = np.zeros((self.DAT.PutSize(), len(self.DAT(0).PutYList())))
        omega = np.zeros((self.DAT.PutSize(), len(self.DAT(0).PutYList())))
        intensity = np.zeros((self.DAT.PutSize(), len(self.DAT(0).PutYList())))
        ones = np.ones((len(self.DAT(0).PutYList())))
        for ecidx in range(0, self.DAT.PutSize()):
            omega[ecidx, :] = np.array(self.DAT(ecidx).PutXList()[:-1])
            q[ecidx, :] = ones*self.DAT(ecidx).PutHeader().\
                PutDoubleVector('XRANGE')[0]
            intensity[ecidx, :] = np.array(self.DAT(ecidx).PutYList())
        self.dataset = {}
        self.dataset['omega'] = omega
        self.dataset['q'] = q
        self.dataset['intensity'] = intensity

    def save_pkl(self):
        with open(self.pklfile, 'wb') as f:
            pickle.dump(self.dataset, f, -1)

    def read_pkl(self):
        with open(self.pklfile, 'rb') as f:
            self.dataset = pickle.load(f, encoding='latin1')
            #self.dataset = pickle.load(f)
            
    def save_hdf5(self):
        with open(self.hdf5file, 'wb') as hf:
            hf.create_dataset('omega', data=self.dataset['omega'])
            hf.create_dataset('q', data=self.dataset['q'])
            hf.create_dataset('intensity', data=self.dataset['intensity'])

    def create_fig(self):
        import matplotlib.pyplot as plt
        self.fig = plt.figure(figsize=(10, 8))

    def getXYZ(self, x, y, z, xb, yb):
        ngridx = 100
        ngridy = 100
        xlin = np.linspace(xb[0], xb[1], ngridx)
        ylin = np.linspace(yb[0], yb[1], ngridy)
        X, Y = np.meshgrid(xlin, ylin)
        Z = griddata((x, y), z, (X, Y), method='linear', fill_value=-1)
        return X, Y, Z

    def plot_sub(self, X, Y, Z, axindx, vmax):
        import matplotlib.pyplot as plt
        ax = self.fig.add_subplot(1, 1, axindx)
        cmap1 = copy.copy(plt.cm.jet)
        cmap1.set_under('w')
        c = ax.pcolor(X, Y, Z, cmap=cmap1, vmin=0, vmax=vmax)
        ax.axis('tight')
        if axindx == 2:
            ax.set_ylabel('hw (meV)')
            #ax.set_yticks([1, 8])
            ax.set_xlabel('q (Angs-1)')
        self.fig.colorbar(c, ax=ax)

    def qemap(self, qmin, qmax):
        import matplotlib.pyplot as plt
        x = np.ravel(self.dataset['q'])
        y = np.ravel(self.dataset['omega'])
        z = np.ravel(self.dataset['intensity'])
        _x = x
        area = np.where((_x <= qmax) & (_x >= qmin))
        x = x[area[0]]
        y = y[area[0]]
        z = z[area[0]]
        self.create_fig()

        X, Y, Z = self.getXYZ(x, y, z, [min(x), max(x)], [min(y), max(y)])
        self.plot_sub(Y, X, Z, 1, np.max(Z))
        plt.show()

    def spect(self, qmin, qmax, dataset, isplot=False):
        x = np.ravel(dataset['q'])
        y = np.ravel(dataset['omega'])
        z = np.ravel(dataset['intensity'])
        #print("CHECK SPECT:, total intennsity:", np.sum(z[z < 1.0e+99]))
        _x = x
        area = np.where((_x <= qmax - 1.0e-14) & (_x > qmin - 1.0e-14))
        #area = np.where((_x <= qmax) & (_x > qmin))
        print(np.unique(x[area[0]]))
        x = x[area[0]]
        y = y[area[0]]
        z = z[area[0]]
        ene = np.sort(np.unique(y))
        z1d = np.zeros_like(ene)
        for eidx, e in enumerate(ene):
            #if eidx % (ene.shape[0]/50) == 0:
            #    print(int(eidx/float(ene.shape[0])*100), " %")
            #cond = np.where((y < e) & (y >= e - float(hw[0])))[0]
            cond = np.where(np.abs(y-e) < 0.0000000001)
            z1d[eidx] = np.sum(z[cond])
        self.ene = ene
        self.spectra = z1d
        #print("CHECK SPECT:, total intensity with q of", np.unique(x), ":", np.sum(self.spectra))
        if isplot:
            import matplotlib.pyplot as plt
            self.create_fig()
            plt.plot(ene, z1d)
            #plt.yscale('log')
            plt.show()

    def spect2(self, qmin, qmax, dataset, isplot=False):
        x = np.ravel(dataset['q'])
        y = np.ravel(dataset['omega'])
        z = np.ravel(dataset['intensity'])
        #print("CHECK SPECT:, total intennsity:", np.sum(z[z < 1.0e+99]))
        _x = x
        area = np.where((_x <= qmax - 1.0e-14) & (_x > qmin - 1.0e-14))
        #area = np.where((_x <= qmax) & (_x > qmin))
        print(np.unique(x[area[0]]))
        x = x[area[0]]
        y = y[area[0]]
        z = z[area[0]]
        #ene = np.sort(np.unique(y))
        #z1d = np.zeros_like(ene)
        _z = z[np.argsort(y)]
        _y = np.sort(y)
        ene, mult = np.unique(_y, return_counts=True)
        z1d = np.zeros_like(ene)
        for eidx, e in enumerate(ene):
            lb = np.sum(mult[:eidx])
            ub = lb + mult[eidx]
            z1d[eidx] = np.sum(_z[lb:ub])
        self.ene = ene
        self.spectra = z1d
        #print("CHECK SPECT:, total intensity with q of", np.unique(x), ":", np.sum(self.spectra))
        if isplot:
            import matplotlib.pyplot as plt
            self.create_fig()
            plt.plot(ene, z1d)
            #plt.yscale('log')
            plt.show()

    def spect3(self, qmin, qmax, dataset, isplot=False):
        x = np.ravel(dataset['q'])
        y = np.ravel(dataset['omega'])
        z = np.ravel(dataset['intensity'])
        #print("CHECK SPECT:, total intennsity:", np.sum(z[z < 1.0e+99]))
        _x = x
        area = np.where((_x <= qmax - 1.0e-7) & (_x > qmin - 1.0e-7))
        #area = np.where((_x <= qmax) & (_x > qmin))
        print(np.unique(x[area[0]]))
        x = x[area[0]]
        y = y[area[0]]
        z = z[area[0]]
        #ene = np.sort(np.unique(y))
        #z1d = np.zeros_like(ene)
        _z = z[np.argsort(y)]
        _y = np.sort(y)
        ene, mult = np.unique(_y, return_counts=True)
        z1d = np.zeros_like(ene)
        for eidx, e in enumerate(ene):
            lb = np.sum(mult[:eidx])
            ub = lb + mult[eidx]
            z1d[eidx] = np.sum(_z[lb:ub])
        self.ene = ene
        self.spectra = z1d
        #print("CHECK SPECT:, total intensity with q of", np.unique(x), ":", np.sum(self.spectra))
        if isplot:
            import matplotlib.pyplot as plt
            self.create_fig()
            plt.plot(ene, z1d)
            #plt.yscale('log')
            plt.show()

    def spect3e(self, qmin, qmax, dataset):
        x = np.ravel(dataset['q'])
        y = np.ravel(dataset['omega'])
        z = np.ravel(dataset['intensity'])
        err = np.ravel(dataset['error'])
        _x = x
        area = np.where((_x <= qmax - 1.0e-7) & (_x > qmin - 1.0e-7))
        x = x[area[0]]
        y = y[area[0]]
        z = z[area[0]]
        err = err[area[0]]
        _z = z[np.argsort(y)]
        _err = err[np.argsort(y)]
        _y = np.sort(y)
        ene, mult = np.unique(_y, return_counts=True)
        z1d = np.zeros_like(ene)
        err1d = np.zeros_like(ene)
        for eidx, e in enumerate(ene):
            lb = np.sum(mult[:eidx])
            ub = lb + mult[eidx]
            z1d[eidx] = np.sum(_z[lb:ub])
            err1d[eidx] = (np.sum(_err[lb:ub]**2))**0.5
        self.ene = ene
        self.spectra = z1d
        self.err = err1d

    def spectm(self, qmin, qmax, dataset):
        for iq, (qmi, qma) in enumerate(zip(qmin, qmax)):
            self.spect3(qmi, qma, dataset)
            if iq == 0:
                _spectra = np.zeros((len(qmin), self.spectra.shape[0]))
                _ene = np.zeros((len(qmin), self.ene.shape[0]))
            _spectra[iq, :] = self.spectra
            _ene[iq, :] = self.ene
        self.spectra = _spectra.flatten()
        self.ene = _ene.flatten()

    def spectme(self, qmin, qmax, dataset):
        for iq, (qmi, qma) in enumerate(zip(qmin, qmax)):
            self.spect3e(qmi, qma, dataset)
            if iq == 0:
                _spectra = np.zeros((len(qmin), self.spectra.shape[0]))
                _ene = np.zeros((len(qmin), self.ene.shape[0]))
                _err = np.zeros((len(qmin), self.err.shape[0]))
            _spectra[iq, :] = self.spectra
            _ene[iq, :] = self.ene
            _err[iq, :] = self.err
        self.spectra = _spectra.flatten()
        self.ene = _ene.flatten()
        self.err = _err.flatten()

    def save_spectra(self, spectrafile, old=False):
        dataset = {}
        dataset['spectra'] = self.spectra
        if not old:
            dataset['energy'] = self.ene
        with open(spectrafile, 'wb') as f:
            pickle.dump(dataset, f, -1)

    def get_all_monidata(self):
        self.spectra = np.zeros((1, 3, len(self.EC.PutYList())))
        self.spectra[0, 0, :] = self.EC.PutXList()[:-1]
        self.spectra[0, 1, :] = self.EC.PutYList()
        self.spectra[0, 2, :] = self.EC.PutEList()

    def run_moni(self):
        self.get_all_monidata()
        self.save_spectra(self.pklfile + ".moni", old=True)


def samplerun_after_dna():
    save_file = "./srlz/run6204_half.srlz"
    pklfile = "./srlz/run6204_half.pkl"
    spectrafile = "./srlz/run6204spectra_half.pkl"
    qmin = 0.55
    qmax = 0.70
    prj = get_qlist(save_file, pklfile=pklfile)
    prj.read_pkl()
    prj.spectra(qmin, qmax)
    prj.save_spectra(spectrafile)


def samplerun_after_sdna():
    save_file = "./srlz/run6204s_half.srlz"
    pklfile = "./srlz/run6204s_half.pkl"
    spectrafile = "./srlz/run6204sspectra_half.pkl"
    qmin = 0.55
    qmax = 0.70
    prj = get_qlist(save_file, pklfile=pklfile)
    prj.read_pkl()
    prj.spectra(qmin, qmax)
    prj.save_spectra(spectrafile)


def samplerun_on_dna():
    save_file = "./srlz/run6204_half.srlz"
    pklfile = "./srlz/run6204_half.pkl"
    prj = get_qlist(save_file, pklfile=pklfile)
    prj.init_ecm()
    #prj.loadDAT()
    prj.get_all_data()
    prj.save_pkl()


def samplerun_on_dna_sdata():
    save_file = "./srlz/run6204s_half.srlz"
    pklfile = "./srlz/run6204s_half.pkl"
    prj = get_qlist(save_file=save_file, pklfile=pklfile)
    prj.init_eca()
    #prj.loadsDAT()
    prj.get_all_sdata()
    prj.save_pkl()


#samplerun_on_dna()
#samplerun_on_dna_sdata()

#samplerun_after_dna()
#samplerun_after_sdna()

