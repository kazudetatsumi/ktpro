#!/usr/bin/env python
# This script read the raw data of DNA through specific libraries and generate 3D arrays of q, omega, counts.
# Each of the arrays have axes of eca(element container array), ec(element container), and energy.
# This script also generate a 1D energy spectrum by integrating the counts over a specific q region.
# The spectrum is saved as 1D arrays of energy and counts.
# Note that the Cmm and Manyo library should be loaded on the computer which serves a DNA4 environment.
#
### We just read pklfiles containing numpy arrays. If you use this script on dna, uncomment the following two lines.
import Cmm
import Manyo
###
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import griddata
import copy

m = 1.674927471*10**(-27)   # [kg]
h = 6.62607015*10**(-34)    # [J. s]
meVtoJ = 1.60218*10**(-22)  # [J/meV]

print((meVtoJ * 2*m/(h**2))*10**(-20))
meVtoangsm2 = (1./0.81787)*0.01  # [Angs-2/meV]
print(meVtoangsm2)


class get_qlist:
    def __init__(self, runno, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile
        self.runno = runno

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
        EC = Cmm.GetHistogramMon(
                runNo=self.runno, useEiConv=True, LambdaParam="6.321/4.15",
                t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                TimeParam="-1.0/-1.0", UseFastChopper=True, isHistogram=False)
        Cmm.MutiplyConstant(dat=EC, factor=1e-09)
        self.DAT = Cmm.GetHistogramHW(
                runNo=self.runno, HwParam="0.000025/-0.05/0.15",
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                useT0ModCorr=False, TimeParam="-1.0/-1.0", UseFastChopper=True,
                tofOffsetFile="none", isHistogram=False)

    def get_sdata(self):
        EC = Cmm.GetHistogramMon(
                runNo=self.runno, useEiConv=True, LambdaParam="6.321/4.15",
                t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                TimeParam="-1.0/-1.0", UseFastChopper=True, isHistogram=False)
        Cmm.MutiplyConstant(dat=EC, factor=1e-09)
        DAT = Cmm.GetHistogramHW(
                runNo=self.runno, HwParam="0.000025/-0.05/0.15",
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                useT0ModCorr=False, TimeParam="-1.0/-1.0", UseFastChopper=True,
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
        self.get_data()
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
            #self.dataset = pickle.load(f, encoding='latin1')
            self.dataset = pickle.load(f)

    def create_fig(self):
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

    def spectra(self, qmin, qmax):
        x = np.ravel(self.dataset['q'])
        y = np.ravel(self.dataset['omega'])
        z = np.ravel(self.dataset['intensity'])
        print("CHCK", np.sort(np.unique(y)[0]))
        print("CHCK", np.sort(np.unique(y)[-1]))
        _x = x
        area = np.where((_x <= qmax) & (_x >= qmin))
        x = x[area[0]]
        y = y[area[0]]
        z = z[area[0]]
        ene = np.linspace(np.min(y), np.max(y), 801)
        z1d = np.zeros((801))
        for eidx in range(800):
            cond = np.where((y < ene[eidx+1]) & (y >= ene[eidx]))[0]
            z1d[eidx] = np.sum(z[cond])
        self.ene = ene
        self.spectra = z1d
        self.create_fig()
        plt.plot(ene, z1d)
        plt.yscale('log')
        plt.show()

    def save_spectra(self, spectrafile):
        dataset = {}
        dataset['spectra'] = self.spectra
        dataset['energy'] = self.ene
        with open(spectrafile, 'wb') as f:
            pickle.dump(dataset, f, -1)


def samplerun_after_dna():
    save_file = "./srlz/run6206.srlz"
    pklfile = "./srlz/run6206.pkl"
    spectrafile = "./srlz/run6206spectra.pkl"
    qmin = 0.55
    qmax = 0.70
    prj = get_qlist(save_file, pklfile=pklfile)
    prj.read_pkl()
    prj.spectra(qmin, qmax)
    prj.save_spectra(spectrafile)


def samplerun_after_sdna():
    save_file = "./srlz/run6206s.srlz"
    pklfile = "./srlz/run6206s.pkl"
    spectrafile = "./srlz/run6206sspectra.pkl"
    qmin = 0.55
    qmax = 0.70
    prj = get_qlist(save_file, pklfile=pklfile)
    prj.read_pkl()
    prj.spectra(qmin, qmax)
    prj.save_spectra(spectrafile)


def samplerun_on_dna():
    runnos = [6202, 6203, 6204, 6205, 6206, 6207]
    for runno in runnos:
        save_file = "./srlz/0000025/run" + str(runno) + ".srlz"
        pklfile = "./srlz/0000025/run" + str(runno) + ".pkl"
        prj = get_qlist(runno, save_file, pklfile=pklfile)
        prj.init_ecm()
        #prj.loadDAT()
        prj.get_all_data()
        prj.save_pkl()


def samplerun_on_dna_sdata():
    save_file = "./srlz/0000025/run6206s.srlz"
    pklfile = "./srlz/0000025/run6206s.pkl"
    prj = get_qlist(save_file=save_file, pklfile=pklfile)
    prj.init_eca()
    #prj.loadsDAT()
    prj.get_all_sdata()
    prj.save_pkl()


def check():
    save_file = "./srlz/run6206_half.srlz"
    pklfile = "./srlz/run6206_half.pkl"
    spectrafile = "./srlz/run6206spectra.pkl"
    qmin = 0.55
    qmax = 0.70
    prj = get_qlist(save_file, pklfile=pklfile)
    prj.read_pkl()
    prj.spectra(qmin, qmax)


samplerun_on_dna()
#samplerun_on_dna_sdata()

#check()

#samplerun_after_dna()
#samplerun_after_sdna()
