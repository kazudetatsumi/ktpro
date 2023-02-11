#!/usr/bin/env python
# This script read the raw data of DNA through specific libraries and generate 3D arrays of q, omega, counts.
# Each of the arrays have axes of eca(element container array), ec(element container), and energy.
# This script also generate a 1D energy spectrum by integrating the counts over a specific q region.
# The spectrum is saved as 1D arrays of energy and counts.
# Note that the Cmm and Manyo library should be loaded on the computer which serves a DNA4 environment.
#
### We just read pklfiles containing numpy arrays. If you use this script on dna, uncomment the following two lines.
#import Cmm
#import Manyo
###
import numpy as np
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import griddata
import copy
import datetime
from get_qlist_nova_class import get_qlist as gq

m = 1.674927471*10**(-27)   # [kg]
h = 6.62607015*10**(-34)    # [J. s]
meVtoJ = 1.60218*10**(-22)  # [J/meV]

print((meVtoJ * 2*m/(h**2))*10**(-20))
meVtoangsm2 = (1./0.81787)*0.01  # [Angs-2/meV]
print(meVtoangsm2)


class Sget_qlist(gq):
    def __init__(self, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile

    def get_org_data(self, binw):
        self.DAT = Cmm.GetHistogramHW(
                runNo=6204, HwParam=binw+"/-0.05/0.15",
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                useT0ModCorr=False, TimeParam="-1.0/-1.0", UseFastChopper=True,
                tofOffsetFile="none", isHistogram=False)
        self.EC = Cmm.GetHistogramMon(
                    runNo=6204, useEiConv=True, LambdaParam="6.321/4.15",
                    t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                    TimeParam="-1.0/-1.0", UseFastChopper=True, isHistogram=False)
        Cmm.MutiplyConstant(dat=self.EC, factor=1e-09)

    def get_org_intensity_array(self):
        self.intensity = np.zeros((self.DAT.PutSize(),
                                 self.DAT(0).PutSize(),
                                 len(self.DAT(0)(0).PutYList())))
        for ecaidx in range(0, self.DAT.PutSize()):
            for ecidx in range(0, self.DAT(0).PutSize()):
                self.intensity[ecaidx, ecidx, :] = np.array(self.DAT(ecaidx, ecidx)
                                                           .PutYList())

    def get_qemapb(self, intensityb):
        self.DATB = Manyo.ElementContainerMatrix(self.DAT)
        ctp = Manyo.CppToPython()
        for ecaidx in range(0, self.DATB.PutSize()):
            for ecidx in range(0, self.DATB(0).PutSize()):
                vecy = ctp.ListToDoubleVector(intensityb[ecaidx, ecidx, :].tolist())
                vece = ctp.ListToDoubleVector((intensityb[ecaidx, ecidx, :]**0.5).tolist())
                self.DATB(ecaidx, ecidx).Replace("Intensity", vecy)
                self.DATB(ecaidx, ecidx).Replace("Error", vece)
                self.DATB(ecaidx, ecidx).SetKeys("EnergyTransfer", "Intensity",
                                                 "Error")
        Cmm.DoMask(dat=self.DATB, filename="maskTY.txt")
        ECM = Cmm.ILambdaCorrDNA(dat=self.DATB, ec=self.EC, useMonEff=True)
        ECM2 = Cmm.SolidAngleCorrDNA(
                dat=ECM, useDetEff=True, useAbsoCorr=False, useEffCorr=False,
                sampletype="sample", sampleDataPath="test_sample_data.dat",
                DetEffDataPath="none")
        Cmm.MutiplyConstant(dat=ECM2, factor=1e-06)
        self.DATBQE = Cmm.CreateQEMap(dat=ECM2, startQ=0.0, endQ=2.0,
                                      deltaQ=0.05)

    def get_all_sdatab(self):
        q = np.zeros((self.DATBQE.PutSize(), len(self.DATBQE(0).PutYList())))
        omega = np.zeros((self.DATBQE.PutSize(), len(self.DATBQE(0).PutYList())))
        intensity = np.zeros((self.DATBQE.PutSize(), len(self.DATBQE(0).PutYList())))
        ones = np.ones((len(self.DATBQE(0).PutYList())))
        for ecidx in range(0, self.DATBQE.PutSize()):
            omega[ecidx, :] = np.array(self.DATBQE(ecidx).PutXList()[:-1])
            q[ecidx, :] = ones*self.DATBQE(ecidx).PutHeader().\
                PutDoubleVector('XRANGE')[0]
            intensity[ecidx, :] = np.array(self.DATBQE(ecidx).PutYList())
        self.dataset = {}
        self.dataset['omega'] = omega
        self.dataset['q'] = q
        self.dataset['intensity'] = intensity

    def save_pkl(self):
        with open(self.pklfile, 'wb') as f:
            pickle.dump(self.spectrab, f, -1)

    def load_pkl(self):
        with open(self.pklfile, 'rb') as f:
            self.spectrab = pickle.load(f, encoding='latin1')

    def get_boot_strap_sampled_spectra(self, nbs, seed=314):
        intensity1d = self.intensity.flatten().astype(int)
        nonzeroidx = np.nonzero(intensity1d)[0]
        x = np.array([idx for idx in nonzeroidx for num_repeat in
                     range(intensity1d[idx])], dtype=int)
        N = x.shape[0]
        np.random.seed(seed)
        intensityb = np.zeros_like(self.intensity)
        for inb in range(nbs):
            intensityb *= 0.
            print(inb)
            Nb = np.random.poisson(lam=N)
            idx = np.random.randint(0, N, Nb)
            xb = x[idx]
            test_idxs = np.unravel_index(xb, self.intensity.shape)
            for _idx0, _idx1, _idx2 in zip(test_idxs[0], test_idxs[1],
                                           test_idxs[2]):
                intensityb[_idx0, _idx1, _idx2] += 1
            print(datetime.datetime.now(), '1 resampling ended')
            self.get_qemapb(intensityb)
            print(datetime.datetime.now(), '1 resampled manyo-data obtained')
            self.get_all_sdatab()
            print(datetime.datetime.now(), '1 corrected resampled manyo-data obtained')
            self.spect(0.55, 0.70)
            if inb == 0:
                print(self.ene.shape)
                print(self.spectra.shape)
                self.spectrab = np.zeros((nbs, 2, self.ene.shape[0]))
            self.spectrab[inb, 0, :] = self.ene
            self.spectrab[inb, 1, :] = self.spectra


def run():
    prj = Sget_qlist(pklfile="run6204spectrab.pkl")
    prj.get_org_data("0.000025")
    #prj.get_org_data("0.001025")
    prj.get_org_intensity_array()
    print(datetime.datetime.now(), 'org data obtained')
    nbs = 1
    prj.get_boot_strap_sampled_spectra(nbs)
    prj.save_pkl()


def check():
    prj = Sget_qlist(pklfile="run6204spectrab.pkl")
    prj.load_pkl()
    for inb in range(prj.spectrab.shape[0]):
        plt.plot(prj.spectrab[inb, 0, :], prj.spectrab[inb, 1, :])
    plt.show()


#run()
#check()


