#!/usr/bin/env python
# This script read the raw neutron count data of DNA without any corrections and resample the count data
# and apply the necessary corrections to deduce resampled QENS data corresponding to dobule differential 
# cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
try:
    import Cmm
    import Manyo
except ModuleNotFoundError as err:
    print(err)
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.interpolate import griddata
import datetime
from get_qlist_nova_class import get_qlist as gq

m = 1.674927471*10**(-27)   # [kg]
h = 6.62607015*10**(-34)    # [J. s]
meVtoJ = 1.60218*10**(-22)  # [J/meV]

#print((meVtoJ * 2*m/(h**2))*10**(-20))
meVtoangsm2 = (1./0.81787)*0.01  # [Angs-2/meV]
#print(meVtoangsm2)


class Sget_qlist(gq):
    def __init__(self, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile

    def get_org_data(self, binw, runNo):
        self.DAT = Cmm.GetHistogramHW(
                runNo=runNo, HwParam=binw+"/-0.05/0.15",
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                useT0ModCorr=False, TimeParam="-1.0/-1.0", UseFastChopper=True,
                tofOffsetFile="none", isHistogram=False)
        self.EC = Cmm.GetHistogramMon(
                    runNo=runNo, useEiConv=True, LambdaParam="6.321/4.15",
                    t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                    TimeParam="-1.0/-1.0", UseFastChopper=True,
                    isHistogram=False)
        Cmm.MutiplyConstant(dat=self.EC, factor=1e-09)

    def get_org_intensity_array(self):
        self.intensity = np.zeros((self.DAT.PutSize(),
                                   self.DAT(0).PutSize(),
                                   len(self.DAT(0)(0).PutYList())))
        for ecaidx in range(0, self.DAT.PutSize()):
            for ecidx in range(0, self.DAT(0).PutSize()):
                self.intensity[ecaidx, ecidx, :] = np.array(
                                                   self.DAT(ecaidx, ecidx)
                                                   .PutYList())

    def get_qemapb(self, intensityb):
        self.DATB = Manyo.ElementContainerMatrix(self.DAT)
        ctp = Manyo.CppToPython()
        for ecaidx in range(0, self.DATB.PutSize()):
            for ecidx in range(0, self.DATB(0).PutSize()):
                vecy = ctp.ListToDoubleVector(intensityb[ecaidx,
                                                         ecidx, :].tolist())
                vece = ctp.ListToDoubleVector((intensityb[ecaidx,
                                               ecidx, :]**0.5).tolist())
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

    def get_qemap(self):
        Cmm.DoMask(dat=self.DAT, filename="maskTY.txt")
        ECM = Cmm.ILambdaCorrDNA(dat=self.DAT, ec=self.EC, useMonEff=True)
        ECM2 = Cmm.SolidAngleCorrDNA(
                dat=ECM, useDetEff=True, useAbsoCorr=False, useEffCorr=False,
                sampletype="sample", sampleDataPath="test_sample_data.dat",
                DetEffDataPath="none")
        Cmm.MutiplyConstant(dat=ECM2, factor=1e-06)
        self.DATQE = Cmm.CreateQEMap(dat=ECM2, startQ=0.0, endQ=2.0,
                                     deltaQ=0.05)

    def get_all_sdatab(self):
        q = np.zeros((self.DATBQE.PutSize(), len(self.DATBQE(0).PutYList())))
        omega = np.zeros((self.DATBQE.PutSize(), len(self.DATBQE(0).PutYList())
                          ))
        intensity = np.zeros((self.DATBQE.PutSize(),
                              len(self.DATBQE(0).PutYList())))
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

    def get_all_sdata(self):
        q = np.zeros((self.DATQE.PutSize(), len(self.DATQE(0).PutYList())))
        omega = np.zeros((self.DATQE.PutSize(), len(self.DATQE(0).PutYList())))
        intensity = np.zeros((self.DATQE.PutSize(),
                              len(self.DATQE(0).PutYList())))
        ones = np.ones((len(self.DATQE(0).PutYList())))
        for ecidx in range(0, self.DATQE.PutSize()):
            omega[ecidx, :] = np.array(self.DATQE(ecidx).PutXList()[:-1])
            q[ecidx, :] = ones*self.DATQE(ecidx).PutHeader().\
                PutDoubleVector('XRANGE')[0]
            intensity[ecidx, :] = np.array(self.DATQE(ecidx).PutYList())
        self.dataset = {}
        self.dataset['omega'] = omega
        self.dataset['q'] = q
        self.dataset['intensity'] = intensity

    def save_pkl(self):
        with open(self.pklfile, 'wb') as f:
            pickle.dump(self.spectrab, f, -1)

    def load_pkl(self, ispython2=False):
        with open(self.pklfile, 'rb') as f:
            if ispython2:
                self.spectrab = pickle.load(f, encoding='latin1')
            else:
                self.spectrab = pickle.load(f)

    def get_boot_strap_sampled_spectra(self, nbs, qmin, qmax, seed=314):
        intensity1d = self.intensity.flatten().astype(int)
        nonzeroidx = np.nonzero(intensity1d)[0]
        x = np.array([idx for idx in nonzeroidx for num_repeat in
                     range(intensity1d[idx])], dtype=int)
        N = x.shape[0]
        np.random.seed(seed)
        intensityb = np.zeros_like(self.intensity)
        randomstates = []
        for inb in range(nbs):
            randomstates.append(np.random.get_state())
            Nb = np.random.poisson(lam=N)
            idx = np.random.randint(0, N, Nb)
        with open('randomstates.pkl', 'wb') as f:
            pickle.dump(randomstates, f, -1)
        for inb in range(nbs):
            intensityb *= 0.
            print(inb)
            #randomstates.append(np.random.get_state())
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
            print(datetime.datetime.now(),
                  '1 corrected resampled manyo-data obtained')
            self.spect(qmin, qmax, isplot=True)
            if inb == 0:
                print(self.ene.shape)
                print(self.spectra.shape)
                self.spectrab = np.zeros((nbs, 2, self.ene.shape[0]))
            self.spectrab[inb, 0, :] = self.ene
            self.spectrab[inb, 1, :] = self.spectra
        #with open('randomstates2.pkl', 'wb') as f:
        #    pickle.dump(randomstates, f, -1)

        #with open('tst.pkl', 'wb') as f:
        #    pickle.dump(intensityb, f, -1)
        #    print('intensityb has been saved in tst.pkl')
        ### Add the original data at the end of the 1st index
        ###self.get_qemapb(self.intensity)
        ###self.get_all_sdatab()
        ###self.spect(0.55, 0.70)
        ###self.spectrab[-1, 0, :] = self.ene
        ###self.spectrab[-1, 1, :] = self.spectra


def run():
    prj = Sget_qlist(pklfile="run6202spectrab.pkl")
    #prj.get_org_data("0.000025")
    prj.get_org_data("0.00025")
    print(datetime.datetime.now(), 'org_data ended')
    prj.get_org_intensity_array()
    print(datetime.datetime.now(), 'org_intensity_array ended')
    nbs = 1
    qmin = 0.55
    qmax = 0.70
    prj.get_boot_strap_sampled_spectra(nbs, qmin, qmax)
    print(datetime.datetime.now(), 'boot_strap_sampled_spectra ended')
    prj.save_pkl()
    intensities = np.sort(np.unique(prj.intensity.flatten()))
    print(intensities[0:5])
    print(intensities[-5:-1])


def run_org():
    prj = Sget_qlist(pklfile="run6202spectraorg.pkl")
    prj.get_org_data("0.000025")
    prj.get_qemap()
    prj.get_all_sdata()
    prj.spect(0.55, 0.70, isplot=True)
    prj.spectrab = np.zeros((1, 2, prj.ene.shape[0]))
    prj.spectrab[0, 0, :] = prj.ene
    prj.spectrab[0, 1, :] = prj.spectra
    prj.save_pkl()


def check():
    prj = Sget_qlist(pklfile="run6202spectraorg.pkl")
    prj.load_pkl()
    for inb in range(prj.spectrab.shape[0]):
        plt.plot(prj.spectrab[inb, 0, :], prj.spectrab[inb, 1, :], lw=0.2)
    prj = Sget_qlist(pklfile="run6202spectrab.pkl")
    prj.load_pkl()
    for inb in range(prj.spectrab.shape[0]):
        plt.plot(prj.spectrab[inb, 0, :], prj.spectrab[inb, 1, :], lw=0.2)
    #plt.xlim(-0.025, 0.075)
    plt.ylim(.0, 8.)
    plt.show()


#run()
#check()
#run_org()