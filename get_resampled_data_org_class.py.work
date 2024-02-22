#!/usr/bin/env python
# This script was folked from get_resampled_data_mpi_class.py, only having
# the part not related to the bootstrap sampling.
# This script reads a set of raw neutron count data of DNA without any
# corrections, i.e., an uncorrected energy container matrix using utsusemi DNA
# facede Cmm, as well as applies usual corrections on the raw data.
# The uncorrected and corrected  spectral data are stored as runXXXspectrab.pkl
# to utilize the analysis scripts for bootstrap sampled data sets.
# Kazuyoshi TATSUMI 2023/06/07

try:
    import Cmm
except ModuleNotFoundError as err:
    print(err)
import numpy as np
import pickle
from get_qlist_nova_class import get_qlist as gq


class Sget_qlist(gq):
    def __init__(self, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile

    def get_org_data(self, binw, runNo, TimeParam="-1.0/-1.0", frac=None):
        self.DAT = Cmm.GetHistogramHW(
                runNo=runNo, HwParam=binw+"/-0.05/0.15",
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                useT0ModCorr=False, TimeParam=TimeParam, UseFastChopper=True,
                tofOffsetFile="none", isHistogram=False)
        if frac:
            TimeParam = self.get_frac_TimeParam(TimeParam, frac)
        self.EC = Cmm.GetHistogramMon(
                    runNo=runNo, useEiConv=True, LambdaParam="6.321/4.15",
                    t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                    TimeParam=TimeParam, UseFastChopper=True,
                    isHistogram=False)
        Cmm.MutiplyConstant(dat=self.EC, factor=1e-09)

    def get_frac_TimeParam(self, TimeParam, frac):
        it = float(TimeParam.split(",")[0])
        ft = float(TimeParam.split(",")[1])
        return str(np.round(ft - (ft-it)*frac)) + "," + str(ft)

    def get_org_intensity_array(self):
        self.intensity = np.zeros((self.DAT.PutSize(),
                                   self.DAT(0).PutSize(),
                                   len(self.DAT(0)(0).PutYList())))
        for ecaidx in range(0, self.DAT.PutSize()):
            for ecidx in range(0, self.DAT(0).PutSize()):
                self.intensity[ecaidx, ecidx, :] = np.array(self.DAT(ecaidx,
                                                            ecidx).PutYList())

    def get_qemap(self, qmin, qmax):
        Cmm.DoMask(dat=self.DAT, filename="maskTY.txt")
        ECM = Cmm.ILambdaCorrDNA(dat=self.DAT, ec=self.EC, useMonEff=True)
        ECM2 = Cmm.SolidAngleCorrDNA(
                dat=ECM, useDetEff=True, useAbsoCorr=False, useEffCorr=False,
                sampletype="sample", sampleDataPath="test_sample_data.dat",
                DetEffDataPath="none")
        Cmm.MutiplyConstant(dat=ECM2, factor=1e-06)
        DATQE = Cmm.CreateQEMap(dat=ECM2, startQ=0.0, endQ=2.0, deltaQ=0.05)
        dataset = self.get_all_sdata(DATQE)
        self.spect(qmin, qmax, dataset)

    def get_all_sdata(self, DATQE):
        q = np.zeros((DATQE.PutSize(), len(DATQE(0).PutYList())))
        omega = np.zeros((DATQE.PutSize(), len(DATQE(0).PutYList())))
        intensity = np.zeros((DATQE.PutSize(), len(DATQE(0).PutYList())))
        ones = np.ones((len(DATQE(0).PutYList())))
        for ecidx in range(0, DATQE.PutSize()):
            omega[ecidx, :] = np.array(DATQE(ecidx).PutXList()[:-1])
            q[ecidx, :] = ones*DATQE(ecidx).PutHeader().PutDoubleVector(
                                                                  'XRANGE')[0]
            intensity[ecidx, :] = np.array(DATQE(ecidx).PutYList())
        dataset = {}
        dataset['omega'] = omega
        dataset['q'] = q
        dataset['intensity'] = intensity
        return dataset

    def save_pkl(self):
        with open(self.pklfile, 'wb') as f:
            pickle.dump(self.spectrab, f, -1)

    def load_pkl(self):
        with open(self.pklfile, 'rb') as f:
            self.spectrab = pickle.load(f)

    def get_org_spectra(self, qmin, qmax):
        self.get_qemap(qmin, qmax)
        self.spectrab = np.zeros((1, 2, self.ene.shape[0]))
        self.spectrab[0, 0, :] = self.ene
        self.spectrab[0, 1, :] = self.spectra
        self.get_all_data()
        self.spect(qmin, qmax, self.dataset)
        self.spectrab = np.concatenate((self.spectrab,
                                        self.spectra.reshape(1, 1, -1)),
                                       axis=1)



def run_org():
    prj = Sget_qlist(pklfile="run6202spectraorg.pkl")
    prj.get_org_data("0.0025", 6202)
    prj.get_qemap()
    prj.get_all_sdata()
    prj.spect(0.55, 0.70, prj.dataset, isplot=False)
    prj.spectrab = np.zeros((1, 2, prj.ene.shape[0]))
    prj.spectrab[0, 0, :] = prj.ene
    prj.spectrab[0, 1, :] = prj.spectra
    prj.save_pkl()


def check():
    import matplotlib.pyplot as plt
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


#run_org()
