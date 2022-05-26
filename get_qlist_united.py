#!/usr/bin/env python
import Cmm
import sys
import get_qlist_nova_class as gqnc
import numpy as np
import matplotlib.pyplot as plt

m = 1.674927471*10**(-27)   # [kg]
h = 6.62607015*10**(-34)    # [J. s]
meVtoJ = 1.60218*10**(-22)  # [J/meV]
meVtoangsm2 = (1./0.81787)*0.01  # [Angs-2/meV]


class Sget_qlist(gqnc.get_qlist, object):
    def __init__(self, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile

        super(Sget_qlist, self).__init__()

    def get_data(self):
        EC = Cmm.GetHistogramMon(
                runNo=self.runno, useEiConv=True, LambdaParam="6.321/4.15",
                t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                TimeParam=self.TimeParam, UseFastChopper=True,
                isHistogram=False)
        Cmm.MutiplyConstant(dat=EC, factor=1e-09)
        self.DAT = Cmm.GetHistogramHW(
                runNo=self.runno, HwParam=self.HwParam,
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                useT0ModCorr=False, TimeParam=self.TimeParam,
                UseFastChopper=True, tofOffsetFile="none", isHistogram=False)

    def get_sdata(self):
        EC = Cmm.GetHistogramMon(
                runNo=self.runno, useEiConv=True, LambdaParam="6.321/4.15",
                t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                TimeParam=self.TimeParam, UseFastChopper=True,
                isHistogram=False)
        Cmm.MutiplyConstant(dat=EC, factor=1e-09)
        DAT = Cmm.GetHistogramHW(
                runNo=self.runno, HwParam=self.HwParam,
                LambdaParam="6.321/4.15", t0_offset=12325.0,
                useT0ModCorr=False, TimeParam=self.TimeParam,
                UseFastChopper=True, tofOffsetFile="none", isHistogram=False)
        Cmm.DoMask(dat=DAT, filename="maskTY.txt")
        ECM = Cmm.ILambdaCorrDNA(dat=DAT, ec=EC, useMonEff=True)
        ECM2 = Cmm.SolidAngleCorrDNA(
                dat=ECM, useDetEff=True, useAbsoCorr=False, useEffCorr=False,
                sampletype="sample", sampleDataPath="test_sample_data.dat",
                DetEffDataPath="none")
        Cmm.MutiplyConstant(dat=ECM2, factor=1e-06)
        self.DAT = Cmm.CreateQEMap(dat=ECM2, startQ=0.0, endQ=2.0, deltaQ=0.05)

    def run_one_thousandth(self):
        for itime in range(750, 875):
            self.TimeParam = str(itime/1000.0 * 5842.0) +\
                             ", " + str((itime+1)/1000.0 * 5842.0)
            qmin = 0.55
            qmax = 0.70
            spectrafile = str(self.fh)+"spectra_"+str(itime)+".pkl"
            sspectrafile = str(self.fh)+"sspectra_"+str(itime)+".pkl"
            self.get_data()
            self.get_all_data()
            self.spect(qmin, qmax)
            self.save_spectra(spectrafile)
            self.get_sdata()
            self.get_all_sdata()
            self.spect(qmin, qmax)
            self.save_spectra(sspectrafile)

    def run_one_oneth(self, onlys=False):
        qmin = 0.55
        qmax = 0.70
        spectrafile = str(self.fh)+"spectra.pkl"
        sspectrafile = str(self.fh)+"sspectra.pkl"
        if not onlys:
            self.get_data()
            self.get_all_data()
            self.spect(qmin, qmax)
            self.save_spectra(spectrafile)
        self.get_sdata()
        self.get_all_sdata()
        self.spect(qmin, qmax)
        self.save_spectra(sspectrafile)

    def merge(self):
        for ibundle in range(0, 200):
            numbers = np.arange(0, 5)
            numbers = np.random.permutation(numbers) + ibundle*5
            self.take_data(numbers)
        for islot in range(0, 5):
            if islot == 0:
                self.spectra = self.sampled_spectra[0]
            else:
                self.spectra += self.sampled_spectra[islot]
            self.save_spectra(str(self.fh)+"spectra_"+str(islot)+"of5.pkl")

    def take_data(self, numbers):
        if np.min(numbers) == 0:
            self.initialize()
        for islot, sfidx in enumerate(numbers):
            self.pklfile = str(self.fh)+"spectra_"+str(sfidx)+".pkl"
            self.read_pkl()
            self.sampled_spectra[islot, :] += self.dataset['spectra']

    def initialize(self):
        self.pklfile = str(self.fh)+"spectra_0.pkl"
        self.read_pkl()
        self.sampled_spectra = np.zeros((5, self.dataset['spectra'].shape[0]))
        self.ene = self.dataset['energy']

    def checkspec(self):
        for islot in [0, 4]:
            self.pklfile = str(self.fh)+"spectra_"+str(islot)+"of5.pkl"
            self.read_pkl()
            plt.plot(self.dataset['energy'], self.dataset['spectra'] /
                     np.sum(self.dataset['spectra']))
        plt.yscale('log')
        plt.show()


def samplerun_on_dna():
    prj = Sget_qlist()
    prj.runno = 7784
    prj.fh = "./srlz/0000025new/run"+str(prj.runno)+"united"
    #prj.fh = "./srlz/tmp/run"+str(prj.runno)+"united_0125"
    prj.TimeParam = "-1.0/-1.0"
    #prj.TimeParam = "0.0, 8764.0"
    #prj.TimeParam = "0.0, 4382.0"
    #prj.TimeParam = "0.0, 1461.0"
    #prj.TimeParam = "0.0, 5843.0"
    #prj.TimeParam = "0.0, 10225.0"
    #prj.TimeParam = "0.0, 7303.0"
    prj.HwParam = "0.000025/-0.05/0.15"
    prj.run_one_oneth()


def samplerun_random_merge():
    prj = Sget_qlist()
    prj.runno = 7784
    prj.fh = "./srlz/run"+str(prj.runno)
    prj.merge()
    prj.fh = "./srlz/run"+str(prj.runno)+"s"
    prj.merge()

def samplerun_on_dna_all(binw):
    fn = "0"+binw.split(".")[1]
    #fts = ["", "075_", "0625_", "05_", "0375_", "025_", "0125_"]
    #tps = ["-1.0/-1.0", "0.0, 8764.0", "0.0, 7303.0", "0.0, 5843.0", "0.0, 4382.0", "0.0, 2921.5", "0.0, 1460.7"] # for 62o2
    fts = [""]
    tps = ["-1.0/-1.0"] 
    for ft, tp in zip(fts, tps):
        prj = Sget_qlist()
        prj.runno = 6203
        prj.HwParam = binw+"/-0.05/0.15"
        prj.fh = "./srlz/"+fn+"s/run"+str(prj.runno)+"united_"+ft
        prj.TimeParam = tp
        prj.run_one_oneth(onlys=True)


def check():
    prj = Sget_qlist()
    prj.runno = 7784
    prj.fh = "./srlz/run"+str(prj.runno)
    prj.checkspec()
    prj.fh = "./srlz/run"+str(prj.runno)+"s"
    prj.checkspec()


#samplerun_after_dna()
#samplerun_random_merge()
#check()


#samplerun_on_dna()
samplerun_on_dna_all("0.0001")
#samplerun_random_merge()

