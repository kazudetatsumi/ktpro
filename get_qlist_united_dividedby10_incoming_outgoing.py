#!/usr/bin/env python
#import Cmm
import sys
import get_qlist_nova_class as gqnc
import numpy as np
import matplotlib.pyplot as plt
import time

m = 1.674927471*10**(-27)   # [kg]
h = 6.62607015*10**(-34)    # [J. s]
meVtoJ = 1.60218*10**(-22)  # [J/meV]
meVtoangsm2 = (1./0.81787)*0.01  # [Angs-2/meV]


class Sget_qlist(gqnc.get_qlist, object):
    def __init__(self, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile

        super(Sget_qlist, self).__init__(save_file=save_file, pklfile=pklfile)

    def get_monidata(self):
        self.EC = Cmm.GetHistogramMon(
                runNo=self.runno, useEiConv=True, LambdaParam="6.321/4.15",
                t0_offset=12325.0, background=0.0, useT0ModCorr=False,
                TimeParam=self.TimeParam, UseFastChopper=True,
                isHistogram=False)

    def get_data(self):
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

    #def get_all_monidata(self):
    #    self.energy =  self.EC.PutXList()
    #    self.spectra =  self.EC.PutYList()
    def get_all_monidata(self):
        self.spectra = np.zeros((1, 3, len(self.EC.PutYList())))
        self.spectra[0, 0, :] = self.EC.PutXList()[:-1]
        self.spectra[0, 1, :] = self.EC.PutYList()
        self.spectra[0, 2, :] = self.EC.PutEList()

    def run_moni(self):
        monispectrafile = str(self.fh)+"monispectra.pkl"
        self.get_monidata()
        self.get_all_monidata()
        self.save_spectra(monispectrafile, old=True)

    def run_one_oneth(self):
        qmin = 0.55
        qmax = 0.70
        spectrafile = str(self.fh)+"spectra.pkl"
        self.get_data()
        self.get_all_data()
        self.spect(qmin, qmax)
        self.save_spectra(spectrafile)

    def run_one_oneth_noselectq(self):
        print(' before get_data', time.ctime())
        self.get_data()
        print(' before get_all_data', time.ctime())
        self.get_all_data()
        for i in range(0, 5):
            qmin = 0.1 + 0.15*i
            qmax = qmin + 0.15
            spectrafile = str(self.fh)+"_"+str(qmin)+"-"+str(qmax)+"_spectra.pkl"
            print(' before spect', time.ctime())
            self.spect(qmin, qmax)
            self.save_spectra(spectrafile)

    def run_one_oneth_noq(self):
        print(' before get_data', time.ctime())
        self.get_data()
        print(' before get_all_data', time.ctime())
        self.get_all_data()
        self.pklfile=str(self.fh)+"_data.pkl"
        self.save_pkl()

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


def samplerun_on_dna(binw):
    fn = "0"+binw.split(".")[1]
    hws = np.linspace(-0.05, 0.15, 11)
    print(hws)
    #fts = ["", "0875_", "075_", "0625_", "05_", "0375_", "025_", "0125_"]
    #fts = ["075_", "0625_", "05_", "0375_"]
    #fts = ["0625_"]
    fts = [""]
    #fts = ["0875_", "075_", "0625_", "05_", "0375_", "025_", "0125_"]
    #fts = ["05_", "0375_"]  # for 62o4 lastio
    #fts = ["0125_", "025_"]  # for 62o4 lastio2
    #tps = ["-1.0/-1.0", "0.0, 10225.0", "0.0, 8764.0", "0.0, 7303.0", "0.0, 5843.0", "0.0, 4382.0", "0.0, 2921.5", "0.0, 1460.7"] # for 62o2
    #tps = ["0.0, 19354.0", "0.0, 16589.0", "0.0, 13824.0", "0.0, 11059.0", "0.0, 8294.0", "0.0, 5529.0", "0.0, 2765.0"] #  for 778_4
    #tps = ["0.0, 10225.0", "0.0, 8764.0", "0.0, 7303.0", "0.0, 5843.0", "0.0, 4382.0"]
    #tps = ["0.0, 2921.5", "0.0, 1460.8"]
    #tps = ["0.0, 8764.0", "0.0, 7303.0", "0.0, 5843.0", "0.0, 4382.0"]
    #tps = ["0.0, 2921.5", "0.0, 1460.7"]
    #tps = ["8763.0, 11684.0", "10225.0, 11684.0"] #for 778_4 lastio
    #tps = ["0.0, 7303.0"]                       
    #tps = ["5110.0, 10220.0", "6387.0, 10220.0"]  #for 62o2 lastio
    #tps = ["7523.0, 10220.0", "8764.0, 10220.0"]  #for 62o2 lastio2
    tps = ["-1.0/-1.0"]
    #tps = ["0.0, 21558.0"]
    for ft, tp in zip(fts, tps):
        prj = Sget_qlist()
        prj.runno = 6207
        prj.fh = "./srlz/"+fn+"io/run"+str(prj.runno)+"united_"+ft
        #prj.fh = "./srlz/"+fn+"lastio/run"+str(prj.runno)+"united_"+ft
        prj.TimeParam = tp
        prj.run_moni()
        for ihw in  range(0, 10):
            prj.HwParam = binw + "/"+str(hws[ihw])+"/"+str(hws[ihw+1]-float(binw))
            prj.fh = "./srlz/"+fn+"io/run"+str(prj.runno)+"united_"+ft+str(ihw)
            #prj.fh = "./srlz/"+fn+"lastio/run"+str(prj.runno)+"united_"+ft+str(ihw)
            print(prj.HwParam)
            prj.run_one_oneth()

def samplerun_on_dna_noselectq(binw):
    fn = "0"+binw.split(".")[1]
    hws = np.linspace(-0.05, 0.15, 11)
    print(hws)
    fts = [""]
    tps = ["-1.0/-1.0"]
    for ft, tp in zip(fts, tps):
        prj = Sget_qlist()
        prj.runno = 6207
        prj.TimeParam = tp
        for ihw in  range(0, 10):
            prj.HwParam = binw + "/"+str(hws[ihw])+"/"+str(hws[ihw+1]-float(binw))
            prj.fh = "./srlz/"+fn+"io/run"+str(prj.runno)+"_"+ft+str(ihw)
            print(prj.HwParam)
            prj.run_one_oneth_noselectq()

def samplerun_on_dna_noq(binw):
    fn = "0"+binw.split(".")[1]
    hws = np.linspace(-0.05, 0.15, 11)
    print(hws)
    fts = [""]
    tps = ["-1.0/-1.0"]
    for ft, tp in zip(fts, tps):
        prj = Sget_qlist()
        prj.runno = 6207
        prj.TimeParam = tp
        for ihw in  range(0, 1):
            prj.HwParam = binw + "/"+str(hws[ihw])+"/"+str(hws[ihw+1]-float(binw))
            prj.fh = "./srlz/"+fn+"io/run"+str(prj.runno)+"_"+ft+str(ihw)
            prj.run_one_oneth_noq()

def merge10(binw):
    fn = "0"+binw.split(".")[1]
    hws = np.linspace(-0.05, 0.15, 11)
    nf = 10
    prj = Sget_qlist()
    prj.runno = 6207
    #fh = "./srlz/"+fn+"io/run"+str(prj.runno)+"united_0375_"
    fh = "./srlz/"+fn+"io/run"+str(prj.runno)+"united_"
    #fh = "./srlz/"+fn+"lastio/run"+str(prj.runno)+"united_0125_"
    endterms = ["spectra.pkl", "sspectra.pkl"]
    endterms = ["spectra.pkl"]
    for term in endterms:
        for ihw in range(0, nf):
            print(ihw)
            prj.fh = fh + str(ihw)
            prj.pklfile = str(prj.fh)+term
            prj.read_pkl()
            print(prj.dataset["energy"].shape)
            print(prj.dataset["energy"][0:3])
            print(prj.dataset["energy"][-3:])
            if ihw == 0:
                prj.ene = prj.dataset["energy"]
                prj.spectra = prj.dataset["spectra"]
            else:
                prj.ene = np.concatenate((prj.ene, prj.dataset["energy"]))
                prj.spectra = np.concatenate((prj.spectra, prj.dataset["spectra"]))
        plt.plot(prj.ene, prj.spectra)
        plt.show()
        prj.dataset["energy"] = prj.ene
        prj.dataset["spectra"] = prj.spectra
        prj.pklfile = fh[:-1]+term
        prj.save_pkl()

def merge10_noselectq(binw):
    fn = "0"+binw.split(".")[1]
    hws = np.linspace(-0.05, 0.15, 11)
    nf = 10
    prj = Sget_qlist()
    prj.runno = 6207
    for i in range(0, 5):
        qmin = 0.1 + 0.25*i
        qmax = qmin + 0.25
        for ihw in range(0, nf):
            prj.fh = "./srlz/"+fn+"io/run"+str(prj.runno)+"_"+str(ihw)
            prj.pklfile = str(prj.fh)+"_"+str(qmin)+"-"+str(qmax)+"_spectra.pkl"
            prj.read_pkl()
            if ihw == 0:
                prj.ene = prj.dataset["energy"]
                prj.spectra = prj.dataset["spectra"]
            else:
                prj.ene = np.concatenate((prj.ene, prj.dataset["energy"]))
                prj.spectra = np.concatenate((prj.spectra, prj.dataset["spectra"]))
        plt.plot(prj.ene, prj.spectra)
        plt.show()
        prj.dataset["energy"] = prj.ene
        prj.dataset["spectra"] = prj.spectra
        prj.pklfile = str(prj.fh)[:-1]+str(qmin)+"-"+str(qmax)+"_spectra.pkl"
        print(prj.pklfile)
        prj.save_pkl()

def samplerun_random_merge():
    prj = Sget_qlist()
    prj.runno = 6207
    prj.fh = "./srlz/run"+str(prj.runno)
    prj.merge()
    prj.fh = "./srlz/run"+str(prj.runno)+"s"
    prj.merge()


def check():
    prj = Sget_qlist()
    prj.runno = 6207
    prj.fh = "./srlz/run"+str(prj.runno)
    prj.checkspec()
    prj.fh = "./srlz/run"+str(prj.runno)+"s"
    prj.checkspec()

def check2():
    save_file = "./srlz/run6207_half.srlz"
    pklfile = "./srlz/run6207_half.pkl"
    spectrafile = "./srlz/run6207spectra.pkl"
    qmin = 0.55
    qmax = 0.70
    prj = Sget_qlist(save_file=save_file, pklfile=pklfile)
    prj.HwParam = "0.00025/-0.05/0.15"
    prj.read_pkl()
    prj.spect(qmin, qmax)


#samplerun_after_dna()
#samplerun_random_merge()
#check2()


#samplerun_on_dna("0.0001")
#samplerun_on_dna_noselectq("0.000025")
samplerun_on_dna_noq("0.000001")
#merge10_noselectq("0.000025")
#merge10("0.0001")
#samplerun_random_merge()

