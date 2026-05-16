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


def samplerun_on_dna_sdata(binw):
    fn = "0"+binw.split(".")[1]
    hws = np.linspace(-0.05, 0.15, 11)
    prj = Sget_qlist()
    prj.runno = 6207
    prj.fh = "./srlz/"+fn+"io/run"+str(prj.runno)+"united_"
    prj.TimeParam = "-1.0/1.0"
    for ihw in range(0, 10):
        prj.HwParam = binw + "/"+str(hws[ihw])+"/"+str(hws[ihw+1]-float(binw))
        prj.pklfile = "./srlz/"+fn+"/run"+str(prj.runno)+"_sdata"+str(ihw)+".pkl"
    prj.get_sdata()
    prj.get_all_sdata()
    prj.save_pkl()


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


samplerun_on_dna_sdata("0.000001")
#merge10_noselectq("0.000025")
#merge10("0.0001")
#samplerun_random_merge()

