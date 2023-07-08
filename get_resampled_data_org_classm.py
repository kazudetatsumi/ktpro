#!/usr/bin/env python
# This script was folked from get_resampled_data_mpi_class.py, only having
# the part not related to the bootstrap sampling.
# This script reads a set of raw neutron count data of DNA without any
# corrections, i.e., an uncorrected energy container matrix using utsusemi DNA
# facede Cmm, as well as applies usual corrections on the raw data.
# The uncorrected and corrected  spectral data are stored as runXXXspectrab.pkl
# to utilize the analysis scripts for bootstrap sampled data sets.
# Kazuyoshi TATSUMI 2023/06/07
# Modified specific routines for multiple of qs. 
# Kazuyoshi TATSUMI 2023/06/12

try:
    import Cmm
except ModuleNotFoundError as err:
    print(err)
import numpy as np
# from get_qlist_nova_class import get_qlist as gq
from get_resampled_data_org_class import Sget_qlist as gq


class Sget_qlist(gq):
    def __init__(self, save_file=None, pklfile=None):
        self.save_file = save_file
        self.pklfile = pklfile

    def get_qemap(self, qmin, qmax):
        #Cmm.DoMask(dat=self.DAT, filename="maskTY.txt")
        Cmm.DoMask(dat=self.DAT, filename="maskTY140218ForAfterRun52.txt")
        ECM = Cmm.ILambdaCorrDNA(dat=self.DAT, ec=self.EC, useMonEff=True)
        ECM2 = Cmm.SolidAngleCorrDNA(
                dat=ECM, useDetEff=True, useAbsoCorr=False, useEffCorr=False,
                sampletype="sample", sampleDataPath="test_sample_data.dat",
                DetEffDataPath="none")
        Cmm.MutiplyConstant(dat=ECM2, factor=1e-06)
        #DATQE = Cmm.CreateQEMap(dat=ECM2, startQ=0.0, endQ=2.0, deltaQ=0.05)
        DATQE = Cmm.CreateQEMap(dat=ECM2, startQ=0.074, endQ=1.85, deltaQ=0.148)
        dataset = self.get_all_sdata(DATQE)
        self.spectm(qmin, qmax, dataset)

    def get_org_spectra(self, qmin, qmax):
        self.get_qemap(qmin, qmax)
        self.spectrab = np.zeros((1, 2, self.ene.shape[0]))
        self.spectrab[0, 0, :] = self.ene
        self.spectrab[0, 1, :] = self.spectra
        self.get_all_data3()
        self.spectm(qmin, qmax, self.dataset)
        self.spectrab = np.concatenate((self.spectrab,
                                        self.spectra.reshape(1, 1, -1)),
                                       axis=1)
        _sh = self.spectrab.shape
        self.spectrab = self.spectrab.reshape((_sh[0], _sh[1], len(qmin), -1))
