#!/usr/bin/env python
# This script reads a set of raw neutron count data of DNA without any
# corrections, i.e., an uncorrected energy container matrix using Manyo
# library, samples counts from the data in a bootstrap manner and applies the
# necessary corrections to draw bootstrap sampled QENS data corresponding to
# double differentia cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
try:
    import Cmm
    import Manyo
except ModuleNotFoundError as err:
    if rank == 0:
        print(err)
import numpy as np
import datetime
#import matplotlib.pyplot as plt
import pickle
import os
#from scipy.interpolate import griddata
from get_qlist_nova_class import get_qlist as gq

#m = 1.674927471*10**(-27)   # [kg]
#h = 6.62607015*10**(-34)    # [J. s]
#meVtoJ = 1.60218*10**(-22)  # [J/meV]
#meVtoangsm2 = (1./0.81787)*0.01  # [Angs-2/meV]


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

    def get_qemapb(self, intensityb, qmin, qmax):
        #DATB = Manyo.ElementContainerMatrix(self.DAT) 
        # To reduce memory usage, self.DAT is overwritten.
        ctp = Manyo.CppToPython()
        for ecaidx in range(0, self.DAT.PutSize()):
            for ecidx in range(0, self.DAT(0).PutSize()):
                vecy = ctp.ListToDoubleVector(intensityb[ecaidx, ecidx, :]
                                              .tolist())
                vece = ctp.ListToDoubleVector((intensityb[ecaidx,
                                               ecidx, :]**0.5).tolist())
                self.DAT(ecaidx, ecidx).Replace("Intensity", vecy)
                self.DAT(ecaidx, ecidx).Replace("Error", vece)
                self.DAT(ecaidx, ecidx).SetKeys("EnergyTransfer", "Intensity",
                                                "Error")
        Cmm.DoMask(dat=self.DAT, filename="maskTY.txt")
        ECM = Cmm.ILambdaCorrDNA(dat=self.DAT, ec=self.EC, useMonEff=True)
        ECM2 = Cmm.SolidAngleCorrDNA(
                dat=ECM, useDetEff=True, useAbsoCorr=False, useEffCorr=False,
                sampletype="sample", sampleDataPath="test_sample_data.dat",
                DetEffDataPath="none")
        Cmm.MutiplyConstant(dat=ECM2, factor=1e-06)
        DATBQE = Cmm.CreateQEMap(dat=ECM2, startQ=0.0, endQ=2.0, deltaQ=0.05)
        dataset = self.get_all_sdatab(DATBQE)
        self.spect2(qmin, qmax, dataset, isplot=False)

    def get_all_sdatab(self, DATBQE):
        q = np.zeros((DATBQE.PutSize(), len(DATBQE(0).PutYList())))
        omega = np.zeros((DATBQE.PutSize(), len(DATBQE(0).PutYList())
                          ))
        intensity = np.zeros((DATBQE.PutSize(),
                              len(DATBQE(0).PutYList())))
        ones = np.ones((len(DATBQE(0).PutYList())))
        for ecidx in range(0, DATBQE.PutSize()):
            omega[ecidx, :] = np.array(DATBQE(ecidx).PutXList()[:-1])
            q[ecidx, :] = ones*DATBQE(ecidx).PutHeader().\
                PutDoubleVector('XRANGE')[0]
            intensity[ecidx, :] = np.array(DATBQE(ecidx).PutYList())
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

    def get_boot_strap_sampled_spectra(self, nbs, qmin, qmax, seed=314,
                                       wnocorr=False, frac=None):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        psize = comm.Get_size()
        #intensity1d = self.orgintensity.flatten().astype(int)
        nonzeroidx = np.nonzero(self.orgintensity1d)[0]
        intdtype = self.orgintensity1d.dtype
        x = np.array([idx for idx in nonzeroidx for num_repeat in
                     range(self.orgintensity1d[idx])], dtype=self.get_intdtype(np.max(nonzeroidx)))
        N = x.shape[0]
        np.random.seed(seed)
        intensityb = np.zeros(self.orgintensityshape, dtype=intdtype)
        with open('randomstates.pkl.' + self.pklfile[3:7] + '.30000', 'rb'
                  ) as f:
            randomstates = pickle.load(f)
        #if os.path.isfile(self.pklfile):
        #    if rank == 0:
        #        with open(self.pklfile, 'rb') as f:
        #            results = pickle.load(f)
        #        randoffset = results.shape[0]
        #        comm.Barrier()
        #        comm.bcast(randoffset, root=0)
        if os.path.isfile(self.pklfile):
            with open(self.pklfile, 'rb') as f:
                results = pickle.load(f)
                randoffset = results.shape[0]
        else:
            randoffset = 0
        for inb in range(rank*(nbs//psize), (rank+1)*(nbs//psize)):
            #print(datetime.datetime.now(), 'chk1', inb)
            intensityb *= 0
            np.random.set_state(randomstates[inb+randoffset])
            if frac:
                Nb = np.random.poisson(lam=int(N*frac))
            else:
                Nb = np.random.poisson(lam=N)
            idx = np.random.randint(0, N, Nb)
            xb = x[idx]
            test_idxs = np.unravel_index(xb, self.orgintensityshape)
            for _idx0, _idx1, _idx2 in zip(test_idxs[0], test_idxs[1],
                                           test_idxs[2]):
                intensityb[_idx0, _idx1, _idx2] += 1
            #print(datetime.datetime.now(), 'chk2')
            self.get_qemapb(intensityb, qmin, qmax)
            #print(datetime.datetime.now(), 'chk345')
            if inb == rank*(nbs//psize):
                ener = np.zeros((nbs//psize, self.ene.shape[0]))
                spectrar = np.zeros((nbs//psize, self.ene.shape[0]))
            ener[inb - rank*nbs//psize, :] = self.ene[:]
            spectrar[inb - rank*nbs//psize, :] = self.spectra[:]
            #print(datetime.datetime.now(), 'chk6')
            if wnocorr:
                self.dataset['intensity'] = intensityb
                print(datetime.datetime.now(), 'chk7')
                self.spect2(qmin, qmax, self.dataset, isplot=False)
                print(datetime.datetime.now(), 'chk8')
                #print('energy differences btw corr and nocorr:',
                #      np.sum(self.ene - ener[inb - rank*nbs//psize, :]))
                if inb == rank*(nbs//psize):
                    #enenocorrr = np.zeros_like(ener)
                    spectranocorrr = np.zeros_like(spectrar)
                #enenocorrr[inb - rank*nbs//psize, :] = self.ene[:]
                spectranocorrr[inb - rank*nbs//psize, :] = self.spectra[:]
        ene1dt = np.zeros(ener.shape[1]*nbs)
        spectra1dt = np.zeros(spectrar.shape[1]*nbs)
        ene1dt = comm.gather(ener.flatten(), root=0)
        spectra1dt = comm.gather(spectrar.flatten(), root=0)
        if rank == 0:
            self.spectrab = np.zeros((nbs, 2, ener.shape[1]))
            self.spectrab[:, 0, :] = np.array(ene1dt).reshape((nbs, -1))
            self.spectrab[:, 1, :] = np.array(spectra1dt).reshape((nbs, -1))
        if wnocorr:
            spectra1dt = comm.gather(spectranocorrr.flatten(), root=0)
            if rank == 0:
                self.spectrab = np.concatenate((self.spectrab,
                                                np.array(spectra1dt).
                                                reshape(nbs, 1, -1)),
                                               axis=1)
        if rank == 0 and 'results' in locals():
            self.spectrab = np.concatenate((results, self.spectrab), axis=0)
        #print(datetime.datetime.now(), 'chk9')

