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
import matplotlib.pyplot as plt
import pickle
import os
from scipy.interpolate import griddata
from get_qlist_nova_class import get_qlist as gq

m = 1.674927471*10**(-27)   # [kg]
h = 6.62607015*10**(-34)    # [J. s]
meVtoJ = 1.60218*10**(-22)  # [J/meV]
meVtoangsm2 = (1./0.81787)*0.01  # [Angs-2/meV]


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

    def get_qemapb(self, intensityb):
        self.DATB = Manyo.ElementContainerMatrix(self.DAT)
        ctp = Manyo.CppToPython()
        for ecaidx in range(0, self.DATB.PutSize()):
            for ecidx in range(0, self.DATB(0).PutSize()):
                vecy = ctp.ListToDoubleVector(intensityb[ecaidx, ecidx, :]
                                              .tolist())
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

    def load_pkl(self):
        with open(self.pklfile, 'rb') as f:
            self.spectrab = pickle.load(f)

    def get_boot_strap_sampled_spectra(self, nbs, qmin, qmax, seed=314,
                                       #restart=False, wnocorr=False,
                                       wnocorr=False,
                                       frac=None):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        psize = comm.Get_size()
        intensity1d = self.intensity.flatten().astype(int)
        nonzeroidx = np.nonzero(intensity1d)[0]
        x = np.array([idx for idx in nonzeroidx for num_repeat in
                     range(intensity1d[idx])], dtype=int)
        N = x.shape[0]
        np.random.seed(seed)
        intensityb = np.zeros_like(self.intensity)
        with open('randomstates.pkl.' + self.pklfile[3:7] + '.30000', 'rb') as f:
            randomstates = pickle.load(f)
        #if restart:
        if os.path.isfile(self.pklfile):
            #if rank == 0:
            #    print('restarting with randomstates.pkl')
            #with open('randomstates.pkl.' + self.pklfile[3:7] + '.30000', 'rb') as f:
            #    randomstates = pickle.load(f)
            #np.random.set_state(randomstates[-1])
            #randoffset = len(randomstates)
            #Nb = np.random.poisson(lam=N)
            #idx = np.random.randint(0, N, Nb)
            with open(self.pklfile, 'rb') as f:
                results = pickle.load(f)
            randoffset = results.shape[0]
        else:
            #randomstates = []
            randoffset = 0
        #for inb in range(nbs):
        #    randomstates.append(np.random.get_state())
        #    Nb = np.random.poisson(lam=N)
        #    idx = np.random.randint(0, N, Nb)
        #if rank == 0:
        #    print('writing randomstates.pkl')
        #    with open('randomstates.pkl', 'wb') as f:
        #        pickle.dump(randomstates, f, -1)
        for inb in range(rank*(nbs//psize), (rank+1)*(nbs//psize)):
            print(datetime.datetime.now(), 'chk1', inb)
        #for inb in range(nbs):
            intensityb *= 0.
            np.random.set_state(randomstates[inb+randoffset])
            if frac:
                Nb = np.random.poisson(lam=int(N*frac))
            else:
                Nb = np.random.poisson(lam=N)
            idx = np.random.randint(0, N, Nb)
            xb = x[idx]
            test_idxs = np.unravel_index(xb, self.intensity.shape)
            for _idx0, _idx1, _idx2 in zip(test_idxs[0], test_idxs[1],
                                           test_idxs[2]):
                intensityb[_idx0, _idx1, _idx2] += 1
            print(datetime.datetime.now(), 'chk2')
            self.get_qemapb(intensityb)
            print(datetime.datetime.now(), 'chk3')
            self.get_all_sdatab()
            #if inb == 0 and rank == 0:
            #    print("chk dataset['omega'] shape from qemap:",
            #          self.dataset['omega'].shape)
            #    print(self.dataset['omega'][0, 0:10])
            print(datetime.datetime.now(), 'chk4')
            self.spect(qmin, qmax, self.dataset, isplot=False)
            print(datetime.datetime.now(), 'chk5')
            if inb == rank*(nbs//psize):
                ener = np.zeros((nbs//psize, self.ene.shape[0]))
                spectrar = np.zeros((nbs//psize, self.ene.shape[0]))
            ener[inb - rank*nbs//psize, :] = self.ene[:]
            spectrar[inb - rank*nbs//psize, :] = self.spectra[:]
            print(datetime.datetime.now(), 'chk6')
            if wnocorr:
                # OLD
                # We use self.dataset doubly at present because self.spect
                # specifies this object name.
                # self.dataset = copy.deepcopy(self.datasetnocorr)
                # OLD
                self.datasetnocorr['intensity'] = intensityb
                print(datetime.datetime.now(), 'chk7')
                self.spect(qmin, qmax, self.datasetnocorr, isplot=False)
                print(datetime.datetime.now(), 'chk8')
                print('energy differences btw corr and nocorr:',
                      np.sum(self.ene - ener[inb - rank*nbs//psize, :]))
                if inb == rank*(nbs//psize):
                    enenocorrr = np.zeros_like(ener)
                    spectranocorrr = np.zeros_like(spectrar)
                enenocorrr[inb - rank*nbs//psize, :] = self.ene[:]
                spectranocorrr[inb - rank*nbs//psize, :] = self.spectra[:]
        ene1dt = np.zeros(ener.shape[1]*nbs)
        spectra1dt = np.zeros(spectrar.shape[1]*nbs)
        comm.Allgather(ener.flatten(), ene1dt)
        comm.Allgather(spectrar.flatten(), spectra1dt)
        self.spectrab = np.zeros((nbs, 2, ener.shape[1]))
        self.spectrab[:, 0, :] = ene1dt.reshape((nbs, -1))
        self.spectrab[:, 1, :] = spectra1dt.reshape((nbs, -1))
        if wnocorr:
            spectra1dt = np.zeros(spectranocorrr.shape[1]*nbs)
            comm.Allgather(spectranocorrr.flatten(), spectra1dt)
            self.spectrab = np.concatenate((self.spectrab,
                                            spectra1dt.reshape(nbs, 1, -1)),
                                           axis=1)
        #if restart:
        if os.path.isfile(self.pklfile):
            self.spectrab = np.concatenate((results, self.spectrab), axis=0)
        print(datetime.datetime.now(), 'chk9')


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
