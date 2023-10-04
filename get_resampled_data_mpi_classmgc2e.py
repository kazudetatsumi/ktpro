#!/usr/bin/env python
# This script reads a set of raw neutron count data of DNA without any
# corrections, i.e., an uncorrected energy container matrix using Manyo
# library, samples counts from the data in a bootstrap manner and applies the
# necessary corrections to draw bootstrap sampled QENS data corresponding to
# double differentia cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
import numpy as np
import datetime
import pickle
import os
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()
from get_resampled_data_mpi_class7 import Sget_qlist as gq


class Sget_qlist(gq):
    def __init__(self, randfile=None, pklfile=None):
        self.randfile = randfile
        self.pklfile = pklfile

    def get_boot_strap_sampled_spectra(self, nbs, qmin, qmax, seed=314,
                                       wnocorr=False, frac=None,
                                       maskfile="maskTY140218ForAfterRun52.txt"
                                       ):
        import gc
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        psize = comm.Get_size()
        nonzeroidx = np.nonzero(self.orgintensity1d)[0]
        intdtype = self.orgintensity1d.dtype
        x = np.array([idx for idx in nonzeroidx for num_repeat in
                     range(self.orgintensity1d[idx])],
                     dtype=self.get_intdtype(np.max(nonzeroidx)))
        N = x.shape[0]
        np.random.seed(seed)
        intensityb = np.zeros(self.orgintensityshape, dtype=intdtype)
        with open(self.randfile, 'rb') as f:
            randomstates = pickle.load(f)
        if os.path.isfile(self.pklfile):
            if rank == 0:
                with open(self.pklfile, 'rb') as f:
                    results = pickle.load(f)
                    randoffset = results.shape[0]
                results = None
                del results
                gc.collect()
            else:
                randoffset = None
            comm.barrier()
            randoffset = comm.bcast(randoffset, root=0)
        else:
            randoffset = 0
        for inb in range(rank*(nbs//psize), (rank+1)*(nbs//psize)):
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
            self.spectme(qmin, qmax, self.get_qemapbe(intensityb, qmin, qmax,
                                                      maskfile=maskfile))
            if inb == rank*(nbs//psize):
                ener = np.zeros((nbs//psize, self.ene.shape[0]))
                spectrar = np.zeros((nbs//psize, self.ene.shape[0]))
                errr = np.zeros((nbs//psize, self.ene.shape[0]))
            ener[inb - rank*nbs//psize, :] = self.ene[:]
            spectrar[inb - rank*nbs//psize, :] = self.spectra[:]
            errr[inb - rank*nbs//psize, :] = self.err[:]
            if wnocorr:
                self.dataset['intensity'] = intensityb
                self.spectm(qmin, qmax, self.dataset)
                print(datetime.datetime.now(), 'chk8', rank)
                if inb == rank*(nbs//psize):
                    spectranocorrr = np.zeros_like(spectrar)
                spectranocorrr[inb - rank*nbs//psize, :] = self.spectra[:]
        del x, self.dataset, intensityb
        gc.collect()
        ene1dt = np.zeros(ener.shape[1]*nbs)
        spectra1dt = np.zeros(spectrar.shape[1]*nbs)
        err1dt = np.zeros(spectrar.shape[1]*nbs)
        ene1dt = comm.gather(ener.flatten(), root=0)
        spectra1dt = comm.gather(spectrar.flatten(), root=0)
        err1dt = comm.gather(errr.flatten(), root=0)
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
        if rank == 0:
            self.spectrab = np.concatenate((self.spectrab,
                                            np.array(err1dt).
                                            reshape(nbs, 1, -1)),
                                           axis=1)
            _sh = self.spectrab.shape
            self.spectrab = self.spectrab.reshape((_sh[0], _sh[1], len(qmin),
                                                  -1))
        if randoffset > 0:
            if rank == 0:
                with open(self.pklfile, 'rb') as f:
                    results = pickle.load(f)
                self.spectrab = np.concatenate((results, self.spectrab),
                                               axis=0)
                del results
                gc.collect()
        print(datetime.datetime.now(), 'chk9', rank)
