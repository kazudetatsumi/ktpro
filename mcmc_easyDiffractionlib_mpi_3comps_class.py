#!/usr/bin/env python
import numpy as np
from mpi4py import MPI
from easyCore.Fitting.Fitting import Fitter
from easyDiffractionLib import Site, Phase, Phases
from easyDiffractionLib.sample import Sample as Job
from easyDiffractionLib.interface import InterfaceFactory as Calculator
from easyDiffractionLib.Jobs import Powder1DCW
from easyDiffractionLib.Jobs import Powder1DTOF
from easyDiffractionLib.Profiles.P1D import Instrument1DCWParameters\
        as CWParams
from easyDiffractionLib.Profiles.P1D import Instrument1DTOFParameters\
        as TOFParams
from easyDiffractionLib.Profiles.P1D import Powder1DParameters
from easyDiffractionLib.elements.Backgrounds.Point import\
        PointBackground, BackgroundPoint
import sys
import pprint
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if rank == 0:
    pprint.pprint(sys.path)


class mcmc():
    def __init__(self, crystfile, exptfile,  K=3, N_sampling=20000,
                 burn_in=10000, L=36, gamma=1.3, d_Mu=1.0, C_Mu=1.0,
                 sigma=0.1, min_Mu=2., max_Mu=20.):
        self.crystfile = crystfile
        self.exptfile = exptfile
        self.K = K
        self.N_sampling = N_sampling
        self.burn_in = burn_in
        self.L = L
        self.gamma = gamma
        self.d_Mu = d_Mu
        self.C_Mu = C_Mu
        if rank == 0:
            print("d =", d_Mu)
            print("C =", C_Mu)
            print("gamma =", gamma)
            print("L =", L)
        self.sigma = sigma
        self.min_Mu = min_Mu
        self.max_Mu = max_Mu
        self.init_easyDiffractionLib()
        self.N = self.meas_x.shape[0]
        np.random.seed(8)
        self.beta = np.zeros(L)
        self.beta[0] = 0
        self.step_Mu = np.zeros(self.L)
        for itemp in range(1, self.L):
            self.beta[itemp] = self.gamma**(itemp - (self.L - 1))
            if self.N*self.beta[itemp] <= 1:
                self.step_Mu[itemp] = self.C_Mu
            else:
                self.step_Mu[itemp] = self.C_Mu/(self.N*self.beta[itemp]
                                                 )**self.d_Mu
        self.count_accept_Mu = np.zeros(self.L)
        self.count_exchange = np.zeros(self.L)
        self.Mu_ar = np.zeros([self.L, self.K, self.N_sampling])
        self.MSE = np.zeros([self.L, self.N_sampling]) + 9999
        for itemp in range(0, L):
            self.Mu_ar[itemp, :, 0] = np.random.uniform(self.min_Mu,
                                                        self.max_Mu, self.K)

    def calc_E(self, _Mu):
        self.job.phases[0].cell.a = _Mu[0]
        self.job.phases[0].cell.b = _Mu[1]
        self.job.phases[0].cell.c = _Mu[2]
        out = (np.square((self.job.create_simulation(self.meas_x) -
                          self.meas_y)/self.meas_e)/self.meas_x.shape[0]).sum()
        return out

    def emcmc(self, isamp, itemp):
        self.Mu_ar[itemp, :, isamp] = self.Mu_ar[itemp, :, isamp-1]
        self.MSE[itemp, isamp] = self.MSE[itemp, isamp-1]
        pre_MSE = self.MSE[itemp, isamp]
        current_Mu = self.Mu_ar[itemp, :, isamp]
        NONNEG = False
        if itemp == 0:
            while not NONNEG:
                next_Mu = np.random.uniform(self.min_Mu, self.max_Mu, self.K)
                NONNEG = all(_Mu > 0. for _Mu in next_Mu)
        else:
            while not NONNEG:
                next_Mu = current_Mu + np.random.normal(0.0,
                                                        self.step_Mu[itemp],
                                                        self.K)
                NONNEG = all(_Mu > 0. for _Mu in next_Mu)
        next_MSE = self.calc_E(next_Mu)
        r = np.exp(self.N * self.beta[itemp] / (2.*self.sigma**2) *
                   (-next_MSE+pre_MSE))
        for ik in range(0, self.K):
            if next_Mu[ik] < self.min_Mu or next_Mu[ik] > self.max_Mu:
                r = -1
        if r >= np.random.uniform(0.0, 1.0) or itemp == 0 or isamp == 1:
            self.Mu_ar[itemp, :, isamp] = next_Mu
            self.MSE[itemp, isamp] = next_MSE
            self.count_accept_Mu[itemp] += 1
        if isamp == self.N_sampling - 1 and itemp == self.L-1:
            minMSE_sample = np.where(self.MSE[itemp] ==
                                     np.min(self.MSE[itemp]))[0][0]
            print("minMSE is found at the lattice constant of ",
                  self.Mu_ar[itemp, :, minMSE_sample],
                  self.calc_E(self.Mu_ar[itemp, :, minMSE_sample]))

    def rmc(self):
        for isamp in range(1, self.N_sampling):
            if isamp % 5 == 0 and rank == 0:
                print(isamp)
            if isamp == self.burn_in-1:
                self.count_accept_Mu = 0*self.count_accept_Mu
                self.count_exchange = 0*self.count_exchange
            for itemp in range(rank, self.L - ((self.L-rank-1) % size) + size
                               - 1, size):
                self.emcmc(isamp, itemp)
                istart = itemp - rank
                iend = istart+size
                self.MSE[istart:iend, isamp] =\
                    np.array(comm.allgather(self.MSE[itemp, isamp]))
                self.Mu_ar[istart:iend, :, isamp] =\
                    np.array(comm.allgather(self.Mu_ar[itemp, :, isamp])
                             ).reshape((size, self.K))
                self.count_accept_Mu[istart:iend] =\
                    np.array(comm.allgather(self.count_accept_Mu[itemp]))
            for itemp in range(0, self.L-1):
                r_exchange = np.exp(
                                    self.N/(2*self.sigma**2) *
                                    (self.beta[itemp+1] - self.beta[itemp]) *
                                    (self.MSE[itemp+1, isamp] -
                                     self.MSE[itemp, isamp])
                                    )
                if isamp % 2 == 0 and itemp % 2 == 0:
                    if r_exchange >= np.random.uniform(0.0, 1.0):
                        self.count_exchange[itemp] += 1
                        tmp = self.MSE[itemp, isamp]
                        self.MSE[itemp, isamp] = self.MSE[itemp+1, isamp]
                        self.MSE[itemp+1, isamp] = tmp
                        tmp_Mu = np.array(self.Mu_ar[itemp, :, isamp])
                        self.Mu_ar[itemp, :, isamp] =\
                            np.array(self.Mu_ar[itemp+1, :, isamp])
                        self.Mu_ar[itemp+1, :, isamp] = tmp_Mu
                elif isamp % 2 == 1 and itemp % 2 == 1:
                    if r_exchange >= np.random.uniform(0.0, 1.0):
                        self.count_exchange[itemp] += 1
                        tmp = self.MSE[itemp, isamp]
                        self.MSE[itemp, isamp] = self.MSE[itemp+1, isamp]
                        self.MSE[itemp+1, isamp] = tmp
                        tmp_Mu = np.array(self.Mu_ar[itemp, :, isamp])
                        self.Mu_ar[itemp, :, isamp] =\
                            np.array(self.Mu_ar[itemp + 1, :, isamp])
                        self.Mu_ar[itemp+1, :, isamp] = tmp_Mu
        if rank == 0:
            print("Accept rates for different beta")
            print(self.count_accept_Mu/(self.N_sampling - self.burn_in+1))
            print("Exchange rates for different beta")
            print(2*self.count_exchange/(self.N_sampling - self.burn_in+1))
            output = np.array([self.beta, self.count_accept_Mu /
                               (self.N_sampling - self.burn_in+1),
                              2*self.count_exchange /
                              (self.N_sampling - self.burn_in+1),
                              np.mean(self.N*self.MSE/self.sigma**2, axis=1)])
            np.savetxt("./results.txt", output.T)

    def init_easyDiffractionLib(self):
        phases = Phases.from_cif_file(self.crystfile)
        self.meas_x, self.meas_y, self.meas_e = np.loadtxt(self.exptfile,
                                                           unpack=True)
        if rank == 0:
            print("number of data points in meas_x", self.meas_x.shape[0])
        calculator = Calculator(interface_name='CrysPy')
        self.job = Powder1DCW("PbSO4", phases=phases, parameters=CWParams(),
                              interface=calculator)
        bkg = PointBackground(linked_experiment='PbSO4')
        bkg.append(BackgroundPoint.from_pars(self.meas_x[0], 200))
        bkg.append(BackgroundPoint.from_pars(self.meas_x[-1], 250))
        self.job.set_background(bkg)
        self.job.parameters.wavelength = 1.912
        self.job.pattern.scale.fixed = False
        self.job.pattern.zero_shift.fixed = False
        self.job.parameters.resolution_u.fixed = False
        self.job.parameters.resolution_v.fixed = False
        self.job.parameters.resolution_w.fixed = False
        self.job.backgrounds[0][0].y.fixed = False
        self.job.backgrounds[0][1].y.fixed = False
        fitter = Fitter(self.job, calculator.fit_func)
        result = fitter.fit(self.meas_x, self.meas_y, weights=1/self.meas_e,
                            method='least_squares',
                            minimizer_kwargs={'diff_step': 1e-5})
        if rank == 0:
            print(self.job.pattern.zero_shift)
            print(self.job.parameters.resolution_u)
            print(self.job.parameters.resolution_v)
            print(self.job.parameters.resolution_w)
            print(self.job.backgrounds[0][0])
            print(self.job.backgrounds[0][1])


def samplerun():
    crystfile = "/home/kazu/desktop/231207/easyDiffraction/examples/PbSO4.cif"
    exptfile = "/home/kazu/desktop/231207/easyDiffraction/examples/D1A@ILL.xye"
    prj = mcmc(crystfile, exptfile)
    prj.rmc()
    #prj.init_easyDiffractionLib()


#samplerun()
