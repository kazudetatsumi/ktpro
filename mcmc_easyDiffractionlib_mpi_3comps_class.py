#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import os
from mpi4py import MPI
from easyCore import np
from easyCore.Fitting.Fitting import Fitter
from easyDiffractionLib import Site, Phase, Phases
from easyDiffractionLib.sample import Sample as Job
from easyDiffractionLib.interface import InterfaceFactory as Calculator
from easyDiffractionLib.Jobs import Powder1DCW
from easyDiffractionLib.Jobs import Powder1DTOF
from easyDiffractionLib.Profiles.P1D import Instrument1DCWParameters as CWParams
from easyDiffractionLib.Profiles.P1D import Instrument1DTOFParameters as TOFParams
from easyDiffractionLib.Profiles.P1D import Powder1DParameters
from easyDiffractionLib.elements.Backgrounds.Point import PointBackground, BackgroundPoint


import sys
import pprint

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
if rank == 0:
    pprint.pprint(sys.path)


class mcmc():
# Test for optimization on the lattice constant
#    def __init__(self, datafile, K=3, N_sampling=20000, burn_in=10000, L=40,
#                 gamma=1.3, d_Mu=1.0, C_Mu=1.0, d_S=1.0, C_S=0.6, d_W=0.9,
#                 C_W=1.0, sigma=0.1, min_Mu=0., max_Mu=2.5, min_S=0.01,
#                 max_S=1.5, min_W=0., max_W=2.0):
    #def __init__(self, crystfile, exptfile,  K=1, N_sampling=20000, burn_in=10000, L=40,
    def __init__(self, crystfile, exptfile,  K=3, N_sampling=20000, burn_in=10000, L=36,
                 gamma=1.3, d_Mu=1.0, C_Mu=1.0, sigma=0.1, min_Mu=2.,
                 max_Mu=20.):
        self.crystfile = crystfile
        self.exptfile = exptfile
        self.K = K
        self.N_sampling = N_sampling
        self.burn_in = burn_in
        self.L = L
        self.gamma = gamma
        self.d_Mu = d_Mu
        self.C_Mu = C_Mu
#        self.d_S = d_S
#        self.C_S = C_S
#        self.d_W = d_W
#        self.C_W = C_W
        self.sigma = sigma
        self.min_Mu = min_Mu
        self.max_Mu = max_Mu
#        self.min_S = min_S
#        self.max_S = max_S
#        self.min_W = min_W
#        self.max_W = max_W
        #rhochi = cryspy.load_file(self.datafile)
        #self.rhochi_dict = rhochi.get_dictionary()
        self.init_easyDiffractionLib()
        self.N = self.meas_x.shape[0]
        np.random.seed(8)
        self.beta = np.zeros(L)
        self.beta[0] = 0
        self.step_Mu = np.zeros(self.L)
#        self.step_S = np.zeros(self.L)
#        self.step_W = np.zeros(self.L)
        for itemp in range(1, self.L):
            self.beta[itemp] = self.gamma**(itemp - (self.L - 1))
            if self.N*self.beta[itemp] <= 1:
                self.step_Mu[itemp] = self.C_Mu
#                self.step_S[itemp] = self.C_S
#                self.step_W[itemp] = self.C_W
            else:
                self.step_Mu[itemp] = self.C_Mu/(self.N*self.beta[itemp]
                                                 )**self.d_Mu
#                self.step_S[itemp] = self.C_S/(self.N*self.beta[itemp]
#                                               )**self.d_S
#                self.step_W[itemp] = self.C_W/(self.N*self.beta[itemp]
#                                               )**self.d_W
        self.count_accept_Mu = np.zeros(self.L)
#        self.count_accept_S = np.zeros(self.L)
#        self.count_accept_W = np.zeros(self.L)
        self.count_exchange = np.zeros(self.L)
        self.Mu_ar = np.zeros([self.L, self.K, self.N_sampling])
#        self.S_ar = np.zeros([self.L, self.K, self.N_sampling])
#        self.W_ar = np.zeros([self.L, self.K, self.N_sampling])
        self.MSE = np.zeros([self.L, self.N_sampling]) + 9999
        for itemp in range(0, L):
            self.Mu_ar[itemp, :, 0] = np.random.uniform(self.min_Mu,
                                                        self.max_Mu, self.K)
#            self.S_ar[itemp, :, 0] = np.random.uniform(self.min_S, self.max_S,
#                                                       self.K)
#            self.W_ar[itemp, :, 0] = np.random.uniform(self.min_W, self.max_W,
#                                                       self.K)

#    def calc_E(self, _Mu, _S, _W):
    def calc_E(self, _Mu):
        # Here I use cryspy routines
        #yhat = np.zeros(self.data.shape[0])
        #for ik in range(0, self.K):
            #yhat += _W[ik]*np.exp(-(self.data[:, 0]-_Mu[ik])**2/(2*_S[ik]**2))
        #self.rhochi_dict['crystal_phase1']['unit_cell_parameters'][0] = _Mu[0]
        #_Mu[np.abs(_Mu) < 1E-10] = 1.E-10
        self.job.phases[0].cell.a = _Mu[0]
        self.job.phases[0].cell.b = _Mu[1]
        self.job.phases[0].cell.c = _Mu[2]
        #out = rhochi_calc_chi_sq_by_dictionary(self.rhochi_dict, dict_in_out={})
        #print(out[0], self.rhochi_dict['crystal_phase1']['unit_cell_parameters'][0])
        out = (np.square((self.job.create_simulation(self.meas_x) - self.meas_y)/self.meas_e)/self.meas_x.shape[0]).sum()
        #plt.plot(self.meas_x, self.meas_y)
        #plt.plot(self.meas_x, self.job.create_simulation(self.meas_x), label=str(_Mu[0])+" "+str(np.log(out)))
        #plt.plot(self.meas_x, self.meas_y - self.job.create_simulation(self.meas_x), label="diff")
        #plt.legend()
        #plt.show()
        #print(out, _Mu)
        return out

    def emcmc(self, isamp, itemp):
        self.Mu_ar[itemp, :, isamp] = self.Mu_ar[itemp, :, isamp-1]
#        self.S_ar[itemp, :, isamp] = self.S_ar[itemp, :, isamp-1]
#        self.W_ar[itemp, :, isamp] = self.W_ar[itemp, :, isamp-1]
        self.MSE[itemp, isamp] = self.MSE[itemp, isamp-1]
        pre_MSE = self.MSE[itemp, isamp]
        current_Mu = self.Mu_ar[itemp, :, isamp]
#        current_S = self.S_ar[itemp, :, isamp]
#        current_W = self.W_ar[itemp, :, isamp]
        NONNEG = False
        if itemp == 0:
            while not NONNEG:
                next_Mu = np.random.uniform(self.min_Mu, self.max_Mu, self.K)
                NONNEG = all(_Mu > 0. for _Mu in next_Mu)
        else:
            while not NONNEG:
                next_Mu = current_Mu + np.abs(np.random.normal(0.0, self.step_Mu[itemp],
                                                               self.K))
                NONNEG = all(_Mu > 0. for _Mu in next_Mu)
        #next_MSE = self.calc_E(next_Mu, current_S, current_W)
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
#        pre_MSE = self.MSE[itemp, isamp]
#        current_Mu = self.Mu_ar[itemp, :, isamp]
#        current_S = self.S_ar[itemp, :, isamp]
#        current_W = self.W_ar[itemp, :, isamp]
#        if itemp == 0:
#            next_S = np.random.uniform(self.min_S, self.max_S, self.K)
#        else:
#            next_S = current_S + np.random.normal(0.0, self.step_S[itemp],
#                                                  self.K)
#        next_MSE = self.calc_E(current_Mu, next_S, current_W)
#        r = np.exp(self.N*self.beta[itemp]/(2.*self.sigma**2)*(-next_MSE+pre_MSE))
#        for ik in range(0, self.K):
#            if next_S[ik] < self.min_S or next_S[ik] > self.max_S:
#                r = -1
#        if r >= np.random.uniform(0.0, 1.0) or itemp == 0 or isamp == 1:
#            self.S_ar[itemp, :, isamp] = next_S
#            self.MSE[itemp, isamp] = next_MSE
#            self.count_accept_S[itemp] += 1
#        pre_MSE = self.MSE[itemp, isamp]
#        current_Mu = self.Mu_ar[itemp, :, isamp]
#        current_S = self.S_ar[itemp, :, isamp]
#        current_W = self.W_ar[itemp, :, isamp]
#        if itemp == 0:
#            next_W = np.random.uniform(self.min_W, self.max_W, self.K)
#        else:
#            next_W = current_W + np.random.normal(0.0, self.step_W[itemp],
#                                                  self.K)
#        next_MSE = self.calc_E(current_Mu, next_S, current_W)
#        r = np.exp(self.N*self.beta[itemp]/(2.*self.sigma**2)*(-next_MSE+pre_MSE))
#        for ik in range(0, self.K):
#            if next_W[ik] < self.min_W or next_W[ik] > self.max_W:
#                r = -1
#        if r >= np.random.uniform(0.0, 1.0) or itemp == 0 or isamp == 1:
#            self.W_ar[itemp, :, isamp] = next_W
#            self.MSE[itemp, isamp] = next_MSE
#            self.count_accept_W[itemp] += 1
## I should undersand the actual spectrum data set to be fitted.
        if isamp == self.N_sampling - 1 and itemp == self.L-1:
##            yhat = np.zeros(self.data[:, 0].shape[0])
            minMSE_sample = np.where(self.MSE[itemp] ==
                                     np.min(self.MSE[itemp]))[0][0]
##            for ik in range(0, self.K):
##                yhat += self.W_ar[itemp, ik, minMSE_sample] *\
##                        np.exp(-(self.data[:, 0] -
##                                 self.Mu_ar[itemp, ik, minMSE_sample])**2
##                               / (2*self.S_ar[itemp, ik, minMSE_sample]**2))
            print("minMSE is found at the lattice constant of ",
                  self.Mu_ar[itemp, :, minMSE_sample],
                  self.calc_E(self.Mu_ar[itemp, :, minMSE_sample]))
            #self.job.phases[0].cell.a = self.Mu_ar[itemp, 0, minMSE_sample]
            #self.job.phases[0].cell.b = self.Mu_ar[itemp, 1, minMSE_sample]
            #self.job.phases[0].cell.c = self.Mu_ar[itemp, 2, minMSE_sample]
            #calc_y_cryspy = self.job.create_simulation(self.meas_x)
            #plt.plot(self.meas_x, calc_y_cryspy)
            #plt.plot(self.meas_x, self.meas_y)
            #plt.plot(self.meas_x, self.meas_y - calc_y_cryspy)
            #plt.show()
##            plt.plot(self.data[:, 0], yhat, ".")
##            plt.plot(self.data[:, 0], self.data[:, 1], ".")
##            plt.show()

    def rmc(self):
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        for isamp in range(1, self.N_sampling):
            #if rank == 0:
            #    print(isamp)
            if isamp % 5 == 0 and rank == 0:
                print(isamp)
            if isamp == self.burn_in-1:
                self.count_accept_Mu = 0*self.count_accept_Mu
#                self.count_accept_S = 0*self.count_accept_S
#                self.count_accept_W = 0*self.count_accept_W
                self.count_exchange = 0*self.count_exchange
            #for itemp in range(0, self.L):
            #    self.emcmc(isamp, itemp)
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
                #self.S_ar[istart:iend, :, isamp] = np.array(comm.allgather(self.S_ar[itemp, :, isamp])).reshape((size, self.K))
                #self.W_ar[istart:iend, :, isamp] = np.array(comm.allgather(self.W_ar[itemp, :, isamp])).reshape((size, self.K))
                self.count_accept_Mu[istart:iend] =\
                    np.array(comm.allgather(self.count_accept_Mu[itemp]))
                #self.count_accept_S[istart:iend] = np.array(comm.allgather(self.count_accept_S[itemp]))
                #self.count_accept_W[istart:iend] = np.array(comm.allgather(self.count_accept_W[itemp]))
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
#                        tmp_S = np.array(self.S_ar[itemp, :, isamp])
#                        self.S_ar[itemp, :, isamp] =\
#                            np.array(self.S_ar[itemp+1, :, isamp])
#                        self.S_ar[itemp+1, :, isamp] = tmp_S
#                        tmp_W = np.array(self.W_ar[itemp, :, isamp])
#                        self.W_ar[itemp, :, isamp] =\
#                            np.array(self.W_ar[itemp+1, :, isamp])
#                        self.W_ar[itemp+1, :, isamp] = tmp_W
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
#                        tmp_S = np.array(self.S_ar[itemp, :, isamp])
#                        self.S_ar[itemp, :, isamp] =\
#                            np.array(self.S_ar[itemp + 1, :, isamp])
#                        self.S_ar[itemp+1, :, isamp] = tmp_S
#                        tmp_W = np.array(self.W_ar[itemp, :, isamp])
#                        self.W_ar[itemp, :, isamp] =\
#                            np.array(self.W_ar[itemp+1, :, isamp])
#                        self.W_ar[itemp+1, :, isamp] = tmp_W
        if rank == 0:
            print("Accept rates for different beta")
            print(self.count_accept_Mu/(self.N_sampling - self.burn_in+1))
#        print(self.count_accept_S/(self.N_sampling - self.burn_in+1))
#        print(self.count_accept_W/(self.N_sampling - self.burn_in+1))
            print("Exchange rates for different beta")
            print(2*self.count_exchange/(self.N_sampling - self.burn_in+1))
            output = np.array([self.beta, self.count_accept_Mu/(self.N_sampling - self.burn_in+1),
#                              self.count_accept_S/(self.N_sampling - self.burn_in+1),
#                              self.count_accept_W/(self.N_sampling - self.burn_in+1),
                              2*self.count_exchange/(self.N_sampling - self.burn_in+1),
                              np.mean(self.N*self.MSE/self.sigma**2, axis=1)])
            np.savetxt("./results.txt", output.T)

    def init_easyDiffractionLib(self):
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        phases = Phases.from_cif_file(self.crystfile)
        self.meas_x, self.meas_y, self.meas_e = np.loadtxt(self.exptfile,
                                                           unpack=True)
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
            #calc_y_cryspy = self.job.create_simulation(self.meas_x)
            #plt.plot(self.meas_x, calc_y_cryspy)
            #plt.plot(self.meas_x, self.meas_y)
            #plt.plot(self.meas_x, self.meas_y - calc_y_cryspy)
            #plt.show()


def samplerun():
    crystfile = "/home/kazu/desktop/231207/easyDiffraction/examples/PbSO4.cif"
    exptfile = "/home/kazu/desktop/231207/easyDiffraction/examples/D1A@ILL.xye"
    prj = mcmc(crystfile, exptfile)
    prj.rmc()
    #prj.init_easyDiffractionLib()


samplerun()
