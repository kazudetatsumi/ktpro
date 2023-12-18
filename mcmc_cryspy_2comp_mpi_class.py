#!/usr/bin/env python
#import matplotlib.pyplot as plt
import numpy as np
import cryspy
import warnings
#import os
from mpi4py import MPI
from cryspy.procedure_rhochi.rhochi_by_dictionary import \
    rhochi_lsq_by_dictionary, rhochi_rietveld_refinement_by_dictionary,\
    rhochi_calc_chi_sq_by_dictionary
warnings.filterwarnings('ignore')
#import sys
#import pprint

#pprint.pprint(sys.path)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

class mcmc():
# Test for optimization on the lattice constant
#    def __init__(self, datafile, K=3, N_sampling=20000, burn_in=10000, L=40,
#                 gamma=1.3, d_Mu=1.0, C_Mu=1.0, d_S=1.0, C_S=0.6, d_W=0.9,
#                 C_W=1.0, sigma=0.1, min_Mu=0., max_Mu=2.5, min_S=0.01,
#                 max_S=1.5, min_W=0., max_W=2.0):
    def __init__(self, datafile, K=1, N_sampling=8000, burn_in=4000, L=36,
                 gamma=1.3, d_Mu=1.0, C_Mu=1.0, d_S=1.0, C_S=1.0, sigma=0.1,
                 min_Mu=0., max_Mu=20., min_S=0., max_S=1.):
        self.datafile = datafile
        self.K = K
        self.N_sampling = N_sampling
        self.burn_in = burn_in
        self.L = L
        self.gamma = gamma
        self.d_Mu = d_Mu
        self.C_Mu = C_Mu
        self.d_S = d_S
        self.C_S = C_S
#        self.d_W = d_W
#        self.C_W = C_W
        self.sigma = sigma
        self.min_Mu = min_Mu
        self.max_Mu = max_Mu
        self.min_S = min_S
        self.max_S = max_S
#        self.min_W = min_W
#        self.max_W = max_W
        rhochi = cryspy.load_file(self.datafile)
        self.rhochi_dict = rhochi.get_dictionary()
        self.N = 800
        np.random.seed(8)
        self.beta = np.zeros(L)
        self.beta[0] = 0
        self.step_Mu = np.zeros(self.L)
        self.step_S = np.zeros(self.L)
#        self.step_W = np.zeros(self.L)
        for itemp in range(1, self.L):
            self.beta[itemp] = self.gamma**(itemp - (self.L - 1))
            if self.N*self.beta[itemp] <= 1:
                self.step_Mu[itemp] = self.C_Mu
                self.step_S[itemp] = self.C_S
#                self.step_W[itemp] = self.C_W
            else:
                self.step_Mu[itemp] = self.C_Mu/(self.N*self.beta[itemp]
                                                 )**self.d_Mu
                self.step_S[itemp] = self.C_S/(self.N*self.beta[itemp]
                                               )**self.d_S
#                self.step_W[itemp] = self.C_W/(self.N*self.beta[itemp]
#                                               )**self.d_W
        self.count_accept_Mu = np.zeros(self.L)
        self.count_accept_S = np.zeros(self.L)
#        self.count_accept_W = np.zeros(self.L)
        self.count_exchange = np.zeros(self.L)
        self.Mu_ar = np.zeros([self.L, self.K, self.N_sampling])
        self.S_ar = np.zeros([self.L, self.K, self.N_sampling])
#        self.W_ar = np.zeros([self.L, self.K, self.N_sampling])
        self.MSE = np.zeros([self.L, self.N_sampling]) + 9999
        for itemp in range(0, L):
            self.Mu_ar[itemp, :, 0] = np.random.uniform(self.min_Mu,
                                                        self.max_Mu, self.K)
            self.S_ar[itemp, :, 0] = np.random.uniform(self.min_S, self.max_S,
                                                       self.K)
#            self.W_ar[itemp, :, 0] = np.random.uniform(self.min_W, self.max_W,
#                                                       self.K)

#    def calc_E(self, _Mu, _S, _W):
    def calc_E(self, _Mu, _S):
        # Here I use cryspy routines
        #yhat = np.zeros(self.data.shape[0])
        #for ik in range(0, self.K):
            #yhat += _W[ik]*np.exp(-(self.data[:, 0]-_Mu[ik])**2/(2*_S[ik]**2))
        self.rhochi_dict['crystal_phase1']['unit_cell_parameters'][0] = _Mu[0]
        self.rhochi_dict['crystal_phase1']['unit_cell_parameters'][1] = _Mu[0]
        self.rhochi_dict['crystal_phase1']['unit_cell_parameters'][2] = _Mu[0]
        self.rhochi_dict['crystal_phase1']['atom_fract_xyz'][0][2] = _S[0]
        out = rhochi_calc_chi_sq_by_dictionary(self.rhochi_dict, dict_in_out={}
                                               )
        #print(out[0], self.rhochi_dict['crystal_phase1']['unit_cell_parameters'][0],
        #      self.rhochi_dict['crystal_phase1']['atom_fract_xyz'][0][2])
        return out[0]

    def emcmc(self, isamp, itemp):
        self.Mu_ar[itemp, :, isamp] = self.Mu_ar[itemp, :, isamp-1]
        self.S_ar[itemp, :, isamp] = self.S_ar[itemp, :, isamp-1]
#        self.W_ar[itemp, :, isamp] = self.W_ar[itemp, :, isamp-1]
        self.MSE[itemp, isamp] = self.MSE[itemp, isamp-1]
        pre_MSE = self.MSE[itemp, isamp]
        current_Mu = self.Mu_ar[itemp, :, isamp]
        current_S = self.S_ar[itemp, :, isamp]
#        current_W = self.W_ar[itemp, :, isamp]
        if itemp == 0:
            next_Mu = np.random.uniform(self.min_Mu, self.max_Mu, self.K)
        else:
            next_Mu = current_Mu + np.random.normal(0.0, self.step_Mu[itemp],
                                                    self.K)
        #next_MSE = self.calc_E(next_Mu, current_S, current_W)
        if next_Mu < 0.:
            print("WARNING, negative next_Mu", next_Mu)
        next_MSE = self.calc_E(next_Mu, current_S)
        if next_Mu < 0.:
            print("WARNING, next_MSE of negative next_Mu", next_MSE)
        r = np.exp(self.N * self.beta[itemp] / (2.*self.sigma**2) *
                   (-next_MSE+pre_MSE))
        for ik in range(0, self.K):
            if next_Mu[ik] < self.min_Mu or next_Mu[ik] > self.max_Mu:
                r = -1
        if r >= np.random.uniform(0.0, 1.0) or itemp == 0 or isamp == 1:
            self.Mu_ar[itemp, :, isamp] = next_Mu
            self.MSE[itemp, isamp] = next_MSE
            self.count_accept_Mu[itemp] += 1
        pre_MSE = self.MSE[itemp, isamp]
        current_Mu = self.Mu_ar[itemp, :, isamp]
        current_S = self.S_ar[itemp, :, isamp]
#        current_W = self.W_ar[itemp, :, isamp]
        if itemp == 0:
            next_S = np.random.uniform(self.min_S, self.max_S, self.K)
        else:
            next_S = current_S + np.random.normal(0.0, self.step_S[itemp],
                                                  self.K)
        next_MSE = self.calc_E(current_Mu, next_S)
        r = np.exp(self.N*self.beta[itemp]/(2.*self.sigma**2)*(-next_MSE+pre_MSE))
        for ik in range(0, self.K):
            if next_S[ik] < self.min_S or next_S[ik] > self.max_S:
                r = -1
        if r >= np.random.uniform(0.0, 1.0) or itemp == 0 or isamp == 1:
            self.S_ar[itemp, :, isamp] = next_S
            self.MSE[itemp, isamp] = next_MSE
            self.count_accept_S[itemp] += 1
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
            print("minMSE is found as", np.min(self.MSE[itemp]),
                  " at the lattice constant of ",
                  self.Mu_ar[itemp, 0, minMSE_sample],
                  self.S_ar[itemp, 0, minMSE_sample],
                  "at isamp ", minMSE_sample)
##            plt.plot(self.data[:, 0], yhat, ".")
##            plt.plot(self.data[:, 0], self.data[:, 1], ".")
##            plt.show()

    def rmc(self):
        for isamp in range(1, self.N_sampling):
            if isamp % 1000 == 0 and rank == 0:
                print(isamp)
            if isamp == self.burn_in-1:
                self.count_accept_Mu = 0*self.count_accept_Mu
                self.count_accept_S = 0*self.count_accept_S
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
                self.S_ar[istart:iend, :, isamp] =\
                    np.array(comm.allgather(self.S_ar[itemp, :, isamp])
                             ).reshape((size, self.K))
                #self.W_ar[istart:iend, :, isamp] = np.array(comm.allgather(self.W_ar[itemp, :, isamp])).reshape((size, self.K))
                self.count_accept_Mu[istart:iend] =\
                    np.array(comm.allgather(self.count_accept_Mu[itemp]))
                self.count_accept_S[istart:iend] =\
                    np.array(comm.allgather(self.count_accept_S[itemp]))
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
                        tmp_S = np.array(self.S_ar[itemp, :, isamp])
                        self.S_ar[itemp, :, isamp] =\
                            np.array(self.S_ar[itemp+1, :, isamp])
                        self.S_ar[itemp+1, :, isamp] = tmp_S
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
                        tmp_S = np.array(self.S_ar[itemp, :, isamp])
                        self.S_ar[itemp, :, isamp] =\
                            np.array(self.S_ar[itemp + 1, :, isamp])
                        self.S_ar[itemp+1, :, isamp] = tmp_S
#                        tmp_W = np.array(self.W_ar[itemp, :, isamp])
#                        self.W_ar[itemp, :, isamp] =\
#                            np.array(self.W_ar[itemp+1, :, isamp])
#                        self.W_ar[itemp+1, :, isamp] = tmp_W
        if rank == 0:
            print(self.count_accept_Mu/(self.N_sampling - self.burn_in+1))
            print(self.count_accept_S/(self.N_sampling - self.burn_in+1))
#        print(self.count_accept_W/(self.N_sampling - self.burn_in+1))
            print(2*self.count_exchange/(self.N_sampling - self.burn_in+1))
            output = np.array([self.beta, self.count_accept_Mu /
                              (self.N_sampling - self.burn_in+1),
                               self.count_accept_S /
                              (self.N_sampling - self.burn_in+1),
#                          self.count_accept_W/(self.N_sampling - self.burn_in+1),
                              2*self.count_exchange /
                              (self.N_sampling - self.burn_in+1),
                              np.mean(self.N*self.MSE/self.sigma**2, axis=1)])
            np.savetxt("./rsults.txt", output.T)

    def test_cryspy(self):
        rhochi = cryspy.load_file(self.datafile)
        rhochi_dict = rhochi.get_dictionary()
        print(rhochi_dict['crystal_phase1']['unit_cell_parameters'][0])
        dict_in_out = {}
        out = rhochi_calc_chi_sq_by_dictionary(rhochi_dict,
                                               dict_in_out=dict_in_out)
        print(out[0])
        rhochi_dict['crystal_phase1']['unit_cell_parameters'][0] = 1.59559
        print(rhochi_dict['crystal_phase1']['unit_cell_parameters'][0])
        out = rhochi_calc_chi_sq_by_dictionary(rhochi_dict,
                                               dict_in_out=dict_in_out)
        print(out[0])


def samplerun():
    datafile = "/home/kazu/desktop/231123/main.rcif"
    prj = mcmc(datafile)
    print(prj.test_cryspy())


#samplerun()
