#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np


class mcmc():
    def __init__(self, datafile, K=3, N_sampling=20000, burn_in=10000, L=40,
                 gamma=1.3, d_Mu=1.0, C_Mu=1.0, d_S=1.0, C_S=0.6, d_W=0.9,
                 C_W=1.0, sigma=0.1, min_Mu=0., max_Mu=2.5, min_S=0.01,
                 max_S=1.5, min_W=0., max_W=2.0, isplot=True):
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
        self.d_W = d_W
        self.C_W = C_W
        self.sigma = sigma
        self.min_Mu = min_Mu
        self.max_Mu = max_Mu
        self.min_S = min_S
        self.max_S = max_S
        self.min_W = min_W
        self.max_W = max_W
        self.isplot = isplot
        self.data = np.loadtxt(self.datafile)
        plt.plot(self.data[:, 0], self.data[:, 1])
        self.N = self.data.shape[0]
        print('N:', self.N)
        np.random.seed(8)
        self.beta = np.zeros(L)
        self.beta[0] = 0
        self.step_Mu = np.zeros(self.L)
        self.step_S = np.zeros(self.L)
        self.step_W = np.zeros(self.L)
        for itemp in range(1, self.L):
            self.beta[itemp] = self.gamma**(itemp - (self.L - 1))
            if self.N*self.beta[itemp] <= 1:
                self.step_Mu[itemp] = self.C_Mu
                self.step_S[itemp] = self.C_S
                self.step_W[itemp] = self.C_W
            else:
                self.step_Mu[itemp] = self.C_Mu/(self.N*self.beta[itemp]
                                                 )**self.d_Mu
                self.step_S[itemp] = self.C_S/(self.N*self.beta[itemp]
                                               )**self.d_S
                self.step_W[itemp] = self.C_W/(self.N*self.beta[itemp]
                                               )**self.d_W
        self.count_accept_Mu = np.zeros(self.L)
        self.count_accept_S = np.zeros(self.L)
        self.count_accept_W = np.zeros(self.L)
        self.count_exchange = np.zeros(self.L)
        self.Mu_ar = np.zeros([self.L, self.K, self.N_sampling])
        self.S_ar = np.zeros([self.L, self.K, self.N_sampling])
        self.W_ar = np.zeros([self.L, self.K, self.N_sampling])
        self.MSE = np.zeros([self.L, self.N_sampling]) + 9999
        for itemp in range(0, L):
            self.Mu_ar[itemp, :, 0] = np.random.uniform(self.min_Mu,
                                                        self.max_Mu, self.K)
            self.S_ar[itemp, :, 0] = np.random.uniform(self.min_S, self.max_S,
                                                       self.K)
            self.W_ar[itemp, :, 0] = np.random.uniform(self.min_W, self.max_W,
                                                       self.K)

    def calc_E(self, _Mu, _S, _W):
        yhat = np.zeros(self.data.shape[0])
        for ik in range(0, self.K):
            yhat += _W[ik]*np.exp(-(self.data[:, 0]-_Mu[ik])**2/(2*_S[ik]**2))
        return np.mean((self.data[:, 1] - yhat)**2)

    def emcmc(self, isamp, itemp):
        self.Mu_ar[itemp, :, isamp] = self.Mu_ar[itemp, :, isamp-1]
        self.S_ar[itemp, :, isamp] = self.S_ar[itemp, :, isamp-1]
        self.W_ar[itemp, :, isamp] = self.W_ar[itemp, :, isamp-1]
        self.MSE[itemp, isamp] = self.MSE[itemp, isamp-1]
        pre_MSE = self.MSE[itemp, isamp]
        current_Mu = self.Mu_ar[itemp, :, isamp]
        current_S = self.S_ar[itemp, :, isamp]
        current_W = self.W_ar[itemp, :, isamp]
        if itemp == 0:
            next_Mu = np.random.uniform(self.min_Mu, self.max_Mu, self.K)
        else:
            next_Mu = current_Mu + np.random.normal(0.0, self.step_Mu[itemp],
                                                    self.K)
        next_MSE = self.calc_E(next_Mu, current_S, current_W)
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
        current_W = self.W_ar[itemp, :, isamp]
        if itemp == 0:
            next_S = np.random.uniform(self.min_S, self.max_S, self.K)
        else:
            next_S = current_S + np.random.normal(0.0, self.step_S[itemp],
                                                  self.K)
        next_MSE = self.calc_E(current_Mu, next_S, current_W)
        r = np.exp(self.N*self.beta[itemp]/(2.*self.sigma**2)*(-next_MSE+pre_MSE))
        for ik in range(0, self.K):
            if next_S[ik] < self.min_S or next_S[ik] > self.max_S:
                r = -1
        if r >= np.random.uniform(0.0, 1.0) or itemp == 0 or isamp == 1:
            self.S_ar[itemp, :, isamp] = next_S
            self.MSE[itemp, isamp] = next_MSE
            self.count_accept_S[itemp] += 1
        pre_MSE = self.MSE[itemp, isamp]
        current_Mu = self.Mu_ar[itemp, :, isamp]
        current_S = self.S_ar[itemp, :, isamp]
        current_W = self.W_ar[itemp, :, isamp]
        if itemp == 0:
            next_W = np.random.uniform(self.min_W, self.max_W, self.K)
        else:
            next_W = current_W + np.random.normal(0.0, self.step_W[itemp],
                                                  self.K)
        next_MSE = self.calc_E(current_Mu, next_S, current_W)
        r = np.exp(self.N*self.beta[itemp]/(2.*self.sigma**2)*(-next_MSE+pre_MSE))
        for ik in range(0, self.K):
            if next_W[ik] < self.min_W or next_W[ik] > self.max_W:
                r = -1
        if r >= np.random.uniform(0.0, 1.0) or itemp == 0 or isamp == 1:
            self.W_ar[itemp, :, isamp] = next_W
            self.MSE[itemp, isamp] = next_MSE
            self.count_accept_W[itemp] += 1
        if isamp == self.N_sampling - 1 and itemp == self.L-1 and self.isplot:
            yhat = np.zeros(self.data[:, 0].shape[0])
            minMSE_sample = np.where(self.MSE[itemp] ==
                                     np.min(self.MSE[itemp]))[0][0]
            for ik in range(0, self.K):
                yhat += self.W_ar[itemp, ik, minMSE_sample] *\
                        np.exp(-(self.data[:, 0] -
                                 self.Mu_ar[itemp, ik, minMSE_sample])**2
                               / (2*self.S_ar[itemp, ik, minMSE_sample]**2))
            plt.plot(self.data[:, 0], yhat, ".")
            plt.plot(self.data[:, 0], self.data[:, 1], ".")
            plt.show()

    def rmc(self):
        for isamp in range(1, self.N_sampling):
            if isamp % 1000 == 0:
                print(isamp)
            if isamp == self.burn_in-1:
                self.count_accept_Mu = 0*self.count_accept_Mu
                self.count_accept_S = 0*self.count_accept_S
                self.count_accept_W = 0*self.count_accept_W
                self.count_exchange = 0*self.count_exchange
            for itemp in range(0, self.L):
                self.emcmc(isamp, itemp)
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
                        tmp_W = np.array(self.W_ar[itemp, :, isamp])
                        self.W_ar[itemp, :, isamp] =\
                            np.array(self.W_ar[itemp+1, :, isamp])
                        self.W_ar[itemp+1, :, isamp] = tmp_W
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
                        tmp_W = np.array(self.W_ar[itemp, :, isamp])
                        self.W_ar[itemp, :, isamp] =\
                            np.array(self.W_ar[itemp+1, :, isamp])
                        self.W_ar[itemp+1, :, isamp] = tmp_W
        print(self.count_accept_Mu/(self.N_sampling - self.burn_in+1))
        print(self.count_accept_S/(self.N_sampling - self.burn_in+1))
        print(self.count_accept_W/(self.N_sampling - self.burn_in+1))
        print(2*self.count_exchange/(self.N_sampling - self.burn_in+1))
        output = np.array([self.beta, self.count_accept_Mu/(self.N_sampling - self.burn_in+1),
                          self.count_accept_S/(self.N_sampling - self.burn_in+1),
                          self.count_accept_W/(self.N_sampling - self.burn_in+1),
                          2*self.count_exchange/(self.N_sampling - self.burn_in+1),
                          np.mean(self.N*self.MSE/self.sigma**2, axis=1)])
        print(output.shape)
        np.savetxt("./results.txt", output.T)


def sample_run():
    datafile = "/home/kazu/desktop/230630/data.txt"
    prj = mcmc(datafile)
    prj.rmc()


#sample_run()
