#!/usr/bin/env python
import pylab as plt
import numpy as np
from numba import jit


class mcmc():
    def __init__(self, datafile, K=2, N_sampling=20000, burn_in=10000, L=40,
                 gamma=1.3, d_Mu=1.0, C_Mu=1.0, d_S=1.0, C_S=0.6, d_W=0.9,
                 C_W=1.0, sigma=0.1, min_Mu=0., max_Mu=2.5, min_S=0.01,
                 max_S=1.5, min_W=0., max_W=2.0):
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
        self.data = np.loadtxt(self.datafile)
        plt.plot(self.data[:, 0], self.data[:, 1])
        self.N = self.data.shape[0]
        np.random.seed(8)
        self.beta = np.zeros(L)
        self.beta[0] = 0
        for itemp in range(1, self.L):
            self.beta[itemp] = self.gamma**(itemp - (self.L - 1))
            # beta = 1/(kB * T)
            # For large temperatures, step sizes are set so as to avoid from
            # diverged.
            if self.N*self.beta[itemp] <= 1:
                self.step_Mu[itemp] = self.C_Mu
                self.step_S[itemp] = self.C_S
                self.step_W[itemp] = self.C_W
            # For others, step sizes(=sigmas of normal distributions) are set
            # proportional to temperatures**d
            else:
                self.step_Mu[itemp] = self.C_Mu/(self.N*self.beta[itemp]
                                                 )**self.d_Mu
                self.step_S[itemp] = self.C_S/(self.N*self.beta[itemp]
                                               )**self.d_S
                self.step_W[itemp] = self.C_W/(self.N*self.beta[itemp]
                                               )**self.d_W
        # initialization of Mu, S, W
        #self.Mu = np.random.uniform(self.min_Mu, self.max_Mu, self.K)
        #self.S = np.random.uniform(self.min_S, self.max_S, self.K)
        #self.W = np.random.uniform(self.min_W, self.max_W, self.K)
        # Arrays for saving sampling results.
        # Numbers of accepted Mu, S,W for different temperatures.
        self.count_accept_Mu = np.zeros(self.L)
        self.count_accept_S = np.zeros(self.L)
        self.count_accept_W = np.zeros(self.L)
        # Numbers of exchanges performed for different temperatures.
        self.count_exchange = np.zeros(self.L)
        # Recording of Mu, S, W values for different temperatures & samplings.
        self.Mu_ar = np.zeros([self.L, self.K, self._sampling])
        self.S_ar = np.zeros([self.L, self.K, self.N_sampling])
        self.W_ar = np.zeros([self.L, self.K, self.N_sampling])
        # Recording of Mean Squared Errors for different temps and sammplings.
        self.MSE = np.zeros([self.L, self.N_sampling]) + 9999
        #self.Mu_ar[:, :, 0] = np.tile(self.Mu, self.L).reshape((self.L, self.K))
        #self.S_ar[:, :, 0] = np.tile(self.S, self.L).reshape((self.L, self.K))
        #self.W_ar[:, :, 0] = np.tile(self.W, self.L).reshape((self.L, self.K))
        for itemp in range(0, L):
            self.Mu_ar[itemp, :, 0] = np.random.uniform(self.min_Mu,
                                                        self.max_Mu, self.K)
            self.S_ar[itemp, :, 0] = np.random.uniform(self.min_S, self.max_S,
                                                       self.K)
            self.W_ar[itemp, :, 0] = np.random.uniform(self.min_W, self.max_W,
                                                       self.K)

    def calc_E(self, _Mu, _S, _W):
        yhat = np.zeros(self.data.shape[0])
        #datav = np.tile(self.data[:, 0], self.K).rehsape((self.K, -1)).T
        #yhat = _W*np.exp(-(datav - _Mu)**2 / (2.*_S**2))
        for ik in range(0, self.K):
            yhat += _W[ik]*np.exp(-(self.data[:, 0]-_Mu[ik])**2/(2*_S[ik]**2))
        return np.mean((self.data[:, 1] - yhat)**2)
        #return np.mean((self.data[:, 1] - np.sum(yhat, axis=1))**2)

    def emcmc(self, isampling, itemp):
        next_Mu = np.zeros(self.K)
        next_S = np.zeros(self.K)
        next_W = np.zeros(self.K)
        self.Mu_ar[itemp, :, isampling] = self.Mu_ar[itemp, :, isampling-1]
        self.S_ar[itemp, :, isampling] = self.S_ar[itemp, :, isampling-1]
        self.W_ar[itemp, :, isampling] = self.W_ar[itemp, :, isampling-1]
        self.MSE[itemp, isampling] = self.MSE[itemp, isampling-1]
        # update Mu
        pre_MSE = self.MSE[itemp, isampling]
        current_Mu = self.Mu_ar[itemp, :, isampling]
        current_S = self.S_ar[itemp, :, isampling]
        current_W = self.W_ar[itemp, :, isampling]

        if itemp == 0:
            next_Mu = np.random.uniform(self.min_Mu, self.max_Mu, self.K)
        else:
            next_Mu = current_Mu + np.random.normal(0.0, self.step_Mu[itemp],
                                                    self.K)                     # np.normal(mean, std, size)

        next_MSE = self.calc_E(next_Mu, current_S, current_W)                    # get MSE after only Mu is changed.
        r = np.exp(self.N*self.beta[itemp]/(2*self.sigma**2)*(-next_MSE+pre_MSE))           # r = exp(N/(kbT)*Delta(MSE)/(2sigma**2))
                                                                              # This is a Metropolis filter.
                                                                              # The target probability distribution for sampling
                                                                              # is p(theta|K, data) proportional to exp(-N*beta*MSE) with
                                                                              # p(K) and P(theta|K) assumed as uniform distributions.

        # prior(uniform distribution)
        for count_k in range(0, self.K):
            if next_Mu[count_k] < self.min_Mu or next_Mu[count_k] > self.max_Mu:            # If any of the changed Mu is out of the predefined Mu range,
                r = -1                                                            # the changed Mu are totally  rejected.

        if r >= np.random.uniform(0.0, 1.0) or count_temp == 0 or\
                isampling == 1:                                              # count_sampling is iterated from 1 to N_sampling - 1.
            self.Mu_ar[itemp, :, isampling] = next_Mu                        # Monte Carlo for accepting the changed Mu.
            self.MSE[itemp, isampling] = next_MSE
            count_accept_Mu[itemp] += 1

        # update S                                                                # Being similar to the case of updating Mu.
        pre_MSE = self.MSE[itemp, isampling]
        current_Mu = self.Mu_ar[itemp, :, isampling]
        current_S = self.S_ar[itemp, :, isampling]
        current_W = self.W_ar[itemp, :, isampling]

        if itemp == 0:
            next_S = np.random.uniform(self.min_S, self.max_S, self.K)
        else:
            next_S = current_S + np.random.normal(0.0, self.step_S[count_temp],
                                                  self.K)

        next_MSE = self.calc_E(current_Mu, next_S, current_W)
        r = np.exp(self.N*self.beta[itemp]/(2.*self.sigma**2)*(-next_MSE+pre_MSE))
        # prior(uniform distribution)
        for ik in range(0, self.K):
            if next_S[ik] < self.min_S or next_S[ik] > self.max_S:
                r = -1

        if r >= np.random.uniform(0.0, 1.0) or itemp == 0 or\
                isampling == 1:
            self.S_ar[itemp, :, isampling] = next_S
            self.MSE[itemp, isampling] = next_MSE
            self.count_accept_S[itemp] += 1

        # update W                                                                 # Being similar to the case of updating Mu or S.
        pre_MSE = self.MSE[itemp, isampling]
        current_Mu = self.Mu_ar[itemp, :, isampling]
        current_S = self.S_ar[itemp, :, isampling]
        current_W = self.W_ar[itemp, :, isampling]

        if itemp == 0:
            next_W = np.random.uniform(self.min_W, self.max_W, self.K)
        else:
            next_W = current_W + np.random.normal(0.0, self.step_W[itemp],
                                                  self.K)

        next_MSE = self.calc_E(current_Mu, next_S, current_W)
        r = np.exp(self.N*self.beta[itemp]/(2.*self.sigma**2)*(-next_MSE+pre_MSE))
        # prior(uniform distribution)
        for ik in range(0, self.K):
            if next_W[ik] < self.min_W or next_W[ik] > self.max_W:
                r = -1

        if r >= np.random.uniform(0.0, 1.0) or count_temp == 0 or\
                isampling == 1:
            self.W_ar[itemp, :, isampling] = next_W
            self.MSE[itemp, isampling] = next_MSE
            self.count_accept_W[itemp] += 1

        if isampling == self.N_sampling - 1 and itemp == self.L-1:
            yhat = np.zeros(self.data[:, 0].shape[0])
            minMSE_sample = np.where(self.MSE[itemp] ==
                                     np.min(self.MSE[itemp]))[0][0]   # selecting the youngest sample index whose MSE has a minimum value
                                                                  # for the final temperature.
            for ik in range(0, self.K):
                yhat += self.W_ar[itemp, ik, minMSE_sample] *\
                        np.exp(-(self.data[:, 0] -
                                 self.Mu_ar[count_temp, count_k, minMSE_sample])**2
                               / (2*self.S_ar[count_temp, count_k, minMSE_sample]**2))
            plt.plot(self.data[:, 0], yhat, ".")                           # Visually comparing the model at the minimum MSE at the final temperature
            plt.plot(self.data[:, 0], self.data[:, 1], ".")                     # with the raw spectrum.
            plt.show()

        return 0

############################
# Replica Monte Calro Method
############################
#@jit
# def main():


for count_iter in range(1, N_sampling):                                              # loop w.r.t. samplings is start here.
    if count_iter % 1000 == 0:
        print(count_iter)

    if count_iter == burn_in-1:
        count_accept_Mu = 0*count_accept_Mu
        count_accept_S = 0*count_accept_S
        count_accept_W = 0*count_accept_W
        count_exchange = 0*count_exchange

    for count_temp in range(0, L):                                                   # loop w.r.t. temperatures is start here
        emcmc(data, count_iter, count_temp, count_accept_Mu, count_accept_S,
              count_accept_W, Mu_ar, S_ar, W_ar, MSE)

    for count_temp in range(0, L-1):
        r_exchange = np.exp(
                            N/(2*sigma**2) *                                         # r_exchange = exp(N*[1/(kbT_1) -1/(kbT_2)]*Delta(MSE)/(2simga**2))
                            (beta[count_temp+1] - beta[count_temp]) *                #            = p_beta_1(theta_2|K, data)*p_beta_2(theta_1|K, data)/
                            (MSE[count_temp+1, count_iter] -                         #              [p_beta_1(theta_1|K, data)*p_beta_2(theta_2|K, data)]
                             MSE[count_temp, count_iter])
                            )
        # print(r_exchange)
        if count_iter % 2 == 0 and count_temp % 2 == 0:                             # results among different temperatures are exchanged
            if r_exchange >= np.random.uniform(0.0, 1.0):                           # for the current sampling index.
                count_exchange[count_temp] += 1

                tmp = MSE[count_temp, count_iter]
                MSE[count_temp, count_iter] = MSE[count_temp+1, count_iter]
                MSE[count_temp+1, count_iter] = tmp

                tmp_Mu = np.array(Mu_ar[count_temp, :, count_iter])
                Mu_ar[count_temp, :, count_iter] =\
                    np.array(Mu_ar[count_temp+1, :, count_iter])
                Mu_ar[count_temp+1, :, count_iter] = tmp_Mu

                tmp_S = np.array(S_ar[count_temp, :, count_iter])
                S_ar[count_temp, :, count_iter] =\
                    np.array(S_ar[count_temp+1, :, count_iter])
                S_ar[count_temp+1, :, count_iter] = tmp_S

                tmp_W = np.array(W_ar[count_temp, :, count_iter])
                W_ar[count_temp, :, count_iter] =\
                    np.array(W_ar[count_temp+1, :, count_iter])
                W_ar[count_temp+1, :, count_iter] = tmp_W

        elif count_iter % 2 == 1 and count_temp % 2 == 1:
            if r_exchange >= np.random.uniform(0.0, 1.0):
                count_exchange[count_temp] += 1

                tmp = MSE[count_temp, count_iter]
                MSE[count_temp, count_iter] = MSE[count_temp+1, count_iter]
                MSE[count_temp+1, count_iter] = tmp

                tmp_Mu = np.array(Mu_ar[count_temp, :, count_iter])
                Mu_ar[count_temp, :, count_iter] =\
                    np.array(Mu_ar[count_temp + 1, :, count_iter])
                Mu_ar[count_temp+1, :, count_iter] = tmp_Mu

                tmp_S = np.array(S_ar[count_temp, :, count_iter])
                S_ar[count_temp, :, count_iter] =\
                    np.array(S_ar[count_temp + 1, :, count_iter])
                S_ar[count_temp+1, :, count_iter] = tmp_S

                tmp_W = np.array(W_ar[count_temp, :, count_iter])
                W_ar[count_temp, :, count_iter] =\
                    np.array(W_ar[count_temp+1, :, count_iter])
                W_ar[count_temp+1, :, count_iter] = tmp_W

print(count_accept_Mu/(N_sampling - burn_in+1))
print(count_accept_S/(N_sampling - burn_in+1))
print(count_accept_W/(N_sampling - burn_in+1))
print(2*count_exchange/(N_sampling - burn_in+1))

# saving data
output = np.array([beta, count_accept_Mu/(N_sampling - burn_in+1),
                  count_accept_S/(N_sampling - burn_in+1),
                  count_accept_W/(N_sampling - burn_in+1),
                  2*count_exchange/(N_sampling - burn_in+1),
                  np.mean(N*MSE/sigma**2, axis=1)])
np.savetxt("./rsults.txt", output.T)
