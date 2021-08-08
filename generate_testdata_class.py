#!/usr/bin/env python
import numpy as np
from scipy.optimize import minimize


class generate_testdata:
    def __init__(self, means, ints, sigmas, ntot):
        self.means = means
        self.ints = ints
        self.sigmas = sigmas
        self.ntot = ntot

    def pdf(self):
        self.pdf = np.zeros((2, 100))
        x = np.linspace(0, 1, 100)
        for sigma, mean, weight in zip(self.sigmas, self.means, self.ints):
            self.pdf[1, :] += weight*np.exp(-(x - mean)**2/sigma**2)\
                         / (2.0*np.pi)**0.5 / sigma
        self.pdf[1, :] = self.pdf[1, :]/np.sum(self.pdf[1, :])/(x[1]-x[0])
        self.pdf[0, :] = x

    def rate(self, x):
        value = np.zeros_like(x)
        for sigma, mean, weight in zip(self.sigmas, self.means, self.ints):
            value += weight*np.exp(-(x - mean)**2/sigma**2)\
                         / (2.0*np.pi)**0.5 / sigma
        return value

    def neg_rate(self, x):
        return -self.rate(x)

    def get_maxrate(self):
        x0 = 0.5
        results = minimize(self.neg_rate, x0)
        print(results)
        self.maxrate = -results.fun

    def generate_homogeneous_point_process(self):
        self.numpoints = np.random.poisson(self.maxrate*self.ntot)
        self.homoxdata = np.random.uniform(0, 1, self.numpoints)

    def thinning(self):
        self.inhomoxdata = self.homoxdata[np.random.uniform(0, 1,
                                          self.numpoints) <
                                          self.rate(self.homoxdata)]

    def mise(self, x):
        bins = np.arange(0, 1.0+x, x)
        histdata = np.histogram(self.inhomoxdata, bins=bins, density=True)
        theta = histdata[0]
        x = histdata[1][:]
        mise = 0.0
        dpdfx = self.pdf[0, 1] - self.pdf[0, 0]
        for idx, xi in enumerate(x[:-1]):
            for jdx, pdf in enumerate(self.pdf[1, :]):
                if self.pdf[0, jdx] >= xi and self.pdf[0, jdx] < x[idx+1]:
                    mise += (theta[idx] - pdf)**2*dpdfx
        return(mise)

    def pdfhist(self, x):
        bins = np.arange(0, 1.0+x, x)
        dpdfx = self.pdf[0, 1] - self.pdf[0, 0]
        hist = np.zeros((2, bins.shape[0]-1))
        for idx, xi in enumerate(bins[0:-1]):
            hist[0, idx] = (bins[idx]+bins[idx+1])/2
            for jdx, pdf in enumerate(self.pdf[0, :]):
                if pdf >= xi and pdf < bins[idx+1]:
                    hist[1, idx] += self.pdf[1, jdx]*dpdfx/x
        return hist

    def fun_gauss(self, x, h):
        return np.exp(-x**2/h**2)/(2.0*np.pi)**0.5 / h

    def kde(self, h):
        kde = np.zeros_like(self.pdf)
        kde[0, :] = self.pdf[0, :]
        dx = kde[0, 1] - kde[0, 0]
        for _x in self.inhomoxdata:
            for ixx, _xx in enumerate(kde[0, :]):
                kde[1, ixx] += self.fun_gauss(_x-_xx, h)
        kde[1, :] = kde[1, :]/np.sum(kde[1, :])/dx
        return kde

    def misekde(self, h):
        kde = self.kde(h)
        misekde = 0.0
        dpdfx = self.pdf[0, 1] - self.pdf[0, 0]
        for jdx, pdf in enumerate(self.pdf[1, :]):
            misekde += (kde[1, jdx] - pdf)**2*dpdfx
        return(misekde)
        

#def samplerun():
#    proj = generate_testdata([0.3, 0.8], [0.2, 0.4], [0.2, 0.2], 400)
#    proj.get_maxrate()
#    print(proj.maxrate)
#    proj.generate_homogeneous_point_process()
#    proj.thinning()
#    proj.pdf()
#    print(proj.inhomoxdata)
#    plt.scatter(proj.inhomoxdata, np.zeros_like(proj.inhomoxdata), marker='|')
##    plt.plot(proj.pdf[0, :], proj.pdf[1, :])
#    plt.show()

#samplerun()
