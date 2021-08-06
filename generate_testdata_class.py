#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
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
        self.maxrate = -results.fun

    def generate_homogeneous_point_process(self):
        self.numpoints = np.random.poisson(self.maxrate*self.ntot)
        self.homoxdata = np.random.uniform(0, 1, self.numpoints)

    def thinning(self):
        self.inhomoxdata = self.homoxdata[np.random.uniform(0, 1,
                                          self.numpoints) <
                                          self.rate(self.homoxdata)]

    def histogram(self, nbins):


def samplerun():
    proj = generate_testdata([0.3, 0.8], [0.2, 0.4], [0.2, 0.2], 400)
    proj.get_maxrate()
    print(proj.maxrate)
    proj.generate_homogeneous_point_process()
    proj.thinning()
    proj.pdf()
    print(proj.inhomoxdata)
    plt.scatter(proj.inhomoxdata, np.zeros_like(proj.inhomoxdata), marker='|')
    plt.plot(proj.pdf[0, :], proj.pdf[1, :])
    plt.show()


samplerun()
