#!/usr/bin/env python
# Test script for Gaussian Process Regression.
# x : set of positions [num_pos, vec] which can be a vector.
# f : Gaussian Process, scalar.
# y : f + noise (scalar), scalar.
# self.f_bar : mean in predictive distribution
# self.cov   : covariance matrix in predictive distribution
# Kazuyoshi TATSUMI 2022/05/02
import numpy as np
import matplotlib.pyplot as plt


class GaussianProcessRegression:
    def __init__(self, x, x_train, y_train, noiselevel):
        self.x = x
        self.x_train = x_train
        self.y_train = y_train
        self.noiselevel = noiselevel
        self.K = self.kernel(self.x_train, self.x_train)
        self.K_astarisc = self.kernel(self.x, self.x_train)
        self.noise_mat = self.noiselevel*np.eye(self.x_train.shape[0])
        self.K_I = self.K + self.noise_mat
        self.L = np.linalg.cholesky(self.K_I)
        self.f_bar = self.K_astarisc @\
            np.linalg.lstsq(self.L.T, np.linalg.lstsq(self.L, self.y_train,
                                                      rcond=None)[0],
                            rcond=None)[0]
        self.K_double_astarisc = self.kernel(self.x, self.x)
        self.V = np.linalg.lstsq(self.L, self.K_astarisc.T, rcond=None)[0]
        self.conv = self.K_double_astarisc - self.V.T @ self.V
        self.std = np.diag(self.conv)**0.5

    def kernel(self, x1, x2):
        K = np.zeros((x1.shape[0], x2.shape[0]))
        for p, xp in enumerate(x1):
            for q, xq in enumerate(x2):
                K[p, q] = np.exp(-0.5*np.sum(((xp - xq)/1.)**2))
        return(K)


def run():
    x = np.linspace(-5., 5., 80)
    # x_train = np.array([-2.3, 1.0, 3.5, -1.0, -4.0])
    # y_train = np.array([1.11, 3.00,  -2.00, 4.0, 1.0])
    x_train = np.array([-2.3])
    y_train = np.array([1.11])
    noiselevel = 0.000000000001
    prj = GaussianProcessRegression(x, x_train, y_train, noiselevel)
    plt.fill_between(x, prj.f_bar + 2.*prj.std, prj.f_bar - 2.*prj.std,
                     fc='lightgray')
    plt.scatter(x_train, y_train)
    plt.plot(x, prj.f_bar)
    plt.plot(x, prj.std+prj.f_bar)
    plt.show()


def searchrun():
    x = np.linspace(-5., 5., 80)
    x_train = np.array([-2.3])
    y_train = testfunc(x_train)
    noiselevel = 0.000000000001
    for itry in range(0, 20):
        prj = GaussianProcessRegression(x, x_train, y_train, noiselevel)
        plt.fill_between(x, prj.f_bar + 2.*prj.std, prj.f_bar - 2.*prj.std,
                         fc='gray')
        plt.scatter(x_train, y_train, marker='x')
        plt.plot(x, testfunc(x))
        plt.plot(x, prj.f_bar)
        plt.show()
        beta = (0.1*np.log(itry*1.))**0.5
        beta =5.*np.log(itry)**0.5
        if itry % 2 == 1:
            activation = prj.f_bar + beta * prj.std
            nextx = x[np.argmax(activation)]
        else:
            activation = prj.f_bar - beta * prj.std
            nextx = x[np.argmin(activation)]
        #activation = prj.f_bar + beta * prj.std
        #nextx = x[np.argmax(activation)]
        x_train = np.append(x_train, nextx)
        y_train = np.append(y_train, testfunc(nextx))


def testfunc(x):
    return 2.*np.exp(-x**2.) + np.exp(-((x-1.5)/1.3)**2.)


searchrun()
