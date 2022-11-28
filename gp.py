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


    def kernelorg(self, x1, x2):
        K = np.zeros((x1.shape[0], x2.shape[0]))
        for p, xp in enumerate(x1):
            for q, xq in enumerate(x2):
                K[p, q] = np.exp(-0.5*np.sum(((xp - xq)/5.)**2))
        return(K)

    def kernel(self, x1, x2):
        K = np.zeros((x1.shape[0], x2.shape[0]))
        test_x1 = np.repeat(np.expand_dims(x1, 1), x2.shape[0], axis=1)
        test_x2 = np.repeat(np.expand_dims(x2, 1), x1.shape[0], axis=1
                            ).transpose(1, 0)
        K = np.exp(-0.5*((test_x1 - test_x2)/0.05)**2)
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


def searchrun2d():
    X, Y = np.meshgrid(np.linspace(-5., 5., 50), np.linspace(-5., 5., 50))
    x = []
    for _x, _y in zip(X.flatten(), Y.flatten()):
        x.append([_x, _y])
    x = np.array(x)

    x_train = np.array([[-2.3, -2.3]])
    y_train = testfunc2d(x_train)
    noiselevel = 0.000000000001
    plt.pcolor(X, Y, testfunc2d(x).reshape(X.shape), shading='auto', cmap='bwr', vmin=-1., vmax=1.)
    #plt.scatter(x_train[0], x_train[1], marker='x')
    plt.axis('equal')
    plt.show()

    for itry in range(0, 101):
        prj = GaussianProcessRegression(x, x_train, y_train, noiselevel)
        if itry % 10 == 0:
            print('itry:',itry)
            plt.subplot(1, 2, 1)
            plt.pcolor(X, Y, prj.f_bar.reshape(X.shape), shading='auto', cmap='bwr', vmin=-1., vmax=1.)
            plt.axis('equal')
            plt.subplot(1, 2, 2)
            plt.pcolor(X, Y, prj.std.reshape(X.shape), shading='auto', cmap='bwr', vmin=-1., vmax=1.)
            plt.axis('equal')
            plt.show()

        if itry == 0:
            beta = 0.
        else:
            beta = .2*np.log(itry)**0.5
        activation = prj.f_bar + beta * prj.std
        for ixtr, xtrain in enumerate(x_train):
            mask = np.sum((x - xtrain)**2., axis=1) > 0.00001
            if ixtr == 0:
                tmask = mask
            tmask = mask*tmask
        nextx = x[activation == np.max(activation[tmask])][0].squeeze()
        #nextxindx = np.unravel_index(activation.argmax(), X.shape)
        #print('next position:', X[nextxindx], Y[nextxindx])

        print('adding to x_train', nextx)
        x_train = np.vstack((x_train, nextx))
        y_train = np.append(y_train, testfunc2d(nextx[np.newaxis, :]))



def testfunc(x):
    return 2.*np.exp(-x**2.) + np.exp(-((x-1.5)/2.0)**2.)


def testfunc2d(x):
    pos2 = np.array([1., 2.])
    pos1 = np.array([0., 0.])
    return np.exp(-np.sum((x-pos1)**2., axis=1)) + 0.5*np.exp(-np.sum((x-pos2)**2., axis=1))

#searchrun2d()
