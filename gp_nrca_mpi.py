#!/usr/bin/env python
# Test script for Gaussian Process Regression.
# x : set of positions [num_pos, vec] which can be a vector.
# f : Gaussian Process, scalar.
# y : f + noise (scalar), scalar.
# self.f_bar : mean in predictive distribution
# self.cov   : covariance matrix in predictive distribution
# Kazuyoshi TATSUMI 2022/05/02
import os
os.environ["OPENBLAS_NUM_THREADS"] = "36"
os.environ["MKL_NUM_THREADS"] = "36"
os.environ["VECLIB_NUM_THREADS"] = "36"
import numpy as np
import matplotlib.pyplot as plt
import pickle
from mpi4py import MPI


class GaussianProcessRegression:
    def __init__(self, x, x_train, y_train, noiselevel):
        self.x = x
        self.x_train = x_train
        self.y_train = y_train
        self.noiselevel = noiselevel

    def new(self):
        self.K = self.kernel(self.x_train, self.x_train)
        self.K_astarisc = self.kernel(self.x, self.x_train)
        self.noise_mat = self.noiselevel*np.eye(self.x_train.shape[0])
        self.K_I = self.K + self.noise_mat
        self.L = np.linalg.cholesky(self.K_I)
        print("calclulating f_bar")
        self.f_bar = self.K_astarisc @\
            np.linalg.lstsq(self.L.T, np.linalg.lstsq(self.L, self.y_train,
                                                      rcond=None)[0],
                            rcond=None)[0]
        print("calculating K_double_astarisc")
        self.K_double_astarisc = self.kernel_mpi(self.x, self.x)
        print("calculating V")
        self.V = np.linalg.lstsq(self.L, self.K_astarisc.T, rcond=None)[0]
        print("caluclating cov")
        self.conv = self.K_double_astarisc - self.V.T @ self.V
        print("calculating std")
        self.std = np.diag(self.conv)**0.5

    def renew(self, x_train, y_train):
        self.x_train = x_train
        self.y_train = y_train
        self.K = self.kernel(self.x_train, self.x_train)
        print("calculating K_astrisc")
        self.K_astarisc = self.kernel(self.x, self.x_train)
        print("calculating noise_mat")
        self.noise_mat = self.noiselevel*np.eye(self.x_train.shape[0])
        self.K_I = self.K + self.noise_mat
        self.L = np.linalg.cholesky(self.K_I)
        self.f_bar = self.K_astarisc @\
            np.linalg.lstsq(self.L.T, np.linalg.lstsq(self.L, self.y_train,
                                                      rcond=None)[0],
                            rcond=None)[0]
        self.V = np.linalg.lstsq(self.L, self.K_astarisc.T, rcond=None)[0]
        self.conv = self.K_double_astarisc - self.V.T @ self.V
        self.std = np.diag(self.conv)**0.5

    def kernel(self, x1, x2):
        K = np.zeros((x1.shape[0], x2.shape[0]))
        for p, xp in enumerate(x1):
            for q, xq in enumerate(x2):
                K[p, q] = np.exp(-0.5*np.sum(((xp - xq)/5.)**2))
        return(K)

    def kernel_d(self, x1, x2):
        rank = MPI.COMM_WORLD.Get_rank()
        K = np.zeros((x1.shape[0], x2.shape[0]))
        for p, xp in enumerate(x1):
            for q, xq in enumerate(x2):
                K[p, q] = np.exp(-0.5*np.sum(((xp - xq)/5.)**2))
        if rank == 0:
            print(K[3:8, 3], x1.shape, x2.shape)
        return(K)

    def kernel_mpi(self, x1, x2):
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
        size = MPI.COMM_WORLD.Get_size()
        #K = np.zeros((x1.shape[0], x2.shape[0]))
        nx1 = x1.shape[0] // size
        K = np.zeros((x1.shape[0]*x2.shape[0]))
        Ks = np.zeros((nx1, x2.shape[0]))
        #for p, xp in enumerate(x1):
        for p in range(rank*nx1, (rank+1)*nx1):
            for q, xq in enumerate(x2):
                Ks[p-rank*nx1, q] = np.exp(-0.5*np.sum(((x1[p] - xq)/5.)**2))
        Ks = Ks.flatten()
        comm.Allgather([Ks, MPI.DOUBLE], [K, MPI.DOUBLE])
        K = K.reshape((x1.shape[0], x2.shape[0]))
        if rank == 0:
            print(K[3:8, 3], x1.shape, x2.shape)
        return(K)

    def savedata(self, base=None):
        dataset = {}
        dataset['x_train'] = self.x_train
        dataset['y_train'] = self.y_train
        if base:
            dataset['base'] = base
        dataset['f_bar'] = self.f_bar
        dataset['std'] = self.std
        with open("./savedata.pkl", 'wb') as f:
            pickle.dump(dataset, f, 4)

    def loaddata(self, base=None):
        with open("./savedata.pkl", 'rb') as f:
            dataset = pickle.load(f)
            self.x_train = dataset['x_train']
            self.y_train = dataset['y_train']
            self.f_bar = dataset['f_bar']
            self.std = dataset['std']
            if base:
                return dataset['base']


def searchrun():
    x = np.linspace(-5., 5., 80)
    x_train = np.array([-2.3])
    y_train = testfunc(x_train)
    noiselevel = 0.000000000001
    prj = GaussianProcessRegression(x, x_train, y_train, noiselevel)
    prj.new()
    for itry in range(0, 20):
        plt.fill_between(x, prj.f_bar + 2.*prj.std, prj.f_bar - 2.*prj.std,
                         fc='gray')
        plt.scatter(x_train, y_train, marker='x')
        plt.plot(x, testfunc(x))
        plt.plot(x, prj.f_bar)
        plt.show()
        #beta = (0.1*np.log(itry*1.))**0.5
        beta = (20.*np.log(itry))**0.5
        print(itry)
        if itry % 3 == 0:
            activation = prj.f_bar + beta * prj.std
            nextx = x[np.argmax(activation)]
        elif itry % 3 == 1:
            activation = prj.f_bar - beta * prj.std
            nextx = x[np.argmin(activation)]
        else:
            activation = prj.std
            nextx = x[np.argmax(activation)]
        x_train = np.append(x_train, nextx)
        y_train = np.append(y_train, testfunc(nextx))
        prj.renew(x_train, y_train)


def searchrun2dnrca():
    prefix = "/home/kazu/ktpro/sample_data_nrca/"
    estimated_density_file = prefix + "estimated.dat"
    density_file = prefix + "dist_dens.dat"
    #estimated_err_file = prefix + "estimated_err.dat"
    density = []
    for iline, line in enumerate(open(density_file)):
        if iline != 0:
            density.append(float(line[:-1].split(',')[2]))
    density = np.array(density)

    X, Y = np.meshgrid(np.linspace(0., 5., 51), np.linspace(0., 5., 51))
    #X, Y = np.meshgrid(np.linspace(0., 5., 11), np.linspace(0., 5., 11))

    #density = density.reshape(X.shape)[::2, ::2].flatten()
    #X = X[::2, ::2]
    #Y = Y[::2, ::2]
    x = []
    for _x, _y in zip(X.flatten(), Y.flatten()):
        x.append([_x, _y])
    x = np.array(x)

    x_train = np.array([[2.5, 2.0]])
    y_train = np.array([testfunc2dnrca(X, Y, density, x_train)])
    noiselevel = 0.000000000001
    plt.pcolor(X, Y, density.reshape(X.shape), shading='auto', cmap='bwr', vmin=-0.0004, vmax=0.0004)
    plt.axis('equal')
    plt.show()

    for itry in range(0, 101):
        prj = GaussianProcessRegression(x, x_train, y_train, noiselevel)
        if itry % 10 == 0:
            print('itry:',itry)
            plt.subplot(1, 2, 1)
            #plt.pcolor(X, Y, prj.f_bar.reshape(X.shape), shading='auto', cmap='bwr', vmin=-0.2004, vmax=0.0004)
            plt.pcolor(X, Y, prj.f_bar.reshape(X.shape), shading='auto', cmap='bwr', vmin=-0.400, vmax=0.400)
            plt.axis('equal')
            plt.subplot(1, 2, 2)
            plt.pcolor(X, Y, prj.std.reshape(X.shape), shading='auto', cmap='Reds',vmin=0.)
            plt.scatter(prj.x_train[:,0], prj.x_train[:,1], marker='x')
            plt.axis('equal')
            plt.savefig('test_nrca_very_large_'+str(itry)+'.png')

        if itry == 0:
            beta = 0.
        else:
            beta = 1.*np.log(itry)**0.5
        for ixtr, xtrain in enumerate(x_train):
            mask = np.sum((x - xtrain)**2., axis=1) > 0.00001
            if ixtr == 0:
                tmask = mask
            tmask = mask*tmask
        if itry % 2 == 0:
            activation = prj.f_bar + beta * prj.std
            nextx = x[activation == np.max(activation[tmask])][0].squeeze()
        else:
            activation = prj.f_bar - beta * prj.std
            nextx = x[activation == np.min(activation[tmask])][0].squeeze()
        #nextxindx = np.unravel_index(activation.argmax(), X.shape)
        #print('next position:', X[nextxindx], Y[nextxindx])

        print('adding to x_train', nextx)
        x_train = np.vstack((x_train, nextx))
        y_train = np.append(y_train, testfunc2dnrca(X, Y, density, nextx[np.newaxis, :]))



def testfunc2dnrca(X, Y, density, xvec):
    test = (X.flatten() - xvec[0, 0])**2 + (Y.flatten() - xvec[0, 1])**2
    return density[test.argmin()]


def testfunc2d(x):
    pos2 = np.array([2., 3.])
    pos1 = np.array([1., 1.])
    return np.exp(-np.sum((x-pos1)**2., axis=1)) + 0.5*np.exp(-np.sum((x-pos2)**2., axis=1))


def testfunc(x):
    return 2.*np.exp(-x**2.) + np.exp(-((x-1.5)/2.0)**2.)


#searchrun2dnrca()
#searchrun()