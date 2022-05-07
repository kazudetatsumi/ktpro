#!/usr/bin/env python
# Test script for Gaussian Process Regression.
# x : set of positions [num_pos, vec] which can be a vector.
# f : Gaussian Process, scalar.
# y : f + noise (scalar), scalar.
# self.f_bar : mean in predictive distribution
# self.cov   : covariance matrix in predictive distribution
# Kazuyoshi TATSUMI 2022/05/02
import os
# os.environ["OPENBLAS_NUM_THREADS"] = "2"
# os.environ["MKL_NUM_THREADS"] = "2"
# os.environ["VECLIB_NUM_THREADS"] = "2"
import numpy as np
import matplotlib.pyplot as plt
import sys
home = os.path.expanduser("~")
sys.path.append(home+"/ktpro")
from gp_nrca import GaussanProcessRegression as gpr


def searchrun2dnrca():
    prefix = home + "/ktpro/sample_data_nrca/"
    # estimated_density_file = prefix + "estimated.dat"
    density_file = prefix + "dist_dens.dat"
    # estimated_err_file = prefix + "estimated_err.dat"
    density = []
    for iline, line in enumerate(open(density_file)):
        if iline != 0:
            density.append(float(line[:-1].split(',')[2]))
    density = np.array(density)

    X, Y = np.meshgrid(np.linspace(0., 5., 51), np.linspace(0., 5., 51))

    # density = density.reshape(X.shape)[::2, ::2].flatten()
    # X = X[::2, ::2]
    # Y = Y[::2, ::2]
    x = []
    for _x, _y in zip(X.flatten(), Y.flatten()):
        x.append([_x, _y])
    x = np.array(x)

    x_train = np.array([[2.5, 2.0]])
    y_train = np.array([testfunc2dnrca(X, Y, density, x_train)])
    noiselevel = 0.000000000001
    plt.pcolor(X, Y, density.reshape(X.shape), shading='auto', cmap='bwr',
               vmin=-0.4, vmax=0.4)
    plt.axis('equal')
    plt.savefig('test_nrca_very_large_org.png')
    # plt.show()

    for itry in range(0, 101):
        prj = gpr(x, x_train, y_train, noiselevel)
        if itry % 10 == 0:
            print('itry:', itry)
            plt.subplot(1, 2, 1)
            # plt.pcolor(X, Y, prj.f_bar.reshape(X.shape), shading='auto',
            #            cmap='bwr', vmin=-0.2004, vmax=0.0004)
            plt.pcolor(X, Y, prj.f_bar.reshape(X.shape), shading='auto',
                       cmap='bwr', vmin=-0.400, vmax=0.400)
            plt.axis('equal')
            plt.subplot(1, 2, 2)
            plt.pcolor(X, Y, prj.std.reshape(X.shape), shading='auto',
                       cmap='Reds', vmin=0.)
            plt.scatter(prj.x_train[:, 0], prj.x_train[:, 1], marker='x')
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
        # nextxindx = np.unravel_index(activation.argmax(), X.shape)
        # print('next position:', X[nextxindx], Y[nextxindx])

        print('adding to x_train', nextx)
        x_train = np.vstack((x_train, nextx))
        y_train = np.append(y_train, testfunc2dnrca(X, Y, density,
                                                    nextx[np.newaxis, :]))


def testfunc2dnrca(X, Y, density, xvec):
    test = (X.flatten() - xvec[0, 0])**2 + (Y.flatten() - xvec[0, 1])**2
    return density[test.argmin()]


def testfunc2d(x):
    pos2 = np.array([2., 3.])
    pos1 = np.array([1., 1.])
    return np.exp(-np.sum((x-pos1)**2., axis=1)
                  ) + 0.5*np.exp(-np.sum((x-pos2)**2., axis=1))


searchrun2dnrca()
