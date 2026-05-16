#!/usr/bin/env python
import numpy as np
from matplotlib import pyplot as plt


def fun_gaus(x, sigma, mean):
    return np.exp(-np.power((x - mean)/sigma, 2.)/2.)


def fun_lore(x, gamma, mean):
    dx = x[1] - x[0]
    l = (1 / 3.14) * 0.5 * gamma / ((x - mean)**2 + (0.5 * gamma)**2) 
    l = l/(np.sum(l)*dx)
    return l


def run():
    n = 100
    gamma = 0.2
    omegas = np.array([0.1, 0.5, 0.7])

    x = np.arange(-100, 200)*0.01
    dx = x[1] - x[0]

    ll = np.zeros_like(x)

    for o in omegas:
        ll += fun_lore(x, gamma, o)

    ll = ll/(np.sum(ll)*dx)

    k = np.random.poisson(n*ll)

    plt.figure(figsize=(16, 8))
    plt.plot(x, ll*n)
    plt.scatter(x, k)


run()
plt.show()
