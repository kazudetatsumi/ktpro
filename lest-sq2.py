#!/usr/bin/env python

import numpy as np
from scipy import optimize
import random

np.random.seed(0)

# True parameters


def f(x, p):
    return p[0]*x + 0.4*np.sin(p[1]*x) + p[2]


p = np.array([1, 40, 2])

print('True p: ', p)

# Generate random data
xdata = np.linspace(0., 1, 120)
ydata = f(xdata, p) + np.random.normal(0., 0.2, len(xdata))

# Fits
pstart = [1, 42, 1]
errFunc = lambda p, x, y: f(x, p) - y

# Error of parameters


def errFit(hess_inv, resVariance):
    return np.sqrt(np.diag(hess_inv * resVariance))

# optimize.leastsq


fit, hess_inv, infodict, errmsg, success = optimize.leastsq(errFunc, pstart, args=(xdata, ydata), full_output=1)
dFit = errFit(hess_inv, (errFunc(fit, xdata, ydata)**2).sum()/(len(ydata)-len(pstart)))
print('leastsq:\n\tp: ', fit, '\n\tdp: ', dFit)

# optimize.curve_fit
fit, cov = optimize.curve_fit(lambda x, p0, p1, p2: f(x, [p0, p1, p2]), xdata, ydata, p0=pstart)
dFit = np.sqrt(np.diag(cov))
print('curve_fit:\n\tp: ', fit, '\n\tdp: ', dFit)

# optimize.minimize
result = optimize.minimize(lambda p, x, y: np.sum(errFunc(p, x, y)**2), pstart,  args=(xdata, ydata))
dFit = errFit(result.hess_inv,  result.fun/(len(ydata)-len(pstart)))
print('minimize:\n\tp: ', result.x, '\n\tdp: ', dFit)

# optimize.least_squares
result = optimize.least_squares(errFunc, pstart, args=(xdata, ydata))
dFit = errFit(np.linalg.inv(np.dot(result.jac.T, result.jac)), (errFunc(result.x, xdata, ydata)**2).sum()/(len(ydata)-len(pstart)))
dFit2 = errFit(np.linalg.inv(2 * np.dot(result.jac.T, result.jac)), (errFunc(result.x, xdata, ydata)**2).sum()/(len(ydata)-len(pstart)))
print('least_squares:\n\tp: ', result.x, '\n\tdp: ', dFit, '\n\tdp2: ', dFit2)
