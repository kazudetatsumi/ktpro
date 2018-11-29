#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.major.size'] = 5
fig, ax = plt.subplots(3, 2, figsize=(6, 10))


def parse_band(bfile):
    f = h5py.File(bfile)
    xdata = f["distance"]
    ydata = f["frequency"]
    return(xdata, ydata)


def plotdata(x, y, i, j):
    ydim = y.shape
    for k in range(0, ydim[1]):
        if j == 1:
            ax[i, j].axvline(x[100], color='k', linewidth=0.025)
        ax[i, j].set_ylim(0, 34)
        ax[i, j].set_xlim(min(x), max(x))
        ax[i, j].set_xticks([min(x), max(x)])
        ax[i, j].set_aspect(0.009)
        ax[i, j].plot(x, y[:, k], color='k', linewidth=1)


def caserun(bfile, M):
    xdata, ydata = parse_band(bfile)
    x = xdata[0, :]
    y = ydata[0, :, :]
    plotdata(x, y, M, 0)
    x = np.concatenate((xdata[2, :], xdata[3, :]), axis=0)
    y = np.concatenate((ydata[2, :, :], ydata[3, :, :]), axis=0)
    plotdata(x, y, M, 1)


def run():
    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band.hdf5"
    caserun(cbfile, 0)
    sbfile = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/band.hdf5"
    caserun(sbfile, 1)
    gbfile = "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/band.hdf5"
    caserun(gbfile, 2)
    plt.savefig("band-alpha-beta-gamma2.eps")

run()
plt.show()
