#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 5


def parse_band(bfile):
    f = h5py.File(bfile)
    xdata = f["distance"]
    ydata = f["frequency"]
    return(xdata, ydata)


def plotdata(x, y, i):
    ydim = y.shape
    plt.subplot(4, 2, i, aspect=0.014)
    plt.axvline(x[100], color='k')
    for j in range(0, ydim[1]):
        plt.ylim(0, 34)
        plt.xlim(min(x), max(x))
        plt.xticks([min(x), max(x)])
        plt.plot(x, y[:, j], color='k', linewidth=1.25)


def caserun(bfile, M):
    xdata, ydata = parse_band(bfile)
    x = xdata[0, :]
    y = ydata[0, :, :]
    plotdata(x, y, 1+M*2)
    x = np.concatenate((xdata[2, :], xdata[3, :]), axis=0)
    y = np.concatenate((ydata[2, :, :], ydata[3, :, :]), axis=0)
    plotdata(x, y, 2+M*2)


def run():
    plt.figure(figsize=(6, 18))
    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band.hdf5"
    caserun(cbfile, 0)
    sbfile = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/band.hdf5"
    caserun(sbfile, 1)
    gbfile = "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/band.hdf5"
    caserun(gbfile, 2)


run()
plt.show()

#plt.savefig("band-alpha-beta-gamma.eps")


