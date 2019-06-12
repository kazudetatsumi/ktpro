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
fig, ax = plt.subplots(2, 2, figsize=(8, 14))


def parse_band(bfile):
    f = h5py.File(bfile)
    xdata = f["distance"]
    ydata = f["frequency"]
    zdata = f["eigenvector"]
    print ydata.shape
    print zdata.shape
    return(xdata, ydata, zdata)


def plotdata(x, y, z, i, j, title):
    ydim = y.shape
    ax[i, j].axvline(x[100], color='k', linewidth=0.025)
    ax[i, j].set_ylim(0, 34)
    ax[i, j].set_xlim(min(x), max(x))
    ax[i, j].set_xticks([min(x), max(x)])
    ax[i, j].set_aspect(0.009)
    ax[i, j].set_title(title)
    sc = ax[i, j].scatter(x, y, c=z, vmin=0, vmax=1, linewidth=0.01, s=1)
    #plt.colorbar(sc)
    ax[i,j].set_facecolor('k')


def caserun(bfile, M, K, title):
    xdata, ydata, zdata = parse_band(bfile)
    x = xdata[0, :]
    y = ydata[0, :, :]
    if title == "beta_Nz":
        z = (abs(zdata[0, :, 41, :]))**2 + (abs(zdata[0, :, 29, :]))**2  
    elif title == "alpha_Nz":
        z = (abs(zdata[0, :, 38, :]))**2 + (abs(zdata[0, :, 62, :]))**2 + (abs(zdata[0, :, 59, :]))**2 + (abs(zdata[0, :, 83, :]))**2 
    elif title == "alpha_Nz_near_void":
        z = (abs(zdata[0, :, 59, :]))**2 + (abs(zdata[0, :, 83, :]))**2  
    #print zdata[0, 0, 41, 0]
    #print (abs(zdata[0, 10, 10, 1]))**2    
    #print z[10,10]
    x2 = np.zeros_like(y)
    for i in range(0, x2.shape[1]):
        x2[:, i] = x
    #print x2.shape
    #print y.shape
    #print z.shape
    x3 = np.reshape(x2, (x2.shape[0]*x2.shape[1]))
    y2 = np.reshape(y, (y.shape[0]*y.shape[1]))
    z2 = np.reshape(z, (z.shape[0]*z.shape[1]))
    #print x3.shape
    #print y2.shape
    #print z2.shape

    plotdata(x3, y2, z2, M, K, title)


def run():
    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band.hdf5"
    caserun(cbfile, 0, 0, "alpha_Nz")
    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band.hdf5"
    caserun(cbfile, 0, 1, "alpha_Nz_near_void")
    sbfile = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/band.hdf5"
    caserun(sbfile, 1, 0, "beta_Nz")
    sbfile = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/band.hdf5"
    caserun(sbfile, 1, 1, "beta_Nz")
    #plt.savefig("band-alpha-beta-gamma2.eps")

run()
plt.show()
