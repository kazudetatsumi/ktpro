#!/usr/bin/env python
# This script plot 2D map of S(q, w) from the OCLIMAX output csv files each of which has a different scattering angle.
# The paths of the csv files are hard coded and u should alter it for proper operation.
# 2020/03/10 Kazuyoshi TATSUMI
import matplotlib.pyplot as plt
import numpy as np
import math
import copy
import h5py
from scipy.interpolate import griddata
fig = plt.figure(figsize=(10, 8))
cmtomeV = 0.123984  # [meV cm]


def get_data(infile):
    f = open(infile)
    q = []
    omega = []
    intensity = []
    for i, line in enumerate(f):
        if i >= 5:
            #intensity.append(float(line.split()[2]))
            #omega.append(float(line.split()[1]))
            #q.append(float(line.split()[0]))
            if len(line.split(",")) >= 3:
               intensity.append(float(line.split(",")[2]))
               omega.append(float(line.split(",")[1]))
               q.append(float(line.split(",")[0]))
    intensity = np.array(intensity)
    omega = np.array(omega)
    q = np.array(q)
    data = np.zeros((q.shape[0],3))
    #data[:, 0] = omega*cmtomeV
    data[:, 0] = omega
    data[:, 1] = q
    data[:, 2] = intensity
    return(data)


def get_data_csv(infile):
    data = np.genfromtxt(infile, dtype=float, delimiter=",")
    return(data)


def save_hdf(data, outfile):
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('alldata', data=data)


def read_hdf(outfile):
    f = h5py.File(outfile)
    return f["alldata"]


def getXYZ(x, y, z, xb, yb):
    ngrid = 1000
    xlin = np.linspace(xb[0], xb[1], ngrid)
    ylin = np.linspace(yb[0], yb[1], ngrid)
    X, Y = np.meshgrid(xlin, ylin)
    Z = griddata((x, y), z, (X, Y), method='linear', fill_value=-1)
    return X, Y, Z


def plot_sub(X, Y, Z, axindx, vmax):
    ax = fig.add_subplot(1, 1, axindx)
    cmap1 = copy.copy(plt.cm.jet)
    cmap1.set_under('w')
    ax.pcolor(X, Y, Z, cmap=cmap1, vmin=0, vmax=vmax)
    ax.axis('tight')
    if axindx == 2:
       ax.set_ylabel('hw (meV)')
       #ax.set_yticks([1, 8])
       ax.set_xlabel('q (Angs-1)')


def plot_scatteringlaw(data):
    X, Y, Z = getXYZ(data[:, 0], data[:, 1], data[:, 2], [min(data[:, 0]),
                     max(data[:, 0])], [min(data[:, 1]), max(data[:, 1])])
    plot_sub(Y, X, Z, 1, max(data[:, 2])/10.0)
    #plot_sub(Y, X, Z, 1, 0.5376)
    print(max(data[:, 2]))
    #X, Y, Z = getXYZ(data[:, 0], data[:, 1], data[:, 2], [min(data[:, 0]),
    #                 max(data[:, 0])], [min(data[:, 1]), max(data[:, 1])])
    #                 #max(data[:, 0])], [minq, maxq])
    #plot_sub(X, Y, Z, 2, max(data[:, 2])/1.0)


def run():
    #infile = "al_file.sqw"
    infile = "out_2Dmesh_coh_50K_sigma2.2meV.csv"
    infile = "Cu_disp_H11_2Dmesh_scqw_300K.csv"
    infile = "Cu_disp_tst_2Dmesh_scqw_300K.csv"
    data = get_data(infile)
    #data = get_data_csv(infile)

    plot_scatteringlaw(data)
    plt.show()
    #plt.savefig("out_2Dmesh_coh_298K.png")
    #plt.savefig("Cu_disp_H11_2Dmesh_scqw_300K.png")
    #plt.savefig("Cu_disp_tst_2Dmesh_scqw_300K.png")


run()
