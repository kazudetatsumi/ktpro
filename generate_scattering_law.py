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
plt.rcParams['font.family'] = 'Times New Roman'
cmtomeV = 0.123984  # [meV cm]
# twomoverhbarsq_windsor=1/(0.81785*100) # [meV-1 angs-2]
E0 = 2419.66*cmtomeV  # [meV]
sites = ["3f"]
# sites = ["12n","6m","4h","3f","12o"]
ntheta = 135-2
nthetaperfile = 21
nthetainlastfile = ntheta % nthetaperfile
nfile = ntheta // nthetaperfile + 1
head = "./"
csvname = "out_direct_inc_293K.csv"
mass = 1.674927471*10**(-27)  # [kg]
hbar = 1.054571817*10**(-34)  # [J s]
twomoverhbarsq = 2*mass/((hbar**2) * (6.2415093418967*(10**21)) * (10**20))
minq = 1  # [angs-1]
maxq = 8  # [angs=1]
plt.rcParams['font.family'] = 'Times New Roman'
fig = plt.figure(figsize=(4, 5))


def initialize(data):
    datasize = data.shape
    nhw = datasize[0]
    hw = data[:, 0]*cmtomeV  # [meV]
    alldata = np.zeros((len(sites), nhw*ntheta, 3))  # 1st axis -> site, 2nd -> samplepoint(q, hw), 3rd ->  0:hw, 1:q, 2:intensity
    alldata[:, :, 0] = np.tile(hw, (len(sites), ntheta))
    return alldata


def get_data(sites, alldata):
    sidx = 0
    for site in sites:
        for fidx in range(1, nfile+1):
            infile = head + site + "/" + csvname+"."+str(fidx)
            data = np.genfromtxt(infile, dtype=float, delimiter=",")
            nhw = data.shape[0]
            ntidx = data.shape[1]
            for tidx in range(1, ntidx):
                theta = tidx + (fidx - 1)*nthetaperfile + 2
                sin = math.sin(math.radians(theta))
                cos = math.cos(math.radians(theta))
                q = (twomoverhbarsq * (2 * E0 - data[:, 0]*cmtomeV -
                     2 * cos * ((E0 * (E0 - data[:, 0]*cmtomeV))**0.5)))**0.5
                intensity = data[:, tidx]
                ii = (fidx-1)*nthetaperfile*nhw + (tidx - 1)*nhw
                ff = ii+nhw
                alldata[sidx, ii:ff, 1] = q
                alldata[sidx, ii:ff, 2] = intensity*sin
        sidx += 1
    return(alldata)


def save_hdf(data, outfile):
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('alldata', data=data)


def read_hdf(outfile):
    f = h5py.File(outfile)
    return f["alldata"]


def getXYZ(x, y, z, xb, yb):
    ngrid = 400
    xlin = np.linspace(xb[0], xb[1], ngrid)
    ylin = np.linspace(yb[0], yb[1], ngrid)
    X, Y = np.meshgrid(xlin, ylin)
    Z = griddata((x, y), z, (X, Y), method='linear', fill_value=-1)
    return X, Y, Z


def plot_sub(X, Y, Z, axindx, vmax):
    ax = fig.add_subplot(2, 1, axindx)
    cmap1 = copy.copy(plt.cm.jet)
    cmap1.set_under('w')
    ax.pcolor(X, Y, Z, cmap=cmap1, vmin=0, vmax=vmax)
    ax.axis('tight')
    if axindx == 2:
        ax.set_xlabel('hw (meV)')
        ax.set_yticks([1, 8])
    ax.set_ylabel('q (Angs-1)')


def plot_scatteringlaw(data):
    X, Y, Z = getXYZ(data[:, 0], data[:, 1], data[:, 2], [min(data[:, 0]),
                     max(data[:, 0])], [min(data[:, 1]), max(data[:, 1])])
    plot_sub(X, Y, Z, 1, max(data[:, 2])/8.0)
    X, Y, Z = getXYZ(data[:, 0], data[:, 1], data[:, 2], [min(data[:, 0]),
                     max(data[:, 0])], [minq, maxq])
    plot_sub(X, Y, Z, 2, max(data[:, 2])/8.0)


def run():
    infile = head+"3f"+"/"+csvname+".1"
    outfile = "alldata.hdf5"
    #datatmp = np.genfromtxt(infile, dtype=float, delimiter=",")
    #alldata = initialize(datatmp)
    #alldata = get_data(sites, alldata)
    #save_hdf(alldata, outfile)
    alldata = read_hdf(outfile)

    plot_scatteringlaw(alldata[0, :, :])
    # plt.show()
    plt.savefig("scattering_law_3f_model_lda.png")


run()
