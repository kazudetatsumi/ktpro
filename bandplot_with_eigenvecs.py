#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import yaml

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.major.size'] = 5
fig, ax = plt.subplots(2, 4, figsize=(16, 14))


def parse_band(bfile):
    f = h5py.File(bfile)
    xdata = f["distance"]
    ydata = f["frequency"]
    zdata = f["eigenvector"]
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
    ax[i, j].set_facecolor('k')


def get_normalvec(pfile, ic, bondlenlim):
    with open(pfile) as f:
        data = yaml.load(f)
        celldata = data["primitive_cell"]

    latvec = np.array(celldata["lattice"])
    numa = len(celldata["points"])
    frac = np.zeros((numa, 3))
    for i in range(0, numa):
        frac[i, :] = np.array(celldata["points"][i]["coordinates"])

    fracs = np.zeros((numa*27, 3))
    l = 0
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                for m in range(0, numa):
                    fracs[l, :] = np.array([i, j, k]) + frac[m, :]
                    l += 1
    cars = np.matmul(fracs, latvec)

    car = np.matmul(frac, latvec)
    ca_cars = np.zeros_like(cars)
    for i in range(0, numa*27):
        ca_cars[i, :] = car[ic, :]

    diff_cars = cars - ca_cars
    dist = (diff_cars[:, 0]**2 + diff_cars[:, 1]**2 + diff_cars[:, 2]**2)**0.5
    nn_ids = []
    for i in range(0, numa*27):
        if dist[i] <= bondlenlim and dist[i] > 0.00001:
            nn_ids += [i]

    diff1 = cars[nn_ids[0], :] - cars[nn_ids[1], :]
    diff2 = cars[nn_ids[0], :] - cars[nn_ids[2], :]
    normalvec = np.cross(diff1, diff2)
    normalvec /= np.linalg.norm(normalvec, 2)

    return normalvec


def get_z(pfile, atomlist, direction, zdata, bl):
    z = np.zeros_like(zdata[0, :, 0, :])
    if direction == "perp":
        for an in atomlist:
            nv = get_normalvec(pfile, an - 1, bl)
            for i in range(0, 3):
                z += abs(zdata[0, :, an * 3 - 3 + i, :] * nv[i])**2
    elif direction == "z":
        for an in atomlist:
            z += (abs(zdata[0, :, an * 3 - 1, :]))**2
    elif direction == "xy":
        for an in atomlist:
            z += (abs(zdata[0, :, an * 3 - 2, :]))**2
            z += (abs(zdata[0, :, an * 3 - 3, :]))**2
    elif direction == "para":
        for an in atomlist:
            nv = get_normalvec(pfile, an - 1, bl)
            for i in range(0, 3):
                z -= abs(zdata[0, :, an * 3 - 3 + i, :] * nv[i])**2
                z += (abs(zdata[0, :, an * 3 - 3 + i, :]))**2
            #z += (abs(zdata[0, :, an * 3 - 2, :]))**2
            #z += (abs(zdata[0, :, an * 3 - 3, :]))**2
    return z


def caserun(bfile, pfile, M, K, title, bondlenlim):
    xdata, ydata, zdata = parse_band(bfile)
    x = xdata[0, :]
    y = ydata[0, :, :]
    if title == "alpha_Nperp":
        atomlist = [13, 20, 21, 28]
        z = get_z(pfile, atomlist, 'perp', zdata, bondlenlim)
    elif title == "alpha_Npara":
        atomlist = [13, 20, 21, 28]
        z = get_z(pfile, atomlist, 'para', zdata, bondlenlim)
    elif title == "alpha_Nperp_II":
        atomlist = [14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27]
        z = get_z(pfile, atomlist, 'perp', zdata, bondlenlim)
    elif title == "alpha_Npara_II":
        atomlist = [14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27]
        z = get_z(pfile, atomlist, 'para', zdata, bondlenlim)
    elif title == "alpha_Nz":
        atomlist = [13, 20, 21, 28]
        z = get_z(pfile, atomlist, 'z', zdata, bondlenlim)
    elif title == "alpha_Nxy":
        atomlist = [13, 20, 21, 28]
        z = get_z(pfile, atomlist, 'xy', zdata, bondlenlim)
    elif title == "beta_Nperp":
        atomlist = [10, 14]
        z = get_z(pfile, atomlist, 'perp', zdata, bondlenlim)
    elif title == "beta_Npara":
        atomlist = [10, 14]
        z = get_z(pfile, atomlist, 'para', zdata, bondlenlim)
    elif title == "beta_Nperp_II":
        atomlist = [7, 8, 9, 11, 12, 13]
        z = get_z(pfile, atomlist, 'perp', zdata, bondlenlim)
    elif title == "beta_Npara_II":
        atomlist = [7, 8, 9, 11, 12, 13]
        z = get_z(pfile, atomlist, 'para', zdata, bondlenlim)
    elif title == "beta_Nz":
        atomlist = [10, 14]
        z = get_z(pfile, atomlist, 'z', zdata, bondlenlim)
    elif title == "beta_Nxy":
        atomlist = [10, 14]
        z = get_z(pfile, atomlist, 'xy', zdata, bondlenlim)
    x2 = np.zeros_like(y)
    for i in range(0, x2.shape[1]):
        x2[:, i] = x
    x3 = np.reshape(x2, (x2.shape[0]*x2.shape[1]))
    y2 = np.reshape(y, (y.shape[0]*y.shape[1]))
    z2 = np.reshape(z, (z.shape[0]*z.shape[1]))

    plotdata(x3, y2, z2, M, K, title)




def run():
    bondlenlim = 2.0
    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band.hdf5"
    cpfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/phonopy.yaml"
    caserun(cbfile, cpfile, 0, 0, "alpha_Nperp", bondlenlim)
    caserun(cbfile, cpfile,  0, 1, "alpha_Npara", bondlenlim)
    caserun(cbfile, cpfile, 0, 2, "alpha_Nz", bondlenlim)
    caserun(cbfile, cpfile, 0, 3, "alpha_Nxy", bondlenlim)
    #caserun(cbfile, cpfile, 0, 2, "alpha_Nperp_II", bondlenlim)
    #caserun(cbfile, cpfile, 0, 3, "alpha_Npara_II", bondlenlim)
    sbfile = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/band.hdf5"
    spfile = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/phonopy.yaml"
    #normalvec = get_normalvec(spfile, 10, bondlenlim)
    #print normalvec
    caserun(sbfile, spfile, 1, 0, "beta_Nperp", bondlenlim)
    caserun(sbfile, spfile, 1, 1, "beta_Npara", bondlenlim)
    caserun(sbfile, spfile, 1, 2, "beta_Nz", bondlenlim)
    caserun(sbfile, spfile, 1, 3, "beta_Nxy", bondlenlim)
    #caserun(sbfile, spfile, 1, 2, "beta_Nperp_II", bondlenlim)
    #caserun(sbfile, spfile, 1, 3, "beta_Npara_II", bondlenlim)
    #plt.savefig("band-alpha-beta-gamma2.eps")

run()
plt.show()
