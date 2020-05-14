#!/usr/bin/env python

import h5py
import numpy as np
import matplotlib.pyplot as plt
import yaml
import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy import stats

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.major.size'] = 5
fig, ax = plt.subplots(2, 2, figsize=(14, 14))


def parse_band(bfile):
    f = h5py.File(bfile)
    xdata = f["distance"]
    ydata = f["frequency"]
    zdata = f["eigenvector"]
    return(xdata, ydata, zdata)


def scatterplotdata(x, y, z, i, j, title):
    ydim = y.shape
    ax[i, j].axvline(x[100], color='k', linewidth=0.025)
    #ax[i, j].set_ylim(0, 34)
    #ax[i, j].set_ylim(200, 650)
    ax[i, j].set_xlim(min(x), max(x))
    ax[i, j].set_xticks([min(x), max(x)])
    #ax[i, j].set_aspect(0.009)
    ax[i, j].set_aspect(0.009/THztoKayser)
    ax[i, j].set_title(title)
    x2 = []
    y2 = []
    z2 = []
    for xx, yy, zz in zip(x, y, z):
        condition = zz > 0.25
        _x = np.extract(condition, xx)
        _y = np.extract(condition, yy)
        _z = np.extract(condition, zz)
        x2 += list(_x)
        y2 += list(_y)
        z2 += list(_z)
    z2 =np.array(z2)
    x2 =np.array(x2)
    y2 =np.array(y2)
    print max(z2)
    #sc = ax[i, j].scatter(x2, y2, c=z2, vmin=0, vmax=1,  s=z2*100, cmap='jet')
    #sc = ax[i, j].scatter(x2, y2, vmin=0, vmax=1,  s=z2*100, cmap='jet')
    print type(z2)
    print type(x2)
    sc = ax[i, j].scatter(x2, y2, c='k', linewidth=1, marker=".", s=(3*z2)**3)
    #plt.colorbar(sc)
    #ax[i, j].set_facecolor('k')
    print "sum_of_z_values:", np.sum(z)


def lineplotdata(x, y, z, i, j, title):
    ydim = y.shape
    #ax[i, j].axvline(x[100], color='k', linewidth=0.025)
    ax[i, j].set_ylim(0, 34*THztoKayser)
    #ax[i, j].set_ylim(150, 600)
    #ax[i, j].set_ylim(100, 700)
    #ax[i, j].set_ylim(200, 650)
    #ax[i, j].set_ylim(350, 550)
    ax[i, j].set_ylim(0, 600)
    ax[i, j].set_xlim(min(x), max(x))
    ax[i, j].set_xticks([min(x), max(x)])
    ax[i, j].set_aspect(0.009/THztoKayser)
    ax[i, j].set_title(title)
    #sc = ax[i, j].scatter(x, y, c=z, vmin=0, vmax=1, linewidth=0.01, s=1, cmap='binary')
    #sc = ax[i, j].scatter(x, y, c=z, vmin=0, vmax=1, linewidth=0.01, s=0.2)
    print "min",np.min(z)
    print "max",np.max(z)
    #z2 = cm.jet((z-np.min(z))/(np.max(z)-np.min(z)))
    z2 = cm.jet(z)
    zsum = np.sum(z, axis=0)
    print np.max(zsum)
    for l in range(0, y.shape[1]):
        #if zsum[l] >= np.max(zsum)*0.30:
        #if zsum[l] >= 30:
            y2 = y[:, l]
            z3 = z2[:, l, :]
            for k in np.arange(x.shape[0]-1):
                if z[k,l] > 0.20:
                    im = ax[i, j].plot([x[k],x[k+1]], [y2[k],y2[k+1]], c=z3[k], linewidth=0.5 )

    ax[i, j].set_facecolor('k')
    #from matplotlib.colors import Normalize
    #norm = Normalize(vmin=0, vmax=1)
    #from matplotlib.cm import ScalarMappable, get_cmap
    #cmap = get_cmap("jet")
    #mappable=ScalarMappable(norm=norm,cmap=cmap)
    #mappable._A = []
    #plt.colorbar(mappable)
    print "sum_of_z_values:",np.sum(z)


def get_normalvec(celldata, ic, bondlenlim):
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
    ca_cars = np.matmul(np.ones((numa*27, 1)), np.reshape(car[ic, :], (1, 3)))

    diff_cars = cars - ca_cars
    dist = (diff_cars[:, 0]**2 + diff_cars[:, 1]**2 + diff_cars[:, 2]**2)**0.5
    t = np.where((dist <= bondlenlim) & (dist > 0.1))[0][:]

    diff1 = cars[t[0], :] - cars[t[1], :]
    diff2 = cars[t[0], :] - cars[t[2], :]
    normalvec = np.cross(diff1, diff2)
    normalvec /= np.linalg.norm(normalvec, 2)

    return normalvec


def get_z(celldata, atomlist, direction, zdata, bl):
    z = np.zeros(zdata[0, :, 0, :].shape)
    if direction == "perp":
        for an in atomlist:
            nv = get_normalvec(celldata, an - 1, bl)
            for i in range(0, 3):
                z += (np.abs(zdata[0, :, an * 3 - 3 + i, :] * nv[i]))**2
    elif direction == "z":
        for an in atomlist:
            z += (np.abs(zdata[0, :, an * 3 - 1, :]))**2
    elif direction == "xy":
        for an in atomlist:
            z += (np.abs(zdata[0, :, an * 3 - 2, :]))**2
            z += (np.abs(zdata[0, :, an * 3 - 3, :]))**2
    elif direction == "para":
        for an in atomlist:
            nv = get_normalvec(celldata, an - 1, bl)
            for i in range(0, 3):
                z -= np.real_if_close((np.abs(zdata[0, :, an * 3 - 3 + i, :] * nv[i]))**2, 0.1)
                z += np.real_if_close(((np.abs(zdata[0, :, an * 3 - 3 + i, :])))**2, 0.1)
    return z


def caserun(xdata, ydata, zdata, celldata, M, K, title, bondlenlim):
    x = xdata[0, :]
    y = ydata[0, :, :]*THztoKayser
    if title == "alpha_Nperp":
        atomlist = [13, 20, 21, 28]
        z = get_z(celldata, atomlist, 'perp', zdata, bondlenlim)
    elif title == "alpha_Npara":
        atomlist = [13, 20, 21, 28]
        z = get_z(celldata, atomlist, 'para', zdata, bondlenlim)
    elif title == "alpha_Nperp_II":
        atomlist = [14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27]
        z = get_z(celldata, atomlist, 'perp', zdata, bondlenlim)
    elif title == "alpha_Npara_II":
        atomlist = [14, 15, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27]
        z = get_z(celldata, atomlist, 'para', zdata, bondlenlim)
    elif title == "alpha_Nz":
        atomlist = [13, 20, 21, 28]
        z = get_z(celldata, atomlist, 'z', zdata, bondlenlim)
    elif title == "alpha_Nxy":
        atomlist = [13, 20, 21, 28]
        z = get_z(celldata, atomlist, 'xy', zdata, bondlenlim)
    elif title == "beta_Nperp":
        #atomlist = [10, 14]
        atomlist = [19, 20, 27, 28]
        z = get_z(celldata, atomlist, 'perp', zdata, bondlenlim)
    elif title == "beta_Npara":
        #atomlist = [10, 14]
        atomlist = [19, 20, 27, 28]
        z = get_z(celldata, atomlist, 'para', zdata, bondlenlim)
    elif title == "beta_Nperp_II":
        #atomlist = [7, 8, 9, 11, 12, 13]
        atomlist = [13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26]
        z = get_z(celldata, atomlist, 'perp', zdata, bondlenlim)
    elif title == "beta_Npara_II":
        #atomlist = [7, 8, 9, 11, 12, 13]
        atomlist = [13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26]
        z = get_z(celldata, atomlist, 'para', zdata, bondlenlim)
    elif title == "beta_Nz":
        atomlist = [10, 14]
        z = get_z(celldata, atomlist, 'z', zdata, bondlenlim)
    elif title == "beta_Nxy":
        atomlist = [10, 14]
        z = get_z(celldata, atomlist, 'xy', zdata, bondlenlim)
    x2 = np.zeros_like(y)
    for i in range(0, x2.shape[1]):
        x2[:, i] = x
    x3 = np.reshape(x2, (x2.shape[0]*x2.shape[1]))
    y2 = np.reshape(y, (y.shape[0]*y.shape[1]))
    z2 = np.reshape(z, (z.shape[0]*z.shape[1]))
    #print x.shape
    print y.shape
    #print x2.shape

    #scatterplotdata(x3, y2, z2, M, K, title) # for scatter plot
    lineplotdata(x, y, z, M, K, title)     # for line plot


def get_files(bfile, pfile):
    xdata, ydata, zdata = parse_band(bfile)
    with open(pfile) as f:
        data = yaml.load(f)
        celldata = data["primitive_cell"]
    return xdata, ydata, zdata, celldata


def getXYZ(x, y, z, n):
    xlin = np.linspace(min(x), max(x), n[0])
    ylin = np.linspace(min(y), max(y), n[1])
    X, Y = np.meshgrid(xlin, ylin)
    Z = griddata((x, y), z, (X, Y), method='linear')
    return X, Y, Z


def binimage(x, y, nb):
    x = np.squeeze(x)
    y = np.squeeze(y)
    ylin = np.linspace(-1.0, 40.0, nb+1)
    nsize = y.shape
    nc = np.zeros((nsize[0], nb))
    values = np.ones(nsize[1])
    for i in range(0, nsize[0]):
        nc[i, :], res2, res3 = stats.binned_statistic(y[i, :], values, "sum", ylin)

    X, Y = np.meshgrid(x, ylin)
    nc = np.transpose(nc)
    for i in range(0, nsize[0]):
        print "CHK", np.sum(nc[:, i])
    plt.figure()
    plt.pcolor(X, Y, nc, cmap=cm.gray)





def run():
    bondlenlim = 2.0
    nybin = 40

    #cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band.hdf5"
    #cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_00501.hdf5"
    #cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_055045.hdf5"
    #cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_0505.hdf5"
    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_4-14_1-14.hdf5"
    #cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_031006.hdf5"
    #cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_032032.hdf5"
    #cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_K-H.hdf5"
    #cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_043015.hdf5"
    cpfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/primitive.yaml"
    cxdata, cydata, czdata, ccelldata = get_files(cbfile, cpfile)
    binimage(cxdata, cydata, nybin)
    caserun(cxdata, cydata, czdata, ccelldata, 0, 0, "alpha_Nperp", bondlenlim)
    #caserun(cxdata, cydata, czdata, ccelldata,  0, 1, "alpha_Npara", bondlenlim)
    #caserun(cxdata, cydata, czdata, ccelldata, 0, 2, "alpha_Nz", bondlenlim)
    #caserun(cxdata, cydata, czdata, ccelldata, 0, 3, "alpha_Nxy", bondlenlim)
    ##caserun(cxdata, cydata, czdata, ccelldata, 0, 1, "alpha_Nperp_II", bondlenlim)
    #caserun(cxdata, cydata, czdata, ccelldata, 0, 3, "alpha_Npara_II", bondlenlim)

    #sbfile = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/band.hdf5"
    #spfile = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/primitive.yaml"
    #sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band.hdf5"
    #sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_00501.hdf5"
    #sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_055045.hdf5"
    #sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_0505.hdf5"
    sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_4-14_1-14.hdf5"
    #sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_031006.hdf5"
    #sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_032032.hdf5"
    #sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_K-H.hdf5"
    #sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_043015.hdf5"
    spfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/primitive.yaml"
    sxdata, sydata, szdata, scelldata = get_files(sbfile, spfile)
    caserun(sxdata, sydata, szdata, scelldata, 1, 0, "beta_Nperp", bondlenlim)
    binimage(sxdata, sydata, nybin)
    #caserun(sxdata, sydata, szdata, scelldata, 1, 1, "beta_Npara", bondlenlim)
    #caserun(sxdata, sydata, szdata, scelldata, 1, 2, "beta_Nz", bondlenlim)
    #caserun(sxdata, sydata, szdata, scelldata, 1, 3, "beta_Nxy", bondlenlim)
    ##caserun(sxdata, sydata, szdata, scelldata, 1, 1, "beta_Nperp_II", bondlenlim)
    #caserun(sxdata, sydata, szdata, scelldata, 1, 3, "beta_Npara_II", bondlenlim)
    #plt.savefig("band-alpha-beta-gamma2.eps")
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=1)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("jet")
    mappable=ScalarMappable(norm=norm,cmap=cmap)
    mappable._A = []
    plt.colorbar(mappable)

THztoKayser = 33.35641
run()
plt.show()
