#!/usr/bin/env python
# extract the reciplocal distribution of group velocities on a reciplocal plane.
# the output format is designed to be readable by VESTA, i.e., I gave up to make the 3D plot by myself.
# Written by Kazuyoshi TATSUMI
import numpy as np
import h5py
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import yaml

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'


def parse_mesh(filename, dn, max_freq):
    f = h5py.File(filename, 'r')
    gv = f[dn]
    gvdata = gv[:, :, 2]**2
    #gvdata = abs(gv[:, :, 2])
    #gvdata = gv[:, :, 2]**2 / (gv[:, :, 0]**2 + gv[:, :, 1]**2 + gv[:, :, 2]**2)
    #gvdata = (gv[:, :, 0])**2 + (gv[:, :, 1])**2
    omega = f["frequency"][:, :]
    mesh = f["mesh"]
    qp = f["qpoint"][:, :]
    qp[qp < 0] += 1
    dim = omega.shape
    _gvdata = np.zeros(dim)
    _gvdata[omega < max_freq] = gvdata[omega < max_freq]
    gvdata_ave = np.sum(_gvdata, axis=1)
    return(mesh, qp, gvdata_ave, omega[:, 0])


def getXYZD(qx, qy, qz, gvdata_ave, n):
    xlin = np.linspace(min(qx), max(qx), n[0]*2)
    ylin = np.linspace(min(qy), max(qy), n[1]*2)
    zlin = np.linspace(min(qz), max(qz), n[2]*2)
    X, Y, Z = np.meshgrid(xlin, ylin, zlin)
    D = griddata((qx, qy, qz), gvdata_ave, (X, Y, Z), method='linear')

    D2 = np.transpose(D, (2, 1, 0))             # The output format is adjusted to the CHGCAR file.
    D3 = D2.reshape(n[0]*n[1]*n[2] / 5, 5)      # The transposing is nescessary for it. 
    return X, Y, Z, D3


def getXYZ(x, y, gvdata_sum, n):
    xlin = np.linspace(min(x), max(x), n[0]*20)
    ylin = np.linspace(min(y), max(y), n[1]*20)
    X, Y = np.meshgrid(xlin, ylin)
    Z = griddata((x, y), gvdata_sum, (X, Y), method='linear')
    return X, Y, Z


def make_diff(file1, file2, dn, max_freq, rlat1, rlat2):
    mesh1, qp1, gvdata_ave1, omega1 = parse_mesh(file1, dn, max_freq)
    mesh2, qp2, gvdata_ave2, omega2 = parse_mesh(file2, dn, max_freq)
    X1, Y1, Z1, D1 = getXYZD(qp1[:, 0], qp1[:, 1], qp1[:, 2], gvdata_ave1, mesh1[:])
    X2, Y2, Z2, D2 = getXYZD(qp2[:, 0], qp2[:, 1], qp2[:, 2], gvdata_ave2, mesh2[:])
    diff = D2 * np.linalg.det(rlat2) - D1 * np.linalg.det(rlat1)
    return diff


def project_on_basalp(f, dn, max_freq, rlat):
    mesh, qp, gvdata_ave, omega = parse_mesh(f, dn, max_freq)
    gvdata_ave = gvdata_ave * np.linalg.det(rlat)
    k = 0
    gvsum = np.zeros(((mesh[0]+1)*(mesh[1]+1)))
    x = np.zeros(((mesh[0]+1)*(mesh[1]+1)))
    y = np.zeros(((mesh[0]+1)*(mesh[1]+1)))
    qpx = qp[:, 0]
    qpy = qp[:, 1]
    for i in range(0, mesh[0]):
        for j in range(0, mesh[1]):
            x[k] = i * 1.0/mesh[0]
            y[k] = j * 1.0/mesh[1]
            condition = abs(qpx - x[k]) + abs(qpy - y[k]) < 0.001
            _gvdata_ave = gvdata_ave[condition]
            gvsum[k] = np.sum(_gvdata_ave)
            k += 1
    for i in range(0, mesh[0]):
        x[k] = i * 1.0/mesh[0]
        y[k] = 1.0
        condition = abs(qpx - x[k]) + abs(qpy - 0.000) < 0.001
        _gvdata_ave = gvdata_ave[condition]
        gvsum[k] = np.sum(_gvdata_ave)
        k += 1
    for i in range(0, mesh[1]):
        x[k] = 1.0
        y[k] = i * 1.0/mesh[1]
        condition = abs(qpx - 0.000) + abs(qpy - y[k]) < 0.001
        _gvdata_ave = gvdata_ave[condition]
        gvsum[k] = np.sum(_gvdata_ave)
        k += 1
    x[k] = 1.0
    y[k] = 1.0
    condition = abs(qpx - 0.000) + abs(qpy - 0.000) < 0.001
    _gvdata_ave = gvdata_ave[condition]
    gvsum[k] = np.sum(_gvdata_ave)
    return x, y, mesh, gvsum


def make_diffonbasalp(file1, file2, dn, max_freq, rlat1, rlat2):
    x, y, mesh1, gvsum1 = project_on_basalp(file1, dn, max_freq, rlat1)
    x, y, mesh2, gvsum2 = project_on_basalp(file2, dn, max_freq, rlat2)
    xy = np.zeros(((mesh1[0]+1)*(mesh1[1]+1), 2))
    xy[:, 0] = x
    xy[:, 1] = y
    #ratio = gvsum2 / gvsum1
    carxy = np.matmul(xy, rlat1[0:2, 0:2])
    #X, Y, Z = getXYZ(carxy[:,0], carxy[:,1], gvsum1, mesh1[:])
    X, Y, Z = getXYZ(carxy[:,0], carxy[:,1], (gvsum2 - gvsum1) / mesh1[2], mesh1[:])
    #X, Y, Z = getXYZ(x, y, gvsum2 - gvsum1 / mesh[2], mesh1[:])
    #im=plt.pcolor(X,Y,Z,vmin=0,cmap=cm.rainbow)
    #plt.colorbar(im)
    #im = plt.contour(X, Y, Z, np.arange(0,1,0.1), colors='k')
    #im = plt.contour(X, Y, Z, np.arange(0,1,0.1))
    #im = plt.contour(X, Y, Z, np.arange(30,50,1), colors='k')
    #im = plt.contour(X, Y, Z )
    im = plt.contour(X, Y, Z, np.arange(0, 80, 2))
    plt.clabel(im, inline=1, fontsize=12, fmt='%d')


def parse_rlat(pf):
    with open(pf) as f:
        data = yaml.load(f)
        rlat = data["primitive_cell"]["reciprocal_lattice"]
    rlat = np.array(rlat)
    return rlat


def run():
    cdir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift"
    sdir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift"
    ssdir = "/home/kazu/bsi3n4_m/phonopy_doubled_334"
    ssdir = "/home/kazu/bsi3n4_m/phonopy_doubled_334_with_alphalat"
    max_freq = 15
    pc = cdir + "/phonopy.yaml"
    ps = sdir + "/phonopy.yaml"
    pss = ssdir + "/phonopy.yaml"
    crlat = parse_rlat(pc)
    srlat = parse_rlat(ps)
    ssrlat = parse_rlat(pss)
    #c = cdir + "/mesh_252535.hdf5"
    c = cdir + "/mesh.hdf5.org"
    s = sdir + "/mesh.hdf5"
    #ss = ssdir + "/mesh_252535.hdf5"
    ss = ssdir + "/mesh.hdf5.org"
    ##For data extraction for each phase
    #mesh, qp, gvdata_ave, omega = parse_mesh(c, "group_velocity", max_freq)
    #X, Y, Z, D3 = getXYZD(qp[:, 0], qp[:, 1], qp[:, 2], gvdata_ave, mesh[:])
    #np.savetxt("array.txt", D3 * np.linalg.det(srlat), fmt="%17.11e")  #
    
    ##Seek the difference between alpha and beta_doubled_cell
    #diff = make_diff(c, ss, "group_velocity", max_freq, crlat, ssrlat)
    make_diffonbasalp(c, ss, "group_velocity", max_freq, crlat, ssrlat)
    #np.savetxt("array.txt", diff, fmt="%17.11e")  #
    plt.show()



run()
