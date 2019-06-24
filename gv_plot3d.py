#!/usr/bin/env python
# extract the reciplocal distribution of group velocities on a reciplocal plane.
# the output format is designed to be readable by VESTA, i.e., I gave up to make the 3D plot by myself.
# Written by Kazuyoshi TATSUMI
import numpy as np
import h5py
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'


def parse_mesh(filename,  max_freq):
    f = h5py.File(filename, 'r')
    gv = f["group_velocity"]
    #gvdata = gv[:, :, 2]**2 / (gv[:, :, 0]**2 + gv[:, :, 1]**2 + gv[:, :, 2]**2)
    #gvdata = (gv[:, :, 0])**2 + (gv[:, :, 1])**2
    #gvdata = gv[:, :, 2]**2
    gvdata = abs(gv[:, :, 2])
    omega = f["frequency"][:, :]
    qp = f["qpoint"][:, :]
    mesh = f["mesh"]
    qp[qp < 0] += 1
    dim = omega.shape
    gvdata_ave = np.zeros(dim[0])
    for i in range(0, dim[0]):
        num_ele = 0
        sum_gvdata = 0
        for j in range(0, dim[1]):
            if omega[i, j] < max_freq and omega[i, j]:
                num_ele += 1
                sum_gvdata += gvdata[i, j]
        gvdata_ave[i] = sum_gvdata
    gvdata_ave[0] = 0

    return mesh, qp, gvdata_ave, omega[:, 0]


def getXYZD(qx, qy, qz, gvdata_ave, n):
    xlin = np.linspace(min(qx), max(qx), n[0])
    ylin = np.linspace(min(qy), max(qy), n[1])
    zlin = np.linspace(min(qz), max(qz), n[2])
    X, Y, Z = np.meshgrid(xlin, ylin, zlin)
    D = griddata((qx, qy, qz), gvdata_ave, (X, Y, Z), method='linear')

    D2 = np.transpose(D, (2, 1, 0))             # The output format is adjusted to the CHGCAR file.
    D3 = D2.reshape(n[0]*n[1]*n[2] / 5, 5)      # The transposing is nescessary for it. 
    return X, Y, Z, D3


def getXYZ(x, y, gvdata_sum, n):
    xlin = np.linspace(min(x), max(x), n[0]*10)
    ylin = np.linspace(min(y), max(y), n[1]*10)
    X, Y = np.meshgrid(xlin, ylin)
    Z = griddata((x, y), gvdata_sum, (X, Y), method='linear')

    return X, Y, Z


def make_diff(file1, file2, max_freq, n, rlat1, rlat2):
    mesh1, qp1, gvdata_ave1, omega1 = parse_mesh(file1, max_freq)
    mesh2, qp2, gvdata_ave2, omega2 = parse_mesh(file2, max_freq)
    X1, Y1, Z1, D1 = getXYZD(qp1[:, 0], qp1[:, 1], qp1[:, 2], gvdata_ave1, n)
    X2, Y2, Z2, D2 = getXYZD(qp2[:, 0], qp2[:, 1], qp2[:, 2], gvdata_ave2, n)
    diff = D2 * np.linalg.det(rlat2) - D1 * np.linalg.det(rlat1)
    return diff


def make_diffonbasalp(file1, file2, max_freq, n, rlat1, rlat2):
    mesh1, qp1, gvdata_ave1, omega1 = parse_mesh(file1, max_freq)
    mesh2, qp2, gvdata_ave2, omega2 = parse_mesh(file2, max_freq)
    print qp1.shape
    gvdata_ave1 = gvdata_ave1 * np.linalg.det(rlat1)
    gvdata_ave2 = gvdata_ave2 * np.linalg.det(rlat2)
    #diff = diff.reshape(mesh1)
    #diff = np.sum(diff, axis=2)
    #diff = diff.reshape(mesh1[0]*mesh1[1])
    #car = np.matmul(qp1, rlat1)
    #carx = car[:, 0].reshape(mesh1)
    #carx = carx[:, :, 0]
    #carx = carx.reshape(mesh1[0]*mesh1[1])
    #cary = car[:, 1].reshape(mesh1)
    #cary = cary[:, :, 0]
    #cary = cary.reshape(mesh1[0]*mesh1[1])
    #X, Y, Z = getXYZ(carx, cary, diff, n)
    #plt.pcolor(X,Y,Z,cmap=cm.rainbow)
    k = 0
    gvsum1 = np.zeros(((mesh1[0]+1)*(mesh1[1]+1)))
    gvsum2 = np.zeros(((mesh2[0]+1)*(mesh2[1]+1)))
    x = np.zeros(((mesh1[0]+1)*(mesh1[1]+1)))
    y = np.zeros(((mesh1[0]+1)*(mesh1[1]+1)))
    #x = np.zeros((mesh1[0]*mesh1[1]))
    #y = np.zeros((mesh1[0]*mesh1[1]))
    for i in range(0, mesh1[0]):
        for j in range(0, mesh1[1]):
            x[k] = i * 1.0/mesh1[0]
            y[k] = j * 1.0/mesh1[1]
            for l in range(0, mesh1[0]*mesh1[1]*mesh1[2]):
                if abs(x[k] - qp1[l, 0]) < 0.001 and abs(y[k] - qp1[l, 1]) < 0.001:
                    gvsum1[k] += gvdata_ave1[l]
                    gvsum2[k] += gvdata_ave2[l]
            k += 1
    for i in range(0, mesh1[0]):
        x[k] = i * 1.0/mesh1[0]
        y[k] = 1.0
        for l in range(0, mesh1[0]*mesh1[1]*mesh1[2]):
            if abs(x[k] - qp1[l, 0]) < 0.001 and abs(qp1[l, 1] - 0.000) < 0.001:
                gvsum1[k] += gvdata_ave1[l]
                gvsum2[k] += gvdata_ave2[l]
        k += 1
    for i in range(0, mesh1[1]):
        x[k] = 1.0
        y[k] = i * 1.0/mesh1[1]
        for l in range(0, mesh1[0]*mesh1[1]*mesh1[2]):
            if abs(y[k] - qp1[l, 1]) < 0.001 and abs(qp1[l, 0] - 0.000) < 0.001:
                gvsum1[k] += gvdata_ave1[l]
                gvsum2[k] += gvdata_ave2[l]
        k += 1
    x[k] = 1.0
    y[k] = 1.0
    for l in range(0, mesh1[0]*mesh1[1]*mesh1[2]):
        if abs(qp1[l, 0] - 0.000) < 0.001 and abs(qp1[l, 1] - 0.000) < 0.001:
            gvsum1[k] += gvdata_ave1[l]
            gvsum2[k] += gvdata_ave2[l]
    xy = np.zeros(((mesh1[0]+1)*(mesh1[1]+1), 2))
    xy[:, 0] = x
    xy[:, 1] = y
    ratio = gvsum2 / gvsum1
    carxy = np.matmul(xy, rlat1[0:2, 0:2])
    #X, Y, Z = getXYZ(carxy[:,0], carxy[:,1], gvsum, n)
    #X, Y, Z = getXYZ(carxy[:,0], carxy[:,1], gvsum2 - gvsum1, n)
    X, Y, Z = getXYZ(x, y, gvsum2 - gvsum1, n)
    #im=plt.pcolor(X,Y,Z,vmin=0,cmap=cm.rainbow)
    #plt.colorbar(im)
    #im = plt.contour(X, Y, Z, np.arange(0,1,0.1), colors='k')
    #im = plt.contour(X, Y, Z, np.arange(0,1,0.1))
    im = plt.contour(X, Y, Z, np.arange(30,50,1), colors='k')
    plt.clabel(im, inline=1, fontsize=8, fmt='%d')



    


def run():
    cdir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift"
    sdir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift"
    ssdir = "/home/kazu/bsi3n4_m/phonopy_doubled_334"
    max_freq = 15
    cn = np.array([50, 50, 70])
    sn = np.array([20, 20, 52])
    ssn = np.array([50, 50, 70])
    crlat = np.array([[0.128074059, 0.073943593, 0.000000000], [0.000000000,  0.147887185,  0.000000000], [0.000000000, 0.000000000, 0.1767063800]])
    srlat = np.array([[0.130517175, 0.075354126, 0.000000000], [0.000000000, 0.150708252, 0.00000000000], [0.000000000, 0.000000000, 0.3417394180]])
    ssrlat = np.array([[0.130556943, 0.075377086, 0.000000000], [0.000000000, 0.150754172, 0.00000000000], [0.000000000, 0.000000000, 0.1709573530]])
    c = cdir + "/mesh_252535.hdf5"
    s = sdir + "/mesh.hdf5"
    ss = ssdir + "/mesh_252535.hdf5"
    #For data extraction for each phase
    #qp, gvdata_ave, omega = parse_mesh(c, max_freq)
    #X, Y, Z, D3 = getXYZD(qp[:, 0], qp[:, 1], qp[:, 2], gvdata_ave, cn)
    #np.savetxt("array.txt", D3 * np.linalg.det(crlat), fmt="%17.11e")  #
    
    #Seek the difference between alpha and beta_doubled_cell
    #diff = make_diff(c, ss, max_freq, cn, crlat, ssrlat)
    make_diffonbasalp(c, ss, max_freq, cn, crlat,ssrlat)
    #np.savetxt("array.txt", diff, fmt="%17.11e")  #
    plt.show()



run()
