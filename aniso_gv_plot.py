#!/usr/bin/env python
# plot the anisotropic information of group velocities on a reciplocal plane.
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
    aniso = gv[:, :, 2]**2 / (gv[:, :, 0]**2 + gv[:, :, 1]**2 + gv[:, :, 2]**2)
    omega = f["frequency"][:, :]
    qp = f["qpoint"][:, :]
    qp[qp < 0] += 1
    dim = omega.shape
    #for i in range(0, dim[0]):
    #    print qp[i, :]
    aniso_ave = np.zeros(dim[0])
    for i in range(0, dim[0]):
        num_ele = 0
        sum_ani = 0
        for j in range(0, dim[1]):
            if omega[i, j] < max_freq:
                num_ele += 1
                sum_ani += aniso[i, j]
        aniso_ave[i] = sum_ani/num_ele
    aniso_ave[0] = 0

    return qp, aniso_ave, omega[:, 0]


def select_slice(qp, aniso_ave, y):
    qxs = []
    qys = []
    qzs = []
    anis = []
    for qx, qy, qz, ani in zip(qp[:, 0], qp[:, 1], qp[:, 2], aniso_ave):
        condition = abs(qy - y) < 0.01
        _qx = np.extract(condition, qx)
        _qz = np.extract(condition, qz)
        _ani = np.extract(condition, ani)
        _qy = np.extract(condition, qy)
        qxs += list(_qx)
        qys += list(_qy)
        qzs += list(_qz)
        anis += list(_ani)
    qx1 = np.array(qxs)
    qy1 = np.array(qys)
    qz1 = np.array(qzs)
    aniso_ave1 = np.array(anis)
    #print qx1[:]
    #print qz1[:]

    return qx1, qz1, aniso_ave1


def getXYZ(qx, qz, ani):
    #x = qx*0.1306
    #y = qz*0.342
    x = qx
    y = qz
    xlin = np.linspace(min(x), max(x), 200)
    ylin = np.linspace(min(y), max(y), 400)
    X, Y = np.meshgrid(xlin, ylin)
    Z = griddata((x, y), ani, (X, Y), method='linear')

    return X, Y, Z


def run():
    cdir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift"
    sdir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift"
    max_freq = 15
    y_s = 0.4
    c = cdir + "/mesh.hdf5"
    s = sdir + "/mesh.hdf5"

    qp, aniso_ave, omega = parse_mesh(s, max_freq)
    qx1, qz1, aniso_ave1 = select_slice(qp, aniso_ave, y_s)
    #qx1, qz1, aniso_ave1 = select_slice(qp, omega, y_s)

    X, Y, Z = getXYZ(qx1, qz1, aniso_ave1)
    fig = plt.figure(figsize=(12, 12))
    #cf = plt.pcolor(X, Y, Z, vmin=0, vmax=10, cmap=cm.rainbow)
    cf = plt.pcolor(X, Y, Z, vmin=0, vmax=1, cmap=cm.rainbow)
    plt.colorbar(cf)
    #im = plt.contour(X, Y, Z)
    plt.show()
    #print anisos1.shape
    #print np.average(anisos1,axis=0)






run()
