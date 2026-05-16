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
    #aniso = gv[:, :, 2]**2 / (gv[:, :, 0]**2 + gv[:, :, 1]**2 + gv[:, :, 2]**2)
    #aniso = ((gv[:, :, 0])**2 + (gv[:, :, 1])**2)**0.5
    aniso = abs(gv[:, :, 2])
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
        #aniso_ave[i] = sum_ani/num_ele
        aniso_ave[i] = sum_ani
    aniso_ave[0] = 0

    return qp, aniso_ave, omega[:, 0]


def select_slice(qp, aniso_ave, y):
    qxs = []
    qys = []
    qzs = []
    anis = []
    for qx, qy, qz, ani in zip(qp[:, 0], qp[:, 1], qp[:, 2], aniso_ave):
        condition = abs(qy - y) < 0.01
        #condition2 = qx <= 0.5
        #condition3 = qz <= 0.5
        #_qx = np.extract(condition & condition2 & condition3, qx)
        #_qz = np.extract(condition & condition2 & condition3, qz)
        #_ani = np.extract(condition & condition2 & condition3, ani)
        #_qy = np.extract(condition & condition2 & condition3, qy)
        _qx = np.extract(condition, qx)
        _qz = np.extract(condition, qz)
        _ani = np.extract(condition, ani)
        _qy = np.extract(condition, qy)
        #_qz = np.extract(condition, qz)
        #_ani = np.extract(condition, ani)
        #_qy = np.extract(condition, qy)
        qxs += list(_qx)
        qys += list(_qy)
        qzs += list(_qz)
        anis += list(_ani)
        if (abs(_qx) < 0.005):
            qxs += list(_qx+1)
            qys += list(_qy)
            qzs += list(_qz)
            anis += list(_ani)
        if (abs(_qz) <  0.005):
            qxs += list(_qx)
            qys += list(_qy)
            qzs += list(_qz+1.0)
            anis += list(_ani)
        if (abs(_qx) <  0.005) & (abs(_qz) <  0.005):
            qxs += list(_qx+1.0)
            qys += list(_qy)
            qzs += list(_qz+1.0)
            anis += list(_ani)
        #elif (abs(_qz) == 0) & (abs(_qx) == 0):
        #    qxs += list(_qx+1.0)
        #    qys += list(_qy)
        #    qzs += list(_qz+1.0)
        #    anis += list(_ani)
    qx1 = np.array(qxs)
    qy1 = np.array(qys)
    qz1 = np.array(qzs)
    aniso_ave1 = np.array(anis)
    print qx1[:]
    print qz1[:]

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
    y_s = 0.0
    c = cdir + "/mesh.hdf5"
    s = sdir + "/mesh.hdf5"

    qp, aniso_ave, omega = parse_mesh(c, max_freq)
    qx1, qz1, aniso_ave1 = select_slice(qp, aniso_ave, y_s)
    #qx1, qz1, aniso_ave1 = select_slice(qp, omega, y_s)

    X, Y, Z = getXYZ(qx1, qz1, aniso_ave1)
    fig = plt.figure(figsize=(12, 12))
    #cf = plt.pcolor(X, Y, Z, vmin=0, vmax=10, cmap=cm.rainbow)
    #cf = plt.pcolor(X, Y, Z, vmin=0, vmax=1, cmap=cm.rainbow)
    cf = plt.pcolor(X, Y, Z,  vmin=0, vmax=700,  cmap=cm.rainbow)
    plt.colorbar(cf)
    #cont = plt.contour(X, Y, Z, levels=np.arange(1,10))
    #cont.clabel(fmt='%2.0f', fontsize=14)
    #plt.show()
    plt.savefig("tst.png")
    #print anisos1.shape
    #print np.average(anisos1,axis=0)






run()
