#!/usr/bin/env python
# extract the reciplocal distribution of group velocities on a reciplocal plane.
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
    aniso = ((gv[:, :, 0])**2 + (gv[:, :, 1])**2)**0.5
    #aniso = gv[:, :, 2]**2
    #aniso = abs(gv[:, :, 2])
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

def selectdata(qp, aniso, xy):
    dim = qp.shape
    q = []
    a = []
    for i in range(0, dim[0]):
        if qp[i, 0] == xy[0] and qp[i, 1] == xy[1]:
            q += list(qp[i, 2])
            a += list(aniso[i])
    q1 = np.array(q)
    a1 = np.array(a)

    return q1, a1


def getXYZD(qx, qy, qz, ani, n):
    xlin = np.linspace(min(qx), max(qx), n[0])
    ylin = np.linspace(min(qy), max(qy), n[1])
    zlin = np.linspace(min(qz), max(qz), n[2])
    X, Y, Z = np.meshgrid(xlin, ylin, zlin)
    D = griddata((qx, qy, qz), ani, (X, Y, Z), method='linear')

    D2 = np.transpose(D, (2, 1, 0))             # The output format is adjusted to the CHGCAR file.
    D3 = D2.reshape(n[0]*n[1]*n[2] / 5, 5)      # The transposing is nescessary for it. 
    return X, Y, Z, D3


def run():
    cdir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift"
    sdir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift"
    ssdir = "/home/kazu/bsi3n4_m/phonopy_doubled_334"
    max_freq = 15
    xy = np.array([0, 0])
    c = cdir + "/mesh.hdf5"
    s = sdir + "/mesh.hdf5"
    ss = ssdir + "/mesh.hdf5"

    qp, aniso_ave, omega = parse_mesh(c, max_freq)
    qz, anisoz = selectdata(qp, aniso_ave, xy)






run()
