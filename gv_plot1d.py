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


def parse_mesh(filename,  max_freq, min_freq):
    f = h5py.File(filename, 'r')
    gv = f["group_velocity"]
    #aniso = gv[:, :, 2]**2 / (gv[:, :, 0]**2 + gv[:, :, 1]**2 + gv[:, :, 2]**2)
    #aniso = ((gv[:, :, 0])**2 + (gv[:, :, 1])**2)**0.5
    #aniso = gv[:, :, 2]**2
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
            if omega[i, j] < max_freq and omega[i, j] > min_freq:
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
    for qx, qy, qz, ani in zip(qp[:, 0], qp[:, 1], qp[:, 2], aniso):
        condition = qx == xy[0] and qy == xy[1]
        _qz = np.extract(condition, qz)
        _ani = np.extract(condition, ani)
        q += list(_qz)
        a += list(_ani)
    q1 = np.array(q)
    a1 = np.array(a)

    return q1, a1


def get_1d(f, max_freq, min_freq, xy):
    qp, aniso_ave, omega = parse_mesh(f, max_freq, min_freq)
    qz, anisoz = selectdata(qp, aniso_ave, xy)
    
    return(qz, anisoz)


def run():
    cdir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift"
    sdir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift"
    ssdir = "/home/kazu/bsi3n4_m/phonopy_doubled_334"
    min_freq = 10
    max_freq = 15
    xy = np.array([0, 0])
    c = cdir + "/mesh.hdf5"
    s = sdir + "/mesh.hdf5"
    ss = ssdir + "/mesh.hdf5"

    #qp, aniso_ave, omega = parse_mesh(ss, max_freq)
    #qz, anisoz = selectdata(qp, aniso_ave, xy)
    cqz, canisoz = get_1d(c, max_freq, min_freq, xy)
    ssqz, ssanisoz = get_1d(ss, max_freq, min_freq, xy)
    sqz, sanisoz = get_1d(s, max_freq, min_freq, xy)

    fig = plt.figure(figsize=(12, 12))
    plt.plot(cqz, canisoz)
    plt.plot(ssqz, ssanisoz)
    #plt.plot(sqz*2, sanisoz)
    plt.ylim(0, max(ssanisoz))
    plt.show()






run()
