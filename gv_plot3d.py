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




def getXYZD(qx, qy, qz, ani):
    #x = qx*0.1306
    #y = qz*0.342
    xlin = np.linspace(min(qx), max(qx), 10)
    ylin = np.linspace(min(qy), max(qy), 10)
    zlin = np.linspace(min(qz), max(qz), 20)
    X, Y, Z = np.meshgrid(xlin, ylin, zlin)
    D = griddata((qx, qy, qz), ani, (X, Y, Z), method='linear')

    return X, Y, Z, D


def run():
    cdir = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift"
    sdir = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift"
    max_freq = 15
    y_s = 0.0
    c = cdir + "/mesh.hdf5"
    s = sdir + "/mesh.hdf5"

    qp, aniso_ave, omega = parse_mesh(c, max_freq)
    print qp[:,0].shape
    print qp[:,1].shape
    print qp[:,2].shape
    print aniso_ave.shape
    X, Y, Z, D = getXYZD(qp[:,0], qp[:,1], qp[:,2], aniso_ave)
    print X.shape
    print Y.shape
    print Z.shape
    print D.shape
    cc = 0
    tmp = ""
    for xx in range(D.shape[0]):
        for yy in range(D.shape[1]):
            for zz in range(D.shape[2]):
                tmp = tmp+" "+str(D[xx, yy, zz])
                cc += 1
                if cc % 5 == 0:
                    cc = 0
                    print tmp
                    tmp = ""
    D2 = D.reshape(400,5)
    np.savetxt("array.txt", D2, fmt="%f")


    #print Z.shape
    #fig = plt.figure(figsize=(12, 12))
    #cf = plt.pcolor(X, Y, Z, vmin=0, vmax=10, cmap=cm.rainbow)
    #cf = plt.pcolor(X, Y, Z, vmin=0, vmax=1, cmap=cm.rainbow)
    #cf = plt.pcolor(X, Y, Z,  vmin=0, vmax=700,  cmap=cm.rainbow)
    #plt.colorbar(cf)
    #cont = plt.contour(X, Y, Z, levels=np.arange(1,10))
    #cont.clabel(fmt='%2.0f', fontsize=14)
    #plt.show()
    #plt.savefig("tst.png")
    #print anisos1.shape
    #print np.average(anisos1,axis=0)






run()
