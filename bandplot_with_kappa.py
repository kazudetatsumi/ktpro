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
fig, ax = plt.subplots(1,2,  figsize=(14, 14))


def parse_band(bfile):
    f = h5py.File(bfile)
    xdata = f["distance"]
    ydata = f["frequency"]
    zdata = f["eigenvector"]
    return(xdata, ydata, zdata)



def parse_kappa(filename,  temp):
    f = h5py.File(filename, 'r')
    mesh = f["mesh"]
    print f.keys()
    kuc = f["kappa_unit_conversion"][()]
    qp = f["qpoint"][:, :]  # (gp, 3)
    print "qp",qp.shape
    if filename == "/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/kappa-m141416.bz.hdf5":
        print "chk", filename,"does not contain proper gv"
        gv, qpm = parse_mesh("/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/mesh_m141416.hdf5")
        qpm[qpm < 0] += 1
        print "qpm", qpm.shape
        print np.where((qp - qpm) > 0.01)
        print qpm[86,:], qp[86,:]
    else:
        gv = f["group_velocity"] #(gp, band)
    qp[qp < 0] += 1
    omega = f["frequency"][:, :] # (gp, band)
    gamma = f['gamma'][:] #KT, (temps, gp, band)
    cv = f['heat_capacity'][:] # (temps, gp, band)
    t = np.array(f['temperature'][:], dtype='double') # (temps)
    condition = abs(t - temp) < 0.001
    cv = np.squeeze(cv[condition, :, :])
    gamma = np.squeeze(gamma[condition, :, :])
    kappa = kuc * cv * gv[:, :, 2] * gv[:, :, 2] / (2 * gamma)  # (gp, band)
    kappa[0, 0:3] = 0
    return(qp, kappa, omega)




def get_files(bfile, pfile):
    xdata, ydata, zdata = parse_band(bfile)
    with open(pfile) as f:
        data = yaml.load(f)
        celldata = data["primitive_cell"]
    return xdata, ydata, zdata, celldata



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
    #plt.figure()
    #plt.pcolor(X, Y, nc, cmap=cm.gray)
    return X, Y, nc


def oned_kappa(qp, kappa, omega):
    qpx = qp[:,0]
    qpy = qp[:,1]
    qpz = qp[:,2]
    condition = abs(qpx - 4/14) + abs(qpy - 1/14) + abs(qpz - 1/8) < 0.01
    _k = kappa[condition,:]
    _o = omega[condition,:]
    print condition.shape
    print "selected kappa",_k, _k.shape
    print "selected omega",_o, _o.shape


def run():
    bondlenlim = 2.0
    nybin = 40
    temp = 300

    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_4-14_1-14.hdf5"
    cpfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/primitive.yaml"
    cxdata, cydata, czdata, ccelldata = get_files(cbfile, cpfile)
    ckfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/noiso/kappa-m141416.bz.hdf5"
    qp1, kappa1, omega1 = parse_kappa(ckfile, temp)
    kappa = oned_kappa(qp1, kappa1, omega1)


    cX,cY,cnc = binimage(cxdata, cydata, nybin)
    sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_4-14_1-14.hdf5"
    spfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/primitive.yaml"
    sxdata, sydata, szdata, scelldata = get_files(sbfile, spfile)
    sX,sY,snc = binimage(sxdata, sydata, nybin)
    im=ax[0].pcolor(cX, cY, cnc, cmap=cm.gray)
    im=ax[1].pcolor(sX, sY, snc, cmap=cm.gray)
    #from matplotlib.colors import Normalize
    #norm = Normalize(vmin=0, vmax=1)
    #from matplotlib.cm import ScalarMappable, get_cmap
    #cmap = get_cmap("jet")
    #mappable=ScalarMappable(norm=norm,cmap=cmap)
    #mappable._A = []
    #plt.colorbar(mappable)

#THztoKayser = 33.35641
run()
plt.show()
