
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
    xdata = np.squeeze(f["distance"])
    ydata = np.squeeze(f["frequency"])
    zdata = np.squeeze(f["eigenvector"])
    return(xdata, ydata, zdata)


def parse_cell(pfile):
    with open(pfile) as f:
        data = yaml.load(f)
        celldata = data["primitive_cell"]
    return celldata


def parse_mesh(filename):
    f = h5py.File(filename, 'r')
    gv = f["group_velocity"]
    qp = f["qpoint"][:, :]
    return gv, qp


def parse_kappa(filename,  temp):
    f = h5py.File(filename, 'r')
    mesh = f["mesh"]
    print f.keys()
    kuc = f["kappa_unit_conversion"][()]
    qp = f["qpoint"][:, :]  # (gp, 3)
    print "qp",qp.shape
    qp[qp < 0] += 1
    if filename == "/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/kappa-m141416.bz.hdf5":
        print "chk", filename,"does not contain proper gv"
        gv, qpm = parse_mesh("/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/mesh_m141416.hdf5")
        qpm[qpm < 0] += 1
        print "qpm", qpm.shape
        print "qpm2",np.where((qp - qpm) > 0.01)
        print "qpm3",qpm[186,:], qp[186,:]
    else:
        gv = f["group_velocity"] #(gp, band, 3)
    omega = f["frequency"][:, :] # (gp, band)
    gamma = f['gamma'][:] #KT, (temps, gp, band)
    cv = f['heat_capacity'][:] # (temps, gp, band)
    t = np.array(f['temperature'][:], dtype='double') # (temps)
    condition = abs(t - temp) < 0.001
    cv = np.squeeze(cv[condition, :, :])
    gamma = np.squeeze(gamma[condition, :, :])
    kappa = kuc * cv * gv[:, :, 2] * gv[:, :, 2] / (2 * gamma)  # (gp, band)
    kappa[0, 0:3] = 0
    print "kappa=",np.sum(kappa) / (mesh[0]*mesh[1]*mesh[2])
    return(qp, kappa, gv[:, :, 2], omega)





def binimage(x, y, k, nb):
    ylin = np.linspace(-1.0, 40.0, nb+1)
    nsize = y.shape
    nc = np.zeros((nsize[0], nb+1)) #(gp, nb)
    for i in range(0, nsize[0]):
        nc[i,0:nb], res2, res3 = stats.binned_statistic(y[i, :], k[i,:], "sum", ylin)

    nc = np.transpose(nc)  #(nb, gp)
    xsize = x.shape
    ysize = ylin.shape
    x = np.append(x, x[xsize[0]-1] + x[xsize[0]-1] - x[xsize[0]-2])
    ylin = np.append(ylin, ylin[ysize[0]-1] + ylin[ysize[0]-1] - ylin[ysize[0]-2])
    return nc, x, ylin



def parse_dgdq(gfile):
    f = h5py.File(gfile)
    avest = f["avestmk"]
    size = avest.shape
    ngp = size[0]
    nband = size[1]
    dgdq = np.zeros((ngp, nband))
    for i in range(1, ngp):
        dgdq[i, :] = avest[i, :] - avest[i-1, :]
    return dgdq


def parse_dgdq_fine(gfile):
    f = h5py.File(gfile)
    avest = f["avestmk"]
    size = avest.shape
    ngp = size[0]
    nband = size[1]
    n = ( size[0] + 2 ) / 3
    dgdq = np.zeros((n, nband))
    j = 0
    for i in range(1, ngp):
        if i == 1 or i == ngp - 1:
           dgdq[j, :] = np.abs(avest[i, :] - avest[i-1, :])
           j += 1
        elif i % 3 == 0:  
           dgdq[j, :] = np.abs((avest[i+1, :] - avest[i-1, :]) / 2.0)
           j += 1
    return dgdq


def oned_kappa(qp, kappa, omega):
    qpx = qp[:, 0]
    qpy = qp[:, 1]
    qpz = qp[:, 2]
    nsize = kappa.shape
    onedk = np.zeros((9, nsize[1]))
    for i in range(0, 9):
        condition = abs(qpx - (4.0/14.0)) + abs(qpy - (1.0/14.0)) + abs(qpz - (i/16.0)) < 0.01
        _k = kappa[condition, :]
        _o = omega[condition, :]
        _o = np.squeeze(_o)
        onedk[i, :] = _k
        for j in range(0, nsize[1]-1):
            if _o[j] > _o[j+1]:
                print "order error in omega values"
    return onedk




def run():
    bondlenlim = 2.0
    nybin = 30
    temp = 300

    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_4-14_1-14.hdf5"
    cxdata, cydata, czdata = parse_band(cbfile)
    cpfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/primitive.yaml"
    ccelldata = parse_cell(cpfile)
    ckfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/noiso/kappa-m141416.bz.hdf5"
    qp1, kappa1, gv1, omega1 = parse_kappa(ckfile, temp)
    kappa1= oned_kappa(qp1, kappa1, omega1)
    cgfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/asi3n4_avestmk_nn.hdf5"
    dgdq1 = parse_dgdq(cgfile)


    cnc, cx, cy = binimage(cxdata, cydata, kappa1, nybin)
    cnd, cx, cy = binimage(cxdata, cydata, dgdq1, nybin)



    sbfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/band_4-14_1-14.hdf5"
    sxdata, sydata, szdata = parse_band(sbfile)
    spfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/primitive.yaml"
    scelldata = parse_cell(spfile)
    skfile = "/home/kazu/bsi3n4_m/phono3py_doubled_112_fc2_334_sym_monk/kappa-m141416.bz.hdf5"
    qp2, kappa2, gv2, omega2 = parse_kappa(skfile, temp)
    kappa2 = oned_kappa(qp2, kappa2, omega2)
    gv2 = oned_kappa(qp2, gv2*gv2, omega2)
    sgfile = "/home/kazu/bsi3n4_m/phonopy_doubled_334/bsi3n4_avestmk_fine.hdf5"
    dgdq2 = parse_dgdq_fine(sgfile)

    snc, sx, sy = binimage(sxdata, sydata, kappa2, nybin)
    #snd, sx, sy = binimage(sxdata, sydata, gv2, nybin)
    snd, sx, sy = binimage(sxdata, sydata, dgdq2*dgdq2, nybin)
    snd[:,0]=0
    snd[:,8]=0
    

    maxsc=np.max(snc)
    maxscd=np.max(snd)
    im=ax[0].pcolor(sx, sy, snc, vmin=0, vmax=maxsc, cmap=cm.gray)
    im=ax[1].pcolor(sx, sy, snd, vmin=0, vmax=maxscd*0.1, cmap=cm.gray)
    ax[0].set_title('kappa beta ')
    ax[1].set_title('geometry beta ')
    plt.figure()
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=maxscd)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("gray")
    mappable=ScalarMappable(norm=norm,cmap=cmap)
    mappable._A = []
    plt.colorbar(mappable)

run()
plt.show()
