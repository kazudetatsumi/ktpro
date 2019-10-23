#!/usr/bin/env python

import h5py
import numpy as np
#import matplotlib.pyplot as plt
import yaml
import math, cmath
#import matplotlib.cm as cm
from scipy.interpolate import griddata
from scipy import stats

#plt.style.use('classic')
#plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams['font.size'] = 12
#plt.rcParams['axes.linewidth'] = 1
#plt.rcParams['ytick.major.width'] = 1
#plt.rcParams['ytick.major.size'] = 5
#fig, ax = plt.subplots(1,2,  figsize=(14, 14))


def parse_band(bfile):
    f = h5py.File(bfile)
    xdata = f["distance"]
    ydata = f["frequency"]
    zdata = f["eigenvector"]
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
    x = np.squeeze(x) #(gp)
    y = np.squeeze(y) #(gp, band)
    print "Y", y.shape
    ylin = np.linspace(-1.0, 40.0, nb+1)
    nsize = y.shape
    nc = np.zeros((nsize[0], nb)) #(gp, nb)
    #values = np.ones(nsize[1])
    for i in range(0, nsize[0]):
        nc[i, :], res2, res3 = stats.binned_statistic(y[i, :], k[i,:], "sum", ylin)

    X, Y = np.meshgrid(x, ylin)
    nc = np.transpose(nc)  #(nb, gp)
    for i in range(0, nsize[0]):
        print "nc_CHK", np.sum(nc[:, i])
    #plt.figure()
    #plt.pcolor(X, Y, nc, cmap=cm.gray)
    return X, Y, nc


def oned_kappa(qp, kappa, omega):
    qpx = qp[:, 0]
    qpy = qp[:, 1]
    qpz = qp[:, 2]
    nsize = kappa.shape
    onedk = np.zeros((9, nsize[1]))
    for i in range(0, 9):
        condition = abs(qpx - (4.0/14.0)) + abs(qpy - (1.0/14.0)) + abs(qpz - (i/16.0)) < 0.01
        print qp[condition, :]
        _k = kappa[condition, :]
        _o = omega[condition, :]
        _o = np.squeeze(_o)
        onedk[i, :] = _k
        for j in range(0, nsize[1]-1):
            if _o[j] > _o[j+1]:
                print "order error in omega values"
    return onedk


def tst(ic, celldata, bondlenlim, ydata, zdata):
    nic = len(ic)
    latvec = np.array(celldata["lattice"])
    numa = len(celldata["points"])
    frac = np.zeros((numa, 3))
    for i in range(0, numa):
        frac[i, :] = np.array(celldata["points"][i]["coordinates"])
    t = np.zeros((27, 3))
    t1 = np.tile(np.array([1]), 9)
    t2 = np.tile(np.array([0]), 9)
    t3 = np.tile(np.array([-1]), 9)
    t[:, 0] = np.r_[t1, t2, t3]
    t[:, 1] = np.tile(np.array([1, 1, 1, 0, 0, 0, -1, -1, -1]), 3)
    t[:, 2] = np.tile(np.array([1, 0, -1]), 9)
    _t = np.tensordot(np.ones((numa, nic)), t, axes=0)  # (numa, ic, 27, 3)
    __t = np.transpose(_t, (2, 0, 1, 3)) # (27, numa, ic, 3)
    f = np.tensordot(np.ones((27, nic)), frac, axes=0) # (27, ic, numa, 3)
    _f = np.transpose(f, (0, 2, 1, 3)) # (27, numa, ic, 3)
    __f = __t + _f   #(27, numa, ic, 3)
    cars = np.matmul(__f, latvec) #(27, numa, ic, 3)

    ics = frac[ic, :] #(ic, 3)
    _ics = np.tensordot(np.ones((27, numa)), ics, axes=0) # (27, numa, ic, 3)
    carics = np.matmul(_ics, latvec) #(27, numa, ic, 3)
    norm = np.linalg.norm(carics - cars, axis=3) # (27, numa, ic)
    norm = np.transpose(norm, (2, 0, 1))         # (ic, 27, numa)
    nn = np.where((norm < bondlenlim) & (norm > 0.0001 )) # tuple specifiyng (centeratom, trans indices, atom indices)
    for z in ic:
       tmpd = np.where(np.array(nn[0]) == z)
       if len(tmpd[0]) != 4: print "Warning, coordination number is not 4"
       crd = len(tmpd[0])


    # hard coding for generating q by hand
    q = np.zeros((9,3)); q[:, 0] = 4.0/14.0; q[:, 1] = 1.0/14.0; q[:, 2] = np.arange(0,9)/16.0

    # vectorize variables with respect to  # (coord, gp, 3,  band, theta)
    #rm = cars[13, ic, 0, :]  # (3)
    rm = carics[0, 0, :, :] # (ic, 3)
    rk = cars[nn[1], nn[2], nn[0], :].reshape(nic, crd, 3) # (ic, coord, 3)

    ni = np.array(nn[2]).reshape(nic, crd) # (ic, coord)
    #print ni[2,:]
    #print nn[2]

    theta = np.tensordot(np.ones((nic, crd, 9, 3, numa*3)), np.linspace(0, 2*math.pi, 100), axes=0) # (ic, coord, gp, 3, band, theta)

    zdata = np.squeeze(zdata)
    Em = np.zeros((nic, 9, 3, numa*3), dtype=complex)  # (ic, gp, 3, band)
    i = 0
    for ii in ic:
        Em[i, :, :, :] = zdata[:, ii*3:(ii+1)*3, :]   #
        i += 1
    _Em = np.tensordot(Em, np.ones((crd, 100)), axes=0) # (ic, gp, 3, band, coord,  theta)
    __Em = np.transpose(_Em, (0, 4, 1, 2, 3, 5)) # (ic, coord, gp, 3, band,   theta)
    

    Ek = np.zeros((nic, crd, 9, 3, numa*3), dtype=complex) # (ic, coord, gp, 3, band)
    j = 0
    for jj in ic:
        for i in range(0, crd):
            Ek[j, i, :, :, :] = zdata[:, ni[j, i]*3:(ni[j, i]+1)*3, :] 
        j += 1
    __Ek = np.tensordot(Ek, np.ones((100)), axes=0) # (ic, coord, gp, 3, band,  theta)

    argm = np.tensordot(q, np.transpose(rm), axes=1) # (gp, ic)
    _argm = np.tensordot(argm, np.ones((crd, 3, numa*3, 100)), axes=0) # (gp, ic, coord, 3, band, theta)
    __argm = np.transpose(_argm, (1, 2, 0, 3, 4, 5)) # (ic, coord, gp, 3, band, theta)

    argk = np.tensordot(q, np.transpose(rk, (2, 0, 1)), axes=1) # (gp, ic, coord)
    _argk = np.tensordot(argk, np.ones((3, numa*3, 100)), axes=0) # (gp, ic, coord, 3, band, theta)
    __argk = np.transpose(_argk, (1, 2, 0, 3, 4, 5)) # (ic, coord, gp, 3, band, theta)


    expm = np.exp((__argm - theta)*1j)
    expk = np.exp((__argk - theta)*1j)
    em = np.real(__Em*expm) # (ic, coord, gp, 3, band, theta)
    ek = np.real(__Ek*expk) # (ic, coord, gp, 3, band, theta)

    _rm = np.tensordot(rm, np.ones((crd, 9, numa*3, 100)), axes=0) # (ic, 3, coord, gp, band, theta)
    __rm = np.transpose(_rm, (0, 2, 3, 1, 4, 5)) # (ic, coord, gp, 3, band, theta) 
    _rk = np.tensordot(rk, np.ones((9, numa*3, 100)), axes=0) # (ic, coord, 3, gp, band, theta)
    __rk = np.transpose(_rk, (0, 1, 3, 2, 4, 5)) # (ic, coord, gp, 3, band, theta)
    stmk = np.abs(np.linalg.norm(__rm - __rk + em - ek, axis=3) - np.linalg.norm(__rm - __rk, axis=3))/np.linalg.norm(__rm - __rk, axis=3) # (ic, coord, gp, band, theta)
    avestmk = np.sum(stmk, axis=(0, 1, 4))/(nic*crd*math.pi) # (gp, band)
    #for i in range(0,9):
    #    plt.plot(stmk[0,0,i,0,:], label=i)
    #plt.legend()
    return avestmk




def tst2(ic, trans, celldata, bondlenlim, ydata, zdata):
    latvec = np.array(celldata["lattice"])
    numa = len(celldata["points"])
    frac = np.zeros((numa, 3))
    for i in range(0, numa):
        frac[i, :] = np.array(celldata["points"][i]["coordinates"])
    t = np.zeros((27, 3))
    t1 = np.tile(np.array([1]), 9)
    t2 = np.tile(np.array([0]), 9)
    t3 = np.tile(np.array([-1]), 9)
    t[:, 0] = np.r_[t1, t2, t3]
    t[:, 1] = np.tile(np.array([1, 1, 1, 0, 0, 0, -1, -1, -1]), 3)
    t[:, 2] = np.tile(np.array([1, 0, -1]), 9)
    _t = np.tensordot(np.ones((numa)), t, axes=0)  # (numa, 27, 3)
    __t = np.transpose(_t, (1, 0, 2)) # (27, numa, 3)
    f = np.tensordot(np.ones((27)), frac, axes=0) # (27, numa, 3)
    _f = __t + f   #(27, numa, 3)
    cars = np.matmul(_f, latvec) #(27, numa, 3)

    ics = frac[ic, :] #(3)
    _ics = np.tensordot(np.ones((27, numa)), ics, axes=0) # (27, numa,  3)
    carics = np.matmul(_ics, latvec) #(27, numa, 3)
    norm = np.linalg.norm(carics - cars, axis=2) # (27, numa)
    nn = np.where((norm < bondlenlim) & (norm > 0.0001 )) # tuple specifiyng (trans indices, atom indices)
    crd = len(nn[0])


    # hard coding for generating q by hand
    q = np.zeros((9,3)); q[:, 0] = 4.0/14.0; q[:, 1] = 1.0/14.0; q[:, 2] = np.arange(0,9)/16.0

    # vectorize variables with respect to  # (coord, gp, 3,  band, theta)
    rm = carics[0, 0, :] + np.matmul(trans, latvec) # (3)
    rk = cars[nn[0], nn[1], :] + np.matmul(trans, latvec) # (coord, 3)


    theta = np.tensordot(np.ones((crd, 9, 3, numa*3)), np.linspace(0, 2*math.pi, 100), axes=0) # (coord, gp, 3, band, theta)

    zdata = np.squeeze(zdata)
    Em = zdata[:, ic*3:(ic+1)*3, :]   # (gp, 3, band)
    _Em = np.tensordot(Em, np.ones((crd, 100)), axes=0) # (gp, 3, band, coord, theta)
    __Em = np.transpose(_Em, (3, 0, 1, 2, 4)) # (coord, gp, 3, band, theta)
    

    Ek = np.zeros((crd, 9, 3, numa*3), dtype=complex) # (coord, gp, 3, band)
    for i in range(0, crd):
        Ek[i, :, :, :] = zdata[:, nn[1][i]*3:(nn[1][i]+1)*3, :]  # (coord, gp, 3, band)
    __Ek = np.tensordot(Ek, np.ones((100)), axes=0) # (coord, gp, 3, band,  theta)

    argm = np.tensordot(q, rm, axes=1) # (gp)
    _argm = np.tensordot(argm, np.ones((crd, 3, numa*3, 100)), axes=0) # (gp, coord, 3, band, theta)
    __argm = np.transpose(_argm, (1, 0, 2, 3, 4)) # (coord, gp, 3, band, theta)

    argk = np.tensordot(q, np.transpose(rk), axes=1) # (gp, coord)
    _argk = np.tensordot(argk, np.ones((3, numa*3, 100)), axes=0) # (gp, coord, 3, band, theta)
    __argk = np.transpose(_argk, (1, 0, 2, 3, 4)) # (coord, gp, 3, band, theta)


    expm = np.exp((__argm - theta)*1j)
    expk = np.exp((__argk - theta)*1j)
    em = np.real(__Em*expm) # (coord, gp, 3, band, theta)
    ek = np.real(__Ek*expk) # (coord, gp, 3, band, theta)

    _rm = np.tensordot(rm, np.ones((crd, 9, numa*3, 100)), axes=0) # (3, coord, gp, band, theta)
    __rm = np.transpose(_rm, (1, 2, 0, 3, 4)) # (coord, gp, 3, band, theta) 
    _rk = np.tensordot(rk, np.ones((9, numa*3, 100)), axes=0) # (coord, 3, gp, band, theta)
    __rk = np.transpose(_rk, (0, 2, 1, 3, 4)) # (coord, gp, 3, band, theta)
    stmk = np.abs(np.linalg.norm(__rm - __rk + em - ek, axis=2) - np.linalg.norm(__rm - __rk, axis=2))/np.linalg.norm(__rm - __rk, axis=2) # (coord, gp, band, theta)
    avestmk = np.sum(stmk, axis=(0, 3))/(crd*math.pi) # (gp, band)
    #for i in range(0,9):
    #    plt.plot(stmk[0,i,0,:], label=i)
    #plt.legend()
    return avestmk





def run():
    bondlenlim = 2.0
    nybin = 90
    temp = 300

    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band_4-14_1-14.hdf5"
    cxdata, cydata, czdata = parse_band(cbfile)
    cpfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/primitive.yaml"
    ccelldata = parse_cell(cpfile)
    ckfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/noiso/kappa-m141416.bz.hdf5"

    numa = len(ccelldata["points"])
    print numa


    trans = np.zeros((7*14*16, 3))
    l = 0
    for i in range(0, 7):
        for j in range(0, 14):
            for k in range(0, 16):
                trans[l, :] = np.array([i, j, k])
                l += 1

    avestmk = np.zeros((9, numa*3))
    for l in range(0, 7*14*16):
    #    print j
    #    import datetime
    #    d = datetime.datetime.today()
    #    print('d:', d)
         print l, avestmk[0, 0]
         for i in range(0,12):
             avestmk += tst2(i, trans[l, :], ccelldata, bondlenlim, cydata, czdata) / (12*7*14*16)
    
    with h5py.File('asi3n4_avestmk.hdf5', 'w') as hf:
        hf.create_dataset('avestmk', data=avestmk)



run()
