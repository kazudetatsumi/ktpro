#!/usr/bin/env python

import h5py
import numpy as np
import yaml
import math, cmath
import sys
from scipy.interpolate import griddata
from scipy import stats



def parse_band(bfile, nb, nrb):
    f = h5py.File(bfile)
    zdata = f["eigenvector"]
    zdata = np.squeeze(zdata)
    zdatasize = zdata.shape
    if zdatasize[0] != nb:
        print "number of bands is not correctly considered!"
        sys.exit()
    zdata2 = np.zeros((nrb, zdatasize[1], zdatasize[2]))
    j = 0
    for i in range(0, nb):
        if i % 10 == 0 or i % 10 == 9 or i % 10 == 1:
            zdata2[j, :, :] = zdata[i, :, :]
            j += 1
    return(zdata2)


def parse_cell(pfile):
    with open(pfile) as f:
        data = yaml.load(f)
        celldata = data["primitive_cell"]
    return celldata


def tst2(ic, trans, celldata, bondlenlim, zdata, nb, nrb):
    latvec = np.array(celldata["lattice"])
    numa = len(celldata["points"])
    frac = np.zeros((numa - 12, 3))
    for i in range(12, numa):
        frac[i - 12, :] = np.array(celldata["points"][i]["coordinates"])
    t = np.zeros((27, 3))
    t1 = np.tile(np.array([1]), 9)
    t2 = np.tile(np.array([0]), 9)
    t3 = np.tile(np.array([-1]), 9)
    t[:, 0] = np.r_[t1, t2, t3]
    t[:, 1] = np.tile(np.array([1, 1, 1, 0, 0, 0, -1, -1, -1]), 3)
    t[:, 2] = np.tile(np.array([1, 0, -1]), 9)
    _t = np.tensordot(np.ones((numa - 12)), t, axes=0)  # (numa, 27, 3)
    __t = np.transpose(_t, (1, 0, 2)) # (27, numa, 3)
    f = np.tensordot(np.ones((27)), frac, axes=0) # (27, numa, 3)
    _f = __t + f   #(27, numa, 3)
    cars = np.matmul(_f, latvec) #(27, numa, 3)

    ics = frac[ic, :] #(3)
    _ics = np.tensordot(np.ones((27, numa - 12)), ics, axes=0) # (27, numa,  3)
    carics = np.matmul(_ics, latvec) #(27, numa, 3)
    norm = np.linalg.norm(carics - cars, axis=2) # (27, numa)
    nn = np.where((norm < bondlenlim) & (norm > 0.0001 )) # tuple specifiyng (trans indices, atom indices)
    crd = len(nn[0])


    # hard coding for generating q by hand
    q = np.zeros((nrb,3)); q[:, 0] = 4.0/14.0; q[:, 1] = 1.0/14.0; 
    j = 0
    for i in range(0, nb):
        if i % 10 == 0 or i % 10 == 9 or i % 10 == 1:
            q[j, 2] = i / 160.0
            j += 1

    # vectorize variables with respect to  # (coord, gp, 3,  band, theta)
    rm = carics[0, 0, :] + np.matmul(trans, latvec) # (3)
    rk = cars[nn[0], nn[1], :] + np.matmul(trans, latvec) # (coord, 3)


    theta = np.tensordot(np.ones((crd, nrb, 3, numa*3)), np.linspace(0, 2*math.pi, 100), axes=0) # (coord, gp, 3, band, theta)

    zdata = np.squeeze(zdata)
    Em = zdata[:, ic*3:(ic+1)*3, :]   # (gp, 3, band)
    _Em = np.tensordot(Em, np.ones((crd, 100)), axes=0) # (gp, 3, band, coord, theta)
    __Em = np.transpose(_Em, (3, 0, 1, 2, 4)) # (coord, gp, 3, band, theta)
    

    Ek = np.zeros((crd, nrb, 3, numa*3), dtype=complex) # (coord, gp, 3, band)
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

    _rm = np.tensordot(rm, np.ones((crd, nrb, numa*3, 100)), axes=0) # (3, coord, gp, band, theta)
    __rm = np.transpose(_rm, (1, 2, 0, 3, 4)) # (coord, gp, 3, band, theta) 
    _rk = np.tensordot(rk, np.ones((nrb, numa*3, 100)), axes=0) # (coord, 3, gp, band, theta)
    __rk = np.transpose(_rk, (0, 2, 1, 3, 4)) # (coord, gp, 3, band, theta)
    stmk = np.abs(np.linalg.norm(__rm - __rk + em - ek, axis=2) - np.linalg.norm(__rm - __rk, axis=2))/np.linalg.norm(__rm - __rk, axis=2) # (coord, gp, band, theta)
    avestmk = np.sum(stmk, axis=(0, 3))/(crd*math.pi) # (gp, band)
    #for i in range(0,9):
    #    plt.plot(stmk[0,i,0,:], label=i)
    #plt.legend()
    return avestmk


def run():
    bondlenlim = 3.0
    temp = 300
    nb = 81
    nrb = (nb - 11) / 10 * 3 + 2*2

    cbfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/band.hdf5"
    czdata = parse_band(cbfile, nb, nrb)
    cpfile = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/primitive.yaml"
    ccelldata = parse_cell(cpfile)

    numa = len(ccelldata["points"])
    print numa


    trans = np.zeros((7*14*160, 3))
    l = 0
    for i in range(0, 7):
        for j in range(0, 14):
            for k in range(0, 160):
                trans[l, :] = np.array([i, j, k])
                l += 1

    avestmk = np.zeros((nrb, numa*3))
    for l in range(0, 7*14*160/4):
    #for l in range(7*14*160/4, 7*14*160*2/4):
    #for l in range(7*14*160*2/4, 7*14*160*3/4):
    #for l in range(7*14*160*3/4, 7*14*160*4/4):
    #for l in range(0, 10):
        #print j
         import datetime
         d = datetime.datetime.today()
         print('d:', d)
         print l, avestmk[0, 0]
         for i in range(0, 16):
             avestmk += tst2(i, trans[l, :], ccelldata, bondlenlim, czdata, nb, nrb) / (16*7*14*160)
    
    with h5py.File('asi3n4_avestmk_nn_fine1.hdf5', 'w') as hf:
        hf.create_dataset('avestmk', data=avestmk)



run()
