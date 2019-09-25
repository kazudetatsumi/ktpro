#!/usr/bin/env python
import numpy as np
import h5py
import yaml


def get_data(fn):
    f = h5py.File(fn)
    q = f["qpoint"]      # (nqps, 3)
    o = f["frequency"] # (nqps, 3natom)
    # we patch negative frequencies with zero.
    o[o < 0] = 0
    mesh = f["mesh"] # (nqps, 3natom)
    x = q[:, 0]
    y = q[:, 1]
    z = q[:, 2]
    xlin = get_coord(mesh[0]) 
    ylin = get_coord(mesh[1]) 
    zlin = get_coord(mesh[2]) 
    nx = xlin.shape[0]
    ny = ylin.shape[0]
    nz = zlin.shape[0]
    no = o.shape[1]
    karr = np.zeros((nx, ny, nz, no))
    karr2 = np.zeros((nx, ny, nz, no))

    # here we check whether we correctly read the whole lines in the input file
    # and generate a matrix "condition" which describes whether the element location is included in the input file or not.

    i = 0

    for _x, _y, _z in zip(x, y, z):
       xx =  np.where(abs(xlin - _x) < 0.0000001)
       yy =  np.where(abs(ylin - _y) < 0.0000001)
       zz =  np.where(abs(zlin - _z) < 0.0000001)
       karr[xx, yy, zz, :] = o[i, :] + 0.00000001
       i += 1

    condition = karr > 0.0000000001
    print condition.shape
    karrnonzero = np.extract(condition, karr)
    ndata = x.shape[0]
    if karrnonzero.shape[0] != x.shape[0]*no:
         print "num of nonzero karr is not num of data", karrnonzero.shape, x.shape[0]*no
    else:
         print "num of nonzero karr matches  num of data", karrnonzero.shape, x.shape[0]*no

    #here we regenerate the data matrix purely.

    i = 0

    for _x, _y, _z in zip(x, y, z):
       xx =  np.where(abs(xlin - _x) < 0.0000001)
       yy =  np.where(abs(ylin - _y) < 0.0000001)
       zz =  np.where(abs(zlin - _z) < 0.0000001)
       karr2[xx, yy, zz, :] = o[i, :] + 0.00000001
       i += 1
    return karr2
    #return karr2[xi:xf, yi:yf], condition[xi:xf, yi:yf]


def get_coord(n):
    x = np.arange(-n/2 + 1, n-n/2)/(n*1.0) 
    return x


def parse_rlat(my):
    with open(my) as f:
        data = yaml.load(f)
        rlat = data["primitive_cell"]["reciprocal_lattice"]
    rlat = np.array(rlat)
    return rlat


def conv_lore(data, gamma, do):
    x = np.arange(min(data), max(data)+do, do)
    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]
    no = data.shape[3]
    data2 = np.zeros((nx, ny, nz,
    for i in range(0, nx):
        for j in range(0, ny):
            for k in range(0, nz):
                for h in range(0, no):




def fun_lore(x, gamma, mean):
    dx = x[1] - x[0]
    l = (1 / 3.14) * 0.5 * gamma / ((x - mean)**2 + (0.5 * gamma)**2) 
    l = l/(np.sum(l)*dx)
    return l


def run():
    gamma = 0.2
    do = 0.01
    #py = "/home/kazu/cscl/phonopy_222/phonopy.yaml"
    #rlat = parse_rlat(py)
    #ra = (np.sum(rlat[0,:]*rlat[0,:]))**0.5
    #print ra

    fn = "/home/kazu/cscl/phonopy_222/mesh.hdf5"
    #q, omega, mesh = get_data(fn)
    data = get_data(fn)

    data2 = conv_lore(data, gamma, do)
    #print x
    #print q[0:21,0]




run()
