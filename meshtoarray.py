#!/usr/bin/env python
import numpy as np
import h5py
import yaml
from matplotlib import pyplot as plt


def get_data(fn):
    f = h5py.File(fn, 'r')
    ##q = f["qpoint"]      # (nqps, 3)
    ##o = f["frequency"] # (nqps, 3natom)
    # we patch negative frequencies with zero.
    ##o[o < 0] = 0
    ##mesh = f["mesh"] # (nqps, 3natom)
    ##x = q[:, 0]
    ##y = q[:, 1]
    ##z = q[:, 2]
    ##xlin = get_coord(mesh[0]) 
    ##ylin = get_coord(mesh[1]) 
    ##zlin = get_coord(mesh[2]) 
    ##nx = xlin.shape[0]
    ##ny = ylin.shape[0]
    ##nz = zlin.shape[0]
    ##no = o.shape[1]
    ##karr = np.zeros((nx, ny, nz, no))
    ##karr2 = np.zeros((nx, ny, nz, no))

    # here we check whether we correctly read the whole lines in the input file
    # and generate a matrix "condition" which describes whether the element location is included in the input file or not.

    ##i = 0

    ##for _x, _y, _z in zip(x, y, z):
    ##   xx =  np.where(abs(xlin - _x) < 0.0000001)
    ##   yy =  np.where(abs(ylin - _y) < 0.0000001)
    ##   zz =  np.where(abs(zlin - _z) < 0.0000001)
    ##   karr[xx, yy, zz, :] = o[i, :] + 0.00000001
    ##   i += 1
    ##   print i

    ##condition = karr > 0.0000000001
    ##print condition.shape
    ##karrnonzero = np.extract(condition, karr)
    ##ndata = x.shape[0]
    ##if karrnonzero.shape[0] != x.shape[0]*no:
    ##     print "num of nonzero karr is not num of data", karrnonzero.shape, x.shape[0]*no
    ##else:
    ##     print "num of nonzero karr matches  num of data", karrnonzero.shape, x.shape[0]*no

    #here we regenerate the data matrix purely.

    ##i = 0

    ##for _x, _y, _z in zip(x, y, z):
    ##   xx =  np.where(abs(xlin - _x) < 0.0000001)
    ##   yy =  np.where(abs(ylin - _y) < 0.0000001)
    ##   zz =  np.where(abs(zlin - _z) < 0.0000001)
    ##   karr2[xx, yy, zz, :] = o[i, :] + 0.00000001
    ##   i += 1
    mesh = f["mesh"]
    n = f["frequency"].shape[1]
    o = np.reshape(f["frequency"][:], [mesh[0], mesh[1], mesh[2], n], order='F')
    o[ o< 0] = 0.0

    return o
    #return karr2[xi:xf, yi:yf], condition[xi:xf, yi:yf]


def get_coord(n):
    if n%2 == 1:
        x = np.arange(-n/2.0 + 1, n-n/2.0)/(n*1.0) 
    if n%2 == 0:
        x = np.arange(-n/2.0 - 1/2.0 + 1, n-n/2.0 + 1/2.0)/(n*1.0) 
    print x.shape[0]
    return x


def parse_rlat(my):
    with open(my) as f:
        data = yaml.load(f)
        rlat = data["primitive_cell"]["reciprocal_lattice"]
    rlat = np.array(rlat)
    return rlat


def conv_lore(data, gamma, do):
    x = np.arange(np.min(data), np.max(data)+5+do, do)
    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]
    no = data.shape[3]
    noo = x.shape[0]
    data2 = np.zeros((nx, ny, nz, noo))
    for i in range(0, nx):
        print i
        for j in range(0, ny):
            for k in range(0, nz):
                for h in range(0, no):
                    data2[i, j, k, :] += fun_lore(x, gamma, data[i, j, k, h])
    data2 = data2/(np.sum(data2)*do)
    return data2, x


#def poisson_sample(data, n):
#    pdata = np.zeros_like(data)




def fun_lore(x, gamma, mean):
    dx = x[1] - x[0]
    l = (1 / 3.14) * 0.5 * gamma / ((x - mean)**2 + (0.5 * gamma)**2) 
    l = l/(np.sum(l)*dx)
    return l


def run():
    gamma = 0.2
    do = 0.02
    n = 10000000
    #py = "/home/kazu/cscl/phonopy_222/phonopy.yaml"
    #rlat = parse_rlat(py)
    #ra = (np.sum(rlat[0,:]*rlat[0,:]))**0.5
    #print ra

    fn = "/home/kazu/cscl/phonopy_222/m200200200/mesh.hdf5"
    #q, omega, mesh = get_data(fn)
    data = get_data(fn)
    data2, x = conv_lore(data, gamma, do)
    data3 = np.random.poisson(data2*n)
    print data3.shape
    

    #print x
    #print q[0:21,0]
    plt.figure(figsize=(16, 8))
    plt.pcolor(data2[:, 0, 0, :])
    #plt.scatter(x, data3[10, 5, 7, :])
    #plt.plot(x, data2[8, 11, 12, :])




run()
plt.show()
