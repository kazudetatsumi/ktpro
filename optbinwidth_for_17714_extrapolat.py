#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes

def get3ddata(f):
    data = np.genfromtxt(f, dtype=None)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    intensity = data[:, 4]
    dx = 0.02
    dy = 0.02
    dz = 0.02
    #xlin = np.arange(min(x), max(x)+dx, dx)
    xlin = np.arange(-1.64, 4.10, dx)
    nx = xlin.shape[0]
    #ylin = np.arange(min(y), max(y)+dy, dy)
    ylin = np.arange(-2.10, 2.78, dy)
    ny = ylin.shape[0]
    #zlin = np.arange(min(z), max(z)+dz, dz)
    zlin = np.arange(-0.84, 0.90, dz)
    nz = zlin.shape[0]
    karr = np.zeros((nx, ny, nz))
    karr2 = np.zeros((nx, ny, nz))

    for _x, _y, _z, _intensity in zip(x, y, z, intensity):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        zz = np.where(abs(zlin - _z) < 0.0000001)
        karr[xx, yy, zz] = _intensity + 0.00000001

    condition = karr > 0.0000000001
    karrnoval = np.extract(condition, karr)
    print karr2.size
    print karrnoval.shape

    for _x, _y, _z, _intensity in zip(x, y, z, intensity):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        zz = np.where(abs(zlin - _z) < 0.0000001)
        karr2[xx, yy, zz] = _intensity

    return karr2, condition


def make_mappable(maxvalue):
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=maxvalue)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("jet")
    mappable=ScalarMappable(norm=norm,cmap=cmap)
    mappable._A = []
    return mappable


def run():
    head = "/home/kazu/desktop/200108/"
    txtfile = head + "out_hw_0.1_0.2.txt.small"
    outfile = head + "data3.hdf5"
    data, condition = get3ddata(txtfile)
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('data3', data=data3)


def run_tst():
    head = "/home/kazu/desktop/200108/"
    txtfile = head + "out_hw_0.1_0.2.txt"
    data, condition = get3ddata(txtfile) 
    plt.figure(figsize=(16, 8)) 
    data2 = np.sum(data, axis=0)
    plt.pcolor(np.transpose(data2), vmax=np.max(data2), cmap='jet')
    mappable = make_mappable(np.max(data2))
    plt.colorbar(mappable)





run_tst()
plt.show()
