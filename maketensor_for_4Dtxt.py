#!/usr/bin/env python
import numpy as np
import h5py


def get4ddata(f):
    data = np.genfromtxt(f, dtype=float, delimiter=',')
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    w = data[:, 3]
    intensity = data[:, 4]  # data[:, 4] are mask if the values are 1e*30

    dx = 0.05
    dy = 0.05
    dz = 0.05
    dw = 0.20
    # xlin = np.arange(min(x), max(x)+dx, dx)
    xlin = np.arange(-0.65, 3.15, dx)
    nx = xlin.shape[0]
    # ylin = np.arange(min(y), max(y)+dy, dy)
    ylin = np.arange(-0.9, 4.45, dy)
    ny = ylin.shape[0]
    # zlin = np.arange(min(z), max(z)+dz, dz)
    zlin = np.arange(-0.8, 0.60, dz)
    nz = zlin.shape[0]
    # wlin = np.arange(min(w), max(w)+dw, dw)
    wlin = np.arange(-8.0, 36.4, dw)
    nw = wlin.shape[0]
    karr2 = np.zeros((nx, ny, nz, nw))

    for _x, _y, _z, _w, _intensity in zip(x, y, z, w,  intensity):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        zz = np.where(abs(zlin - _z) < 0.0000001)
        ww = np.where(abs(wlin - _w) < 0.0000001)
        karr2[xx, yy, zz, ww] = _intensity

    condition = karr2 < 0.9e+30

    return karr2, condition


def gen_hdf5(num_txtfiles, head):
    for i in range(0, num_txtfiles):
        txtfile = head + "Output4D_00_" + str((i+1)*60) + ".txt"
        outfile = head + "Output4D_00_" + str((i+1)*60) + ".hdf5"
        data4, condition = get4ddata(txtfile)
        with h5py.File(outfile, 'w') as hf:
            hf.create_dataset('data4', data=data4)
            hf.create_dataset('condition', data=condition)


def run():
    num_txtfiles = 1
    head = "./"
    gen_hdf5(num_txtfiles, head)


run()
