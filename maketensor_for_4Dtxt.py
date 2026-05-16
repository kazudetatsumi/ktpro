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
    #wlin = np.arange(-8.0, -7.8, dw)
    nw = wlin.shape[0]
    karr2 = np.ones((nx, ny, nz, nw))*(-1)

    if nx*ny*nz*nw != intensity.shape[0]:
        print("number of elements is wrong!", nx*ny*nz*nw, intensity.shape[0])

    for _x, _y, _z, _w, _intensity in zip(x, y, z, w,  intensity):
        xx = np.where(abs(xlin - _x) < 0.0001)
        yy = np.where(abs(ylin - _y) < 0.0001)
        zz = np.where(abs(zlin - _z) < 0.0001)
        ww = np.where(abs(wlin - _w) < 0.0001)
        #print(xx[0][0], yy[0][0], zz[0][0], ww[0][0])
        if xx[0].shape[0] == 0 or yy[0].shape[0] == 0 or zz[0].shape[0] == 0  or ww[0].shape[0] == 0:
            print("missing data!!:",_x, _y, _z, _w, _intensity)
        else:
            karr2[xx, yy, zz, ww] = _intensity

    if abs(np.sum(karr2) - np.sum(intensity)) > 1:
        print("total intensity is not preserved!!", np.sum(karr2), np.sum(intensity))
    print(np.where(karr2 < 0))



    condition = karr2 < 0.9e+30
    karr2 = karr2*condition

    return karr2, condition


def get4ddata_dummy():

    xlin = np.arange(0, 7, 1)
    nx = xlin.shape[0]
    # ylin = np.arange(min(y), max(y)+dy, dy)
    ylin = np.arange(0, 6, 1)
    ny = ylin.shape[0]
    # zlin = np.arange(min(z), max(z)+dz, dz)
    zlin = np.arange(0, 5, 1)
    nz = zlin.shape[0]
    # wlin = np.arange(min(w), max(w)+dw, dw)
    wlin = np.arange(0, 4, 1)
    nw = wlin.shape[0]
    karr2 = np.zeros((nx, ny, nz, nw))
    intensity = np.zeros((nx*ny*nz*nw))
    x = np.zeros((nx*ny*nz*nw))
    y = np.zeros((nx*ny*nz*nw))
    z = np.zeros((nx*ny*nz*nw))
    w = np.zeros((nx*ny*nz*nw))

    i = 0
    for _w in wlin:
        for _x in xlin:
            for _y in ylin:
                for _z in zlin:
                    intensity[i] = _z*1000 + _y * 100 + _z * 10 + _w 
                    x[i] = _x
                    y[i] = _y
                    z[i] = _z
                    w[i] = _w
                    #if _w % 4 == 0 or _z % 2 == 0 :
                    #    intensity[i] = 1e+30
                    i = i + 1
    
                    

    for _x, _y, _z, _w, _intensity in zip(x, y, z, w,  intensity):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        zz = np.where(abs(zlin - _z) < 0.0000001)
        ww = np.where(abs(wlin - _w) < 0.0000001)
        karr2[xx, yy, zz, ww] = _intensity
    #karr3 = np.reshape(intensity, (nz, ny, nx, nw))
    karr3 = np.reshape(intensity, (nw, nx, ny, nz))
    karr3 = np.transpose(karr3, (1, 2, 3, 0))

    condition = karr2 < 0.9e+30
    karr2 = karr2*condition
    print(np.sum(karr3[:,:,:,:]-karr2[:,:,:,:]))
    return karr2, condition


def gen_hdf5(num_txtfiles, head):
    #for i in range(0, num_txtfiles):
    for i in range(8, 9):
        #txtfile = head + "Output4D_00_" + str((i+2)*60) + "_test.txt"
        txtfile = head + "Output4D_00_" + str((i+2)*60) + ".txt"
        outfile = head + "Output4D_00_" + str((i+2)*60) + ".hdf5"
        data4, condition = get4ddata(txtfile)
        #outfile = head + "Output4D_00_dummy.hdf5"
        #data4, condition = get4ddata_dummy()
        with h5py.File(outfile, 'w') as hf:
            hf.create_dataset('data4', data=data4)
            hf.create_dataset('condition', data=condition)


def run():
    num_txtfiles = 1
    head = "./"
    gen_hdf5(num_txtfiles, head)


run()
