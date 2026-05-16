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
    dx = 0.05
    dy = 0.05
    dz = 0.05
    #xlin = np.arange(min(x), max(x)+dx, dx)
    xlin = np.arange(-1.65, 4.1, dx)
    nx = xlin.shape[0]
    #ylin = np.arange(min(y), max(y)+dy, dy)
    ylin = np.arange(-2.1, 2.8, dy)
    ny = ylin.shape[0] #zlin = np.arange(min(z), max(z)+dz, dz)
    zlin = np.arange(-0.85, 0.9, dz)
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


def gen_hdf5(num_txtfiles, head):
    for i in range(0, num_txtfiles):
        txtfile = head + "out_hw_" + str(i) + "_" + str(i+1) + ".txt"
        outfile = head + "out_hw_" + str(i) + "_" + str(i+1) + ".hdf5"
        data3, condition = get3ddata(txtfile)
        with h5py.File(outfile, 'w') as hf:
            hf.create_dataset('data3', data=data3)
            hf.create_dataset('condition', data=condition)


def unit_hdf5(num_txtfiles, head):
    for i in range(0, num_txtfiles):
        outfile = head + "out_hw_" + str(i) + "_" + str(i+1) + ".hdf5"
        f = h5py.File(outfile, 'r')
        if i == 0:
            tdata = np.zeros((f["data3"].shape[0], f["data3"].shape[1], f["data3"].shape[2], num_txtfiles))
            tcondition = np.zeros((f["data3"].shape[0], f["data3"].shape[1], f["data3"].shape[2], num_txtfiles))
        tdata[:,:,:,i] = f["data3"]
        tcondition[:,:,:,i] = f["condition"]
    outfile = head + "out_hw_all.hdf5"
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('data4', data=tdata)
        hf.create_dataset('condition', data=tcondition)


def run():
    num_txtfiles = 41
    head = "/home/kazu/desktop/200109/fine/fine/"
    #gen_hdf5(num_txtfiles, head)
    #unit_hdf5(num_txtfiles, head)

    outfile = head + "out_hw_all.hdf5"

    f=h5py.File(outfile)
    data4 = f["data4"]
    plt.figure(figsize=(16, 8))
    #data3 = np.sum(data4[:,59:61,:,:], axis=1)
    #data2 = np.sum(data3[:,14:16,:], axis=1)
    print data4.shape
    #data2 = data4[:,60,15,:]
    data2 = data4[:,60,15,:]
    #print np.unravel_index(np.argmax(data4), data4.shape)
    plt.pcolor(np.transpose(data2), vmax=np.max(data2/50), cmap='jet')
    mappable = make_mappable(np.max(data2/50))
    plt.colorbar(mappable)
    #plt.plot(np.log(data4[33,42,:,0]))

    

   


def run_tst():
    head = "/home/kazu/desktop/200109/fine/fine/"
    txtfile = head + "out_hw_1_2.txt"
    data, condition = get3ddata(txtfile)
    plt.figure(figsize=(16, 8))
    data2 = np.sum(data, axis=0)
    plt.pcolor(np.transpose(data2), vmax=np.max(data2), cmap='jet')
    mappable = make_mappable(np.max(data2))
    plt.colorbar(mappable)


run()
#run_tst()
plt.show()
