#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import h5py
import ctypes

plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams["xtick.top"] = True
#plt.rcParams["xtick.labeltop"] = True


def get3ddata(f):
    data = np.genfromtxt(f, dtype=float)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    intensity = data[:, 4]
    error = data[:, 5]  # data[:, 5] are mask if the values are -1.
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
    #karr = np.zeros((nx, ny, nz))
    earr = np.ones((nx, ny, nz))*(-1.0)
    karr2 = np.zeros((nx, ny, nz))

    for _x, _y, _z, _intensity, _error  in zip(x, y, z, intensity, error):
        xx = np.where(abs(xlin - _x) < 0.0000001)
        yy = np.where(abs(ylin - _y) < 0.0000001)
        zz = np.where(abs(zlin - _z) < 0.0000001)
        #karr[xx, yy, zz] = _intensity + 0.00000001
        karr2[xx, yy, zz] = _intensity
        earr[xx, yy, zz] = _error

    #condition = karr > 0.0000000001  
    condition = earr >= 0.0  
    #karrnoval = np.extract(condition, karr)
    #print(karr2.size)
    #print(karrnoval.shape)

    #for _x, _y, _z, _intensity, _error  in zip(x, y, z, intensity, error):
    #    xx = np.where(abs(xlin - _x) < 0.0000001)
    #    yy = np.where(abs(ylin - _y) < 0.0000001)
    #    zz = np.where(abs(zlin - _z) < 0.0000001)
    #    karr2[xx, yy, zz] = _intensity
    #    earr2[xx, yy, zz] = _error

    return karr2, condition


def make_mappable(maxvalue):
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=maxvalue)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("jet")
    mappable = ScalarMappable(norm=norm, cmap=cmap)
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


def unite_hdf5(num_txtfiles, head):
    for i in range(0, num_txtfiles):
        outfile = head + "out_hw_" + str(i) + "_" + str(i+1) + ".hdf5"
        f = h5py.File(outfile, 'r')
        if i == 0:
            tdata = np.zeros((f["data3"].shape[0], f["data3"].shape[1], f["data3"].shape[2], num_txtfiles))
            tcondition = np.zeros((f["data3"].shape[0], f["data3"].shape[1], f["data3"].shape[2], num_txtfiles))
        tdata[:, :, :, i] = f["data3"]
        tcondition[:, :, :, i] = f["condition"]
    outfile = head + "out_hw_all.hdf5"
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('data4', data=tdata)
        hf.create_dataset('condition', data=tcondition)


def save_eliminated_data_hdf5(head, data, condition):
    outfile = head + "eliminated_data.hdf5"
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('data4', data=data)
        hf.create_dataset('condition', data=condition)


def plot_crosssection(xi, xe, yi, ye, zi, ze, data4):
    fig = plt.figure(figsize=(18, 9))
    fig.suptitle("crosssections of 4D INS data #17714", fontsize="x-large")
    axindx=0
    for y in yi, ye:
        for z in zi, ze:
            axindx += 1
            ax = fig.add_subplot(3, 4, axindx)
            ##slice in one-step widths
            ax.pcolor(np.transpose(data4[:, y, z, :]), vmax=np.max(data4[:, y, z, :])/4, cmap='jet')
            ##slice in optimial bin widths
            #data = np.sum(np.sum(data4[:, y-1:y+1, z-6:z+6,:],axis=1),axis=1)
            #ax.pcolor(np.transpose(data), vmax=np.max(data)/4, cmap='jet')
            ax.text(2, 38, 'qy='+str(y)+', qz='+str(z), color='white')
            ax.set_xlabel('qx')
            ax.set_ylabel('E')
            ax.axvline(x=xi, color='white', lw=0.5)
            ax.axvline(x=xe, color='white', lw=0.5)
            ax.axhline(y=35, color='white', lw=0.5)
            ax.xaxis.set_label_coords(0.5,1.145)
            ax.tick_params(direction="in", color="white", top=True, labeltop=True, labelbottom=False)
            ax.axis('tight')
    for x in xi, xe:
        for z in zi, ze:
            axindx += 1
            ax = fig.add_subplot(3, 4, axindx)
            ##slice in one-step widths
            ax.pcolor(np.transpose(data4[x, :, z, :]), vmax=np.max(data4[x, :, z, :])/4, cmap='jet')
            ##slice in optimial bin widths
            #data = np.sum(np.sum(data4[x-1:x+1, :, z-6:z+6,:],axis=2),axis=0)
            #ax.pcolor(np.transpose(data), vmax=np.max(data)/4, cmap='jet')
            ax.text(2, 38, 'qx='+str(x)+', qz='+str(z), color='white')
            ax.set_xlabel('qy')
            ax.set_ylabel('E')
            ax.axvline(x=yi, color='white', lw=0.5)
            ax.axvline(x=ye, color='white', lw=0.5)
            ax.axhline(y=35, color='white', lw=0.5)
            ax.xaxis.set_label_coords(0.5,1.145)
            ax.tick_params(direction="in", color="white", top=True, labeltop=True, labelbottom=False)
            ax.axis('tight')
    for x in xi, xe:
        for y in yi, ye:
            axindx += 1
            ax = fig.add_subplot(3, 4, axindx)
            ##slice in one-step widths
            ax.pcolor(np.transpose(data4[x, y, :, :]), vmax=np.max(data4[x, y, :, :])/4, cmap='jet')
            ##slice in optimial bin widths
            #data = np.sum(np.sum(data4[x-1:x+1, y-1:y+1, :, :],axis=1),axis=0)
            #ax.pcolor(np.transpose(data), vmax=np.max(data)/4, cmap='jet')
            ax.text(2, 38, 'qx='+str(x)+', qy='+str(y), color='white')
            ax.set_xlabel('qz')
            ax.set_ylabel('E')
            ax.axvline(x=zi, color='white', lw=0.5)
            ax.axvline(x=ze, color='white', lw=0.5)
            ax.axhline(y=35, color='white', lw=0.5)
            ax.xaxis.set_label_coords(0.5,1.145)
            ax.tick_params(direction="in", color="white", top=True, labeltop=True, labelbottom=False)
            ax.axis('tight')
    fig.subplots_adjust(top=0.90)

def run():
    num_txtfiles = 41
    #head = "/home/kazu/desktop/200120/"
    head = "./"
    #gen_hdf5(num_txtfiles, head)
    #unite_hdf5(num_txtfiles, head)

    outfile = head + "out_hw_all.hdf5"

    f = h5py.File(outfile)
    data4 = f["data4"]
    condition = f["condition"]
    xi = 60
    xe = 84
    yi = 37
    ye = 68
    zi = 8
    ze = 27
    ei = 0
    ee = 35
    #plot_crosssection(xi, xe, yi, ye, zi, ze, data4)
    #plot_crosssection(xi, xe, yi, ye, zi, ze, condition)

    #save_eliminated_data_hdf5(head, data4[xi:xe, yi:ye, zi:ze, ei:ee], condition[xi:xe, yi:ye, zi:ze, ei:ee])
    save_eliminated_data_hdf5(head, data4[:, :, :, 5:], condition[:, :, :, 5:])


def run_tst():
    head = "/home/kazu/desktop/200120/"
    txtfile = head + "out_hw_0_1.txt"
    data, condition = get3ddata(txtfile)
    #plt.figure(figsize=(16, 8))
    #data2 = np.sum(data, axis=0)
    #plt.pcolor(np.transpose(data2), vmax=np.max(data2), cmap='jet')
    #mappable = make_mappable(np.max(data2))
    #plt.colorbar(mappable)


run()
#run_tst()
#plt.show()
#plt.savefig("crosssection.eps")
