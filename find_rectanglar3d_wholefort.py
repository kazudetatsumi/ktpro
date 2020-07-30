#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import ctypes
lib = ctypes.CDLL("/home/kazu/ktpro/rectanglar3d.so")


def calc_rectanglar_f90(uselight, ei, ee, av, ab, condition):
    class result(ctypes.Structure):
        _fields_ = [("lb_0", ctypes.c_int), ("ub_0", ctypes.c_int), ("lb_1",
                    ctypes.c_int), ("ub_1", ctypes.c_int), ("lb_2",
                    ctypes.c_int), ("ub_2", ctypes.c_int), ("lb_3",
                    ctypes.c_int), ("ub_3", ctypes.c_int)]
    lib.rectanglar.restype = result
    lib.rectanglar.argtypes = [ctypes.POINTER(ctypes.c_bool),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_double),
                               ctypes.POINTER(ctypes.c_double),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               ctypes.POINTER(ctypes.c_int),
                               np.ctypeslib.ndpointer(dtype=np.int32, ndim=4)]
    Nmax0 = condition.shape[0]
    Nmax1 = condition.shape[1]
    Nmax2 = condition.shape[2]
    Nmax3 = condition.shape[3]
    result = lib.rectanglar(
                        ctypes.byref(ctypes.c_bool(uselight)),
                        ctypes.byref(ctypes.c_int(ei)),
                        ctypes.byref(ctypes.c_int(ee)),
                        ctypes.byref(ctypes.c_double(av)),
                        ctypes.byref(ctypes.c_double(ab)),
                        ctypes.byref(ctypes.c_int(Nmax0)),
                        ctypes.byref(ctypes.c_int(Nmax1)),
                        ctypes.byref(ctypes.c_int(Nmax2)),
                        ctypes.byref(ctypes.c_int(Nmax3)),
                        condition)
    result_lb_0 = result.lb_0
    result_lb_1 = result.lb_1
    result_lb_2 = result.lb_2
    result_lb_3 = result.lb_3
    result_ub_0 = result.ub_0
    result_ub_1 = result.ub_1
    result_ub_2 = result.ub_2
    result_ub_3 = result.ub_3
    lb = np.array([result_lb_0, result_lb_1, result_lb_2, result_lb_3])
    ub = np.array([result_ub_0, result_ub_1, result_ub_2, result_ub_3])
    return lb, ub


def plot_crosssection(xi, xe, yi, ye, zi, ze, ei, ee,  data4):
    fig = plt.figure(figsize=(18, 9))
    fig.suptitle("crosssections of 4D INS data #17714", fontsize="x-large")
    axindx = 0
    for y in yi, ye:
        for z in zi, ze:
            axindx += 1
            ax = fig.add_subplot(3, 4, axindx)
            # slice in one-step widths
            ax.pcolor(np.transpose(data4[:, y, z, :]),
                      vmax=np.max(data4[:, y, z, :])/4, cmap='jet')
            # slice in optimial bin widths
            # data = np.sum(np.sum(data4[:, y-1:y+1, z-6:z+6,:],axis=1),axis=1)
            # ax.pcolor(np.transpose(data), vmax=np.max(data)/4, cmap='jet')
            ax.text(2, 38, 'qy='+str(y)+', qz='+str(z), color='white')
            ax.set_xlabel('qx')
            ax.set_ylabel('E')
            ax.axvline(x=xi, color='white', lw=0.5)
            ax.axvline(x=xe, color='white', lw=0.5)
            ax.axhline(y=ei, color='white', lw=0.5)
            ax.axhline(y=ee, color='white', lw=0.5)
            ax.xaxis.set_label_coords(0.5, 1.145)
            ax.axis('tight')
            ax.tick_params(direction="in", color="white",
                           top=True, labeltop=True, labelbottom=False)
    for x in xi, xe:
        for z in zi, ze:
            axindx += 1
            ax = fig.add_subplot(3, 4, axindx)
            # slice in one-step widths
            ax.pcolor(np.transpose(data4[x, :, z, :]),
                      vmax=np.max(data4[x, :, z, :])/4, cmap='jet')
            # slice in optimial bin widths
            # data = np.sum(np.sum(data4[x-1:x+1, :, z-6:z+6,:],axis=2),axis=0)
            # ax.pcolor(np.transpose(data), vmax=np.max(data)/4, cmap='jet')
            ax.text(2, 38, 'qx='+str(x)+', qz='+str(z), color='white')
            ax.set_xlabel('qy')
            ax.set_ylabel('E')
            ax.axvline(x=yi, color='white', lw=0.5)
            ax.axvline(x=ye, color='white', lw=0.5)
            ax.axhline(y=ei, color='white', lw=0.5)
            ax.axhline(y=ee, color='white', lw=0.5)
            ax.xaxis.set_label_coords(0.5, 1.145)
            ax.axis('tight')
            ax.tick_params(direction="in", color="white",
                           top=True, labeltop=True, labelbottom=False)
    for x in xi, xe:
        for y in yi, ye:
            axindx += 1
            ax = fig.add_subplot(3, 4, axindx)
            # slice in one-step widths
            ax.pcolor(np.transpose(data4[x, y, :, :]),
                      vmax=np.max(data4[x, y, :, :])/4, cmap='jet')
            # slice in optimial bin widths
            # data = np.sum(np.sum(data4[x-1:x+1, y-1:y+1, :, :],axis=1),axis=0)
            # ax.pcolor(np.transpose(data), vmax=np.max(data)/4, cmap='jet')
            ax.text(2, 38, 'qx='+str(x)+', qy='+str(y), color='white')
            ax.set_xlabel('qz')
            ax.set_ylabel('E')
            ax.axvline(x=zi, color='white', lw=0.5)
            ax.axvline(x=ze, color='white', lw=0.5)
            ax.axhline(y=ei, color='white', lw=0.5)
            ax.axhline(y=ee, color='white', lw=0.5)
            ax.xaxis.set_label_coords(0.5, 1.145)
            ax.axis('tight')
            ax.tick_params(direction="in", color="white",
                           top=True, labeltop=True, labelbottom=False)
    fig.subplots_adjust(top=0.90)


def run():
    # the lower and upper energy bondaries
    #ei = 81
    #ee = 207
    ei = 10
    ee = 70
    # parameters softening the conditions to select volume and boundaries
    # a_v = 0.9, ab = 0.999 
    av = 0.970
    av = 0.970
    ab = av
    # use argmaxvlight
    uselight = False 
    #maskfile = "/home/kazu/desktop/200204/coarse/hourbyhour/1h/out_hw_all.hdf5"
    maskfile = "/home/kazu/desktop/200204/fine/out_hw_all.hdf5"
    # maskfile = "/home/kazu/desktop/200204/fine/hourbyhour/1h/out_hw_all.hdf5"
    #maskfile = "/home/kazu/desktop/200522/Ei42/veryfineq/Output4D_00_840.hdf5"
    #maskfile = "/home/kazu/desktop/200522/Ei42/fineq/14m/Output4D_00_840.hdf5"
    # maskfile = "/home/kazu/desktop/200312/for_cscl/coarse/out_hw_all.hdf5"
    # maskfile = "/home/kazu/desktop/200312/for_cscl/out_hw_all.hdf5"
    print(maskfile)
    f = h5py.File(maskfile)
    condition = np.array(f["condition"][:, :, :, :], dtype=np.int32)
    lb, ub = calc_rectanglar_f90(uselight, ei, ee, av, ab, condition)
    print(lb)
    print(ub)
    xi = lb[0]
    xe = ub[0]
    yi = lb[1]
    ye = ub[1]
    zi = lb[2]
    ze = ub[2]
    ei = lb[3]
    ee = ub[3]

    '''
    xi = 122
    xe = 172
    yi = 67
    ye = 136
    zi = 16
    ze = 55
    ei = 0
    ee = 35*2

    '''
    plot_crosssection(xi, xe, yi, ye, zi, ze, ei, ee, condition*1.0)
    plt.show()

run()
