#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.rcParams['font.family'] = 'Arial'
#plt.rcParams["xtick.top"] = True
#plt.rcParams["xtick.labeltop"] = True
fig = plt.figure(figsize=(9, 9))
fig.suptitle("crosssections of 4D INS histograms with optimized bin widths", fontsize="x-large")


def plotter(xi, xe, yi, ye, lx, ly, data, vn, hn, cn, dev, dev2, hvlofs):
    from matplotlib.cm import get_cmap
    cmap = get_cmap("gist_gray")
    cmap.set_under(color='black')
    ax = fig.add_subplot(vn,  hn, cn)
    print(np.max(data))
    ax.pcolor(np.transpose(data), vmin=-int(np.max(data[xi:xe, yi:ye])/dev/dev2), vmax=int(np.max(data[xi:xe, yi:ye])/dev), cmap=cmap)
    ax.set_xlabel(lx)
    ax.set_ylabel(ly)
    print("vmin", -1.0*int(np.max(data[xi:xe, yi:ye])/dev/dev2))
    print("vmax", int(np.max(data[xi:xe, yi:ye]/dev)))
    if hvlofs:
        frac = 0.01
        yoff = float(data.shape[1])*frac
        xoff = float(data.shape[0])*frac
    else:
        frac = 0
        yoff = 0
        xoff = 0
    ax.axvline(x=xi+xoff, color='white', lw=1.5,
               ymin=float(yi)/float(data.shape[1]),
               ymax=float(ye)/float(data.shape[1]))
    ax.axvline(x=xe-xoff, color='white', lw=1.5,
               ymin=float(yi)/float(data.shape[1]),
               ymax=float(ye)/float(data.shape[1]))
    ax.axhline(y=yi+yoff, color='white', lw=1.5,
               xmin=float(xi)/float(data.shape[0])+frac,
               xmax=float(xe)/float(data.shape[0])-frac)
    ax.axhline(y=ye-yoff, color='white', lw=1.5,
               xmin=float(xi)/float(data.shape[0])+frac,
               xmax=float(xe)/float(data.shape[0])-frac)
    ax.xaxis.set_label_coords(0.5, 1.1)
    ax.tick_params(direction="in", color="white", labeltop=False, labelbottom=False, bottom=False, top=False, left=False, right=False, labelleft=False)
    ax.axis('tight')


def plot_crosssection_more_useful(cnum, xi, xe, yi, ye, zi, ze, ei, ee, condition, data4, hvlofs, devs, dev2, cpos):
    plotter(xi, xe, yi, ye, 'qx', 'qy', data4[:, :, cpos[2], cpos[3]],     3, 3, cnum, devs[0], dev2, hvlofs)
    plotter(yi, ye, ei, ee, 'qy', 'E',  data4[cpos[0], :, cpos[2], :],                   3, 3, cnum+3, devs[1], dev2, hvlofs)
    plotter(xi, xe, ei, ee, 'qx', 'E',  data4[:, cpos[1], cpos[2], :],                    3, 3, cnum+6, devs[2], dev2, hvlofs)


def preprocess(outfile):
    f = h5py.File(outfile, 'r')
    data4 = np.array(f["data4"])
    print(data4.shape)
    condition = np.array(f["condition"])
    #data4 = np.where(condition == False, -1.0, data4)
    maxc = np.max(condition)
    for ix in range(0, condition.shape[0]):
        for iy in range(0, condition.shape[1]):
            for iz in range(0, condition.shape[2]):
                for iw in range(0, condition.shape[3]):
                    if condition[ix, iy, iz, iw] > 0.1:
                        data4[ix, iy, iz, iw] = data4[ix, iy, iz, iw] / condition[ix, iy, iz, iw] * maxc
                    else:
                        data4[ix, iy, iz, iw] = -10000.0
    return(data4, condition)
    #plot_crosssection_more_useful(xi, xe, yi, ye, zi, ze, ei, ee, condition, data4)
    #plt.savefig(figfile)


def run2(cnum, outfile, lims, ns, hvlofs, devs, dev2, cpos):
    lims_int = np.zeros(lims.shape, dtype=int)
    for i in range(0, lims.shape[0]):
        for j in range(0, lims.shape[1]):
            lims_int[i, j] = int(round(lims[i, j]/ns[i]))
    print(lims_int)
    data4, condition = preprocess(outfile)
    plot_crosssection_more_useful(cnum, lims_int[0, 0], lims_int[0, 1],
                                  lims_int[1, 0], lims_int[1, 1],
                                  lims_int[2, 0], lims_int[2, 1],
                                  lims_int[3, 0], lims_int[3, 1],
                                  condition, data4, hvlofs, devs, dev2, cpos)


def run():
    head = "/home/kazu/desktop/"
    tail = "hist.hdf5"

    cnum = 1
    outfile = head + "200204/fine/hourbyhour/10h/" + tail
    lims = np.array([[120, 172], [61, 145], [16, 53], [20, 70]])*1.0
    ns = np.array([4, 4, 5, 3])*1.
    hvlofs = False                                                # whether hvlines are displaced by a small value from the axis. This should be true when the line overlap on the axes.
    devs = np.array([1.2, 1, 1])
    dev2 = 70
    cpos = np.array([24, 36, 7, 0])
    run2(cnum, outfile, lims, ns, hvlofs, devs, dev2, cpos)

    cnum = 2
    outfile = head + "200522/Ei42/veryfineq/14m/" + tail
    lims = np.array([[114, 200], [69, 122], [11, 20], [81, 207]])*1.0
    ns = np.array([3, 2, 2, 4])*1.
    hvlofs = False                                              
    devs = np.array([200, 1, 1])
    dev2 = 10
    cpos = np.array([58, 39, 8, 10])
    run2(cnum, outfile, lims, ns, hvlofs, devs, dev2, cpos)

    cnum = 3
    outfile = head + "200522/Ei24/fineq/26m/" + tail
    lims = np.array([[0, 228], [1, 202], [0, 8], [150, 289]])*1.0
    ns = np.array([2, 3, 2, 4])*1.
    hvlofs = True                                          
    devs = np.array([1, 1, 1])
    dev2 = 8
    cpos = np.array([75, 39, 2, 6])
    run2(cnum, outfile, lims, ns, hvlofs, devs, dev2, cpos)


run()
plt.show()
