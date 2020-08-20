#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.rcParams['font.family'] = 'Arial'
#plt.rcParams["xtick.top"] = True
#plt.rcParams["xtick.labeltop"] = True
fig = plt.figure(figsize=(12, 10))
fig.suptitle("crosssections of 4D INS histograms with optimized bin widths", fontsize="x-large")


def plotter(xlim, ylim, lx, ly, data, vn, hn, cn, dev, dev2, xr, yr, hvlofs, cposinfo):
    from matplotlib.cm import get_cmap
    cmap = get_cmap("gist_gray")
    cmap.set_under(color='white')
    ax = fig.add_subplot(vn,  hn, cn)
    print(np.max(data))

    u = np.linspace(xr[0], xr[1], data.shape[0]+1)
    v = np.linspace(yr[0], yr[1], data.shape[1]+1)
    X,Y = np.meshgrid(u, v)
    zm = np.ma.masked_where(data >= 0, data).T

    #ax.pcolor(np.transpose(data), vmin=-int(np.max(data[xi:xe, yi:ye])/dev/dev2), vmax=int(np.max(data[xi:xe, yi:ye])/dev), cmap=cmap)
    ax.pcolor(X, Y, np.transpose(data), vmin=0, vmax=int(np.max(data[xlim[0]:xlim[1], ylim[0]:ylim[1]])/dev), cmap=cmap)
    ax.pcolor(X, Y, zm, color='gray', cmap=cmap)

    ax.set_xlabel(lx)
    ax.set_ylabel(ly)
    if hvlofs:
        frac = 0.01
        yoff = float(data.shape[1])*frac
        xoff = float(data.shape[0])*frac
    else:
        frac = 0
        yoff = 0
        xoff = 0
    xival = (xlim[0]+xoff)/float(data.shape[0])*(xr[1] - xr[0]) + xr[0]
    xeval = (xlim[1]-xoff)/float(data.shape[0])*(xr[1] - xr[0]) + xr[0]
    yival = (ylim[0]+yoff)/float(data.shape[1])*(yr[1] - yr[0]) + yr[0]
    yeval = (ylim[1]-yoff)/float(data.shape[1])*(yr[1] - yr[0]) + yr[0]
    ax.axvline(x=xival, color='white', lw=1.5,
               ymin=float(ylim[0])/float(data.shape[1]),
               ymax=float(ylim[1])/float(data.shape[1]))
    ax.axvline(x=xeval, color='white', lw=1.5,
               ymin=float(ylim[0])/float(data.shape[1]),
               ymax=float(ylim[1])/float(data.shape[1]))
    ax.axhline(y=yival, color='white', lw=1.5,
               xmin=float(xlim[0])/float(data.shape[0])+frac,
               xmax=float(xlim[1])/float(data.shape[0])-frac)
    ax.axhline(y=yeval, color='white', lw=1.5,
               xmin=float(xlim[0])/float(data.shape[0])+frac,
               xmax=float(xlim[1])/float(data.shape[0])-frac)
    ax.xaxis.set_label_coords(0.5, -0.09)
    ax.yaxis.set_label_coords(-0.09, 0.5)
    ax.tick_params(direction="out", color="black", pad=0., labeltop=False,
                   labelbottom=True, bottom=True, top=False, left=True,
                   right=False, labelleft=True)
    ax.axis('tight')
    ax.text((xr[1]+xr[0])*.6,yr[1]*1.02, cposinfo)


def plot_crosssection_more_useful(cnum, lims_int,
                                  condition, data4, ranges, hvlofs, devs,
                                  dev2, cpos):
    print("cpos",cpos)
    print("data4.shape", data4.shape)
    cposval = cpos*1.0/np.array(data4.shape)*1.0 * (ranges[:, 1]-ranges[:, 0])\
                                                               + ranges[:, 0]
    cposinfo = '@$\mathregular{(q_c, E)}$' + "=({:.0f}, {:.0f})".format(cposval[2], cposval[3])
    plotter(lims_int[0,:], lims_int[1,:], '$\mathregular{q_a(rlu)}$', '$\mathregular{q_b(rlu)}$', data4[:, :, cpos[2], cpos[3]],
            3, 3, cnum, devs[0], dev2, ranges[0, :], ranges[1, :],
            hvlofs, cposinfo)
    cposinfo = '@$\mathregular{(q_a, q_c)}$'+"=({:.1f}, {:.0f})".format(cposval[0], cposval[2])
    plotter(lims_int[1,:], lims_int[3,:], '$\mathregular{q_b(rlu)}$', 'E(meV)',  data4[cpos[0], :, cpos[2], :],
            3, 3, cnum+3, devs[1], dev2, ranges[1, :], ranges[3, :],
            hvlofs, cposinfo)
    cposinfo = '@$\mathregular{(q_b, q_c)}$'+"=({:.1f}, {:.0f})".format(cposval[1],cposval[2])
    plotter(lims_int[0,:],lims_int[3,:], '$\mathregular{q_a(rlu)}$', 'E(meV)',  data4[:, cpos[1], cpos[2], :],
            3, 3, cnum+6, devs[2], dev2, ranges[0, :], ranges[3, :],
            hvlofs, cposinfo)


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


def run2(cnum, outfile, lims, ns, ranges, hvlofs, devs, dev2, cpos):
    lims_int = np.zeros(lims.shape, dtype=int)
    for i in range(0, lims.shape[0]):
        for j in range(0, lims.shape[1]):
            lims_int[i, j] = int(round(lims[i, j]/ns[i]))
    print("lims_int:", lims_int)
    data4, condition = preprocess(outfile)
    plot_crosssection_more_useful(cnum, lims_int,
                                  condition, data4, ranges, hvlofs, devs, dev2, cpos)


def run():
    head = "/home/kazu/desktop/"
    tail = "hist.hdf5"

    cnum = 1
    outfile = head + "200204/fine/hourbyhour/10h/" + tail
    lims = np.array([[120, 172], [61, 145], [16, 53], [20, 70]])*1.0
    ns = np.array([4, 4, 5, 3])*1.
    ranges = np.array([[-1.65, 4.1], [-2.1, 2.8], [-0.85, 0.9], [0, 40.5]])
    hvlofs = False                                                # whether hvlines are displaced by a small value from the axis. This should be true when the line overlap on the axes.
    devs = np.array([1.2, 1, 1])
    dev2 = 70
    cpos = np.array([27, 36, 7, 0])
    run2(cnum, outfile, lims, ns, ranges, hvlofs, devs, dev2, cpos)

    cnum = 2
    outfile = head + "200522/Ei42/veryfineq/14m/" + tail
    lims = np.array([[114, 200], [69, 122], [11, 20], [81, 207]])*1.0
    ns = np.array([3, 2, 2, 4])*1.
    ranges = np.array([[-0.675, 3.075], [-0.925, 4.375], [-0.8, 0.55], [-8.0, 36.2]])
    hvlofs = False                                              
    devs = np.array([200, 1, 1])
    dev2 = 10
    cpos = np.array([58, 39, 8, 10])
    run2(cnum, outfile, lims, ns, ranges, hvlofs, devs, dev2, cpos)

    cnum = 3
    outfile = head + "200522/Ei24/fineq/26m/" + tail
    lims = np.array([[0, 228], [1, 202], [0, 8], [150, 289]])*1.0
    ns = np.array([2, 3, 2, 4])*1.
    ranges = np.array([[0.01, 2.29], [-0.67, 1.35], [-0.16, 0.16], [-2.0, 21.12]])
    hvlofs = True                                          
    devs = np.array([1, 1, 1])
    dev2 = 8
    cpos = np.array([75, 39, 2, 6])
    run2(cnum, outfile, lims, ns, ranges, hvlofs, devs, dev2, cpos)
    plt.subplots_adjust(wspace=0.25, hspace=0.27)


run()
#plt.show()
plt.savefig("no17714_ei42_ei24_crosssections_tst.pdf")
