#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py

plt.rcParams['font.family'] = 'Arial'
plt.rcParams["font.size"] = 14
#plt.rcParams["xtick.top"] = True
#plt.rcParams["xtick.labeltop"] = True
fig = plt.figure(figsize=(16, 16))
fig.suptitle("crosssections of 4D INS data")


def plotter(xlim, ylim, lx, ly, data, vn, hn, cn, dev, xr, yr, hvlofs, cposinfo):
    from matplotlib.cm import get_cmap
    cmap = get_cmap("jet")
    cmap.set_under(color='k')
    ax = fig.add_subplot(vn,  hn, cn)
    print(np.max(data))

    u = np.linspace(xr[0], xr[1], data.shape[0]+1)
    v = np.linspace(yr[0], yr[1], data.shape[1]+1)
    X,Y = np.meshgrid(u, v)
    zm = np.ma.masked_where(data >= 0, data).T
    print("fuck",np.min(data))

    if int(np.max(data[xlim[0]:xlim[1], ylim[0]:ylim[1]])/dev) != 0:
        #ax.pcolor(X, Y, np.transpose(data), vmin=0, vmax=int(np.max(data[xlim[0]:xlim[1], ylim[0]:ylim[1]])/dev), cmap=cmap)
        c=ax.pcolor(X, Y, np.transpose(data), vmin=0, vmax=int(np.max(data[xlim[0]:xlim[1], ylim[0]:ylim[1]])/dev), cmap=cmap)
    else:
        #ax.pcolor(X, Y, np.transpose(data), vmin=0, vmax=np.max(data)/dev, cmap=cmap)
        c=ax.pcolor(X, Y, np.transpose(data), vmin=0, vmax=np.max(data)/dev, cmap=cmap)

    fig.colorbar(c, ax=ax)
    #c=ax.pcolor(X, Y, zm, color='k', cmap='jet')
    #ax.pcolor(X, Y, zm, color='#feff9e', cmap=cmap)

    ax.set_xlabel(lx, fontsize=18)
    ax.set_ylabel(ly, fontsize=18)
    if hvlofs:
        frac = 0.01
        yoff = float(data.shape[1])*frac
        xoff = float(data.shape[0])*frac
    else:
        frac = 0
        yoff = 0
        xoff = 0
    if cn % 2 == 0:
        xival = (xlim[0]+xoff)/float(data.shape[0])*(xr[1] - xr[0]) + xr[0]
        xeval = (xlim[1]-xoff)/float(data.shape[0])*(xr[1] - xr[0]) + xr[0]
        yival = (ylim[0]+yoff)/float(data.shape[1])*(yr[1] - yr[0]) + yr[0]
        yeval = (ylim[1]-yoff)/float(data.shape[1])*(yr[1] - yr[0]) + yr[0]
        ax.axvline(x=xival, color='#ff1ee5', lw=1.5,
                   ymin=float(ylim[0])/float(data.shape[1]),
                   ymax=float(ylim[1])/float(data.shape[1]))
        ax.axvline(x=xeval, color='#ff1ee5', lw=1.5,
                   ymin=float(ylim[0])/float(data.shape[1]),
                   ymax=float(ylim[1])/float(data.shape[1]))
        ax.axhline(y=yival, color='#ff1ee5', lw=1.5,
                   xmin=float(xlim[0])/float(data.shape[0])+frac,
                   xmax=float(xlim[1])/float(data.shape[0])-frac)
        ax.axhline(y=yeval, color='#ff1ee5', lw=1.5,
                   xmin=float(xlim[0])/float(data.shape[0])+frac,
                   xmax=float(xlim[1])/float(data.shape[0])-frac)
    if cn != 1 and cn != 2:
        ax.set_yticks([0,10,20,30,40])
    ax.xaxis.set_label_coords(0.5, -0.09)
    ax.yaxis.set_label_coords(-0.13, 0.5)
    ax.tick_params(direction="out", color="black", pad=0., labeltop=False,
                   labelbottom=True, bottom=True, top=False, left=True,
                   right=False, labelleft=True)
    ax.axis('tight')
    ax.text(xr[0]+(xr[1]-xr[0])*.2,yr[1]*1.02, cposinfo)


def plot_crosssection_more_useful(cnum, lims_int,
                                  condition, data4, ranges, hvlofs, devs,
                                  cpos):
    print("cpos",cpos)
    print("data4.shape", data4.shape)
    cposval = cpos*1.0/np.array(data4.shape)*1.0 * (ranges[:, 1]-ranges[:, 0])\
                                                               + ranges[:, 0]
    print("cposval",cposval)
    cposinfo = '@$\mathregular{(q_c, E)}$' + "=({:.0f}, {:.0f})".format(cposval[2], cposval[3])
    data2d = np.sum(np.sum(data4[:, :, cpos[2]:cpos[2]+2, cpos[3]:cpos[3]+4],axis=3),axis=2)
    condition2d = np.sum(np.sum(condition[:, :, cpos[2]:cpos[2]+2, cpos[3]:cpos[3]+4],axis=3),axis=2)
    condition2d[condition2d > 0] = 1
    print(np.max(condition2d))
    print(np.min(condition2d))
    data2d = preprocess(data2d, condition2d)
    plotter(lims_int[0,:], lims_int[1,:], '$\mathregular{q_a(rlu)}$', '$\mathregular{q_b(rlu)}$', data2d,
            3, 2, cnum, devs[0], ranges[0, :], ranges[1, :],
            hvlofs, cposinfo)
    cposinfo = '@$\mathregular{(q_a, q_c)}$'+"=({:.1f}, {:.0f})".format(cposval[0], cposval[2])
    data2d = np.sum(np.sum(data4[cpos[0]:cpos[0]+3, :, cpos[2]:cpos[2]+2, :],axis=2),axis=0)
    condition2d = np.sum(np.sum(condition[cpos[0]:cpos[0]+3, :, cpos[2]:cpos[2]+2, :],axis=2),axis=0)
    condition2d[condition2d > 0] = 1
    print(np.max(condition2d))
    print(np.min(condition2d))
    data2d = preprocess(data2d, condition2d)
    plotter(lims_int[1,:], lims_int[3,:], '$\mathregular{q_b(rlu)}$', 'E(meV)',  data2d,
            3, 2, cnum+2, devs[1], ranges[1, :], ranges[3, :],
            hvlofs, cposinfo)
    cposinfo = '@$\mathregular{(q_b, q_c)}$'+"=({:.1f}, {:.0f})".format(cposval[1],cposval[2])
    data2d = np.sum(np.sum(data4[:, cpos[1]:cpos[1]+2, cpos[2]:cpos[2]+2, :], axis=2),axis=1)
    condition2d = np.sum(np.sum(condition[:, cpos[1]:cpos[1]+2, cpos[2]:cpos[2]+2, :], axis=2),axis=1)
    condition2d[condition2d > 0] = 1
    print(np.max(condition2d))
    print(np.min(condition2d))
    data2d = preprocess(data2d, condition2d)
    plotter(lims_int[0,:],lims_int[3,:], '$\mathregular{q_a(rlu)}$', 'E(meV)',  data2d,
            3, 2, cnum+4, devs[2], ranges[0, :], ranges[3, :],
            hvlofs, cposinfo)


def preprocess(data, condition):
    data=data*condition*1.0
    maxc = np.max(condition)
    data[ condition > 0.1 ] = data[ condition > 0.1 ] / condition[ condition > 0.1] * maxc
    data[ condition <= 0.1 ] = -10000.0
    return(data)


def run2(cnum, outfile, lims, ns, ranges, hvlofs, devs,  cpos):
    lims_int = np.zeros(lims.shape, dtype=int)
    for i in range(0, lims.shape[0]):
        for j in range(0, lims.shape[1]):
            lims_int[i, j] = int(round(lims[i, j]/ns[i]))
    print("lims_int:", lims_int)
    f = h5py.File(outfile, 'r')
    data4 = np.array(f["data4"])
    condition = np.array(f["condition"])
    plot_crosssection_more_useful(cnum, lims_int,
                                  condition, data4, ranges, hvlofs, devs, cpos)


def run():
    cnum = 2
    outfile = "/home/kazu/desktop//200522/Ei42/veryfineq/14m/Output4D_00_840.hdf5"
    lims = np.array([[114, 200], [69, 122], [11, 20], [81, 207]])*1.0
    ns = np.array([1, 1, 1, 1])*1.
    ranges = np.array([[-0.675, 3.075], [-0.925, 4.375], [-0.8, 0.55], [-8.0, 36.2]])
    hvlofs = False                                              
    devs = np.array([40, 3, 3])
    cpos = np.array([175, 78, 17, 41])
    run2(cnum, outfile, lims, ns, ranges, hvlofs, devs, cpos)
    plt.subplots_adjust(wspace=0.35, hspace=0.37)


run()
#plt.savefig("Ei42_crosssections.pdf")
plt.savefig("Ei42_advsoft_crosssections.png")
plt.show()
