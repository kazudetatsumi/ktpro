#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams['font.family'] = 'Arial'
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)


class CROSS:

    def __init__(self, outfile, orthotope_lims, wholeranges, devs,
                 cpos, hvlofs=False, binwidths=np.ones(4)):
        self.outfile = outfile
        self.orthotope_lims = orthotope_lims
        self.wholeranges = wholeranges
        self.devs = devs
        self.cpos = cpos
        self.hvlofs = hvlofs
        self.binwidths = binwidths

    def read_and_preprocess_data(self):
        f = h5py.File(self.outfile, 'r')
        data4 = np.array(f["data4"])
        condition = np.array(f["condition"])
        maxc = np.max(condition)
        data4[condition > 0.1] = data4[condition > 0.1] /\
            condition[condition > 0.1] * maxc
        data4[condition <= 0.1] = -10000.0
        return(data4, condition)

    def lims_int(self):
        lims_int = np.zeros(self.orthotope_lims.shape, dtype=int)
        for i in range(0, self.orthotope_lims.shape[0]):
            for j in range(0, self.orthotope_lims.shape[1]):
                lims_int[i, j] = int(round(self.orthotope_lims[i, j]
                                     / self.binwidths[i]))
        return(lims_int)

    def plot_crosssection(self, cnum, xyyx=False):
        data4, condition = self.read_and_preprocess_data()
        lims_int = self.lims_int()
        if xyyx:
            data4 = np.transpose(data4, (1, 0, 2, 3))
            self.cpos[[0, 1]] = self.cpos[[1, 0]]
            lims_int[[0, 1], :] = lims_int[[1, 0], :]
            self.wholeranges[[0, 1], :] = self.wholeranges[[1, 0], :]
        cposval = self.cpos*1.0/np.array(data4.shape)*1.0 *\
            (self.wholeranges[:, 1]-self.wholeranges[:, 0]) +\
            self.wholeranges[:, 0]
        cposinfo = '$(q_c, E)$' + "=({:.0f}, {:.0f})".format(
                   cposval[2], cposval[3])
        self.plotter(lims_int[0, :], lims_int[1, :], '$q_a(rlu)$', '$q_b(rlu)$',
                     data4[:, :, self.cpos[2], self.cpos[3]],
                     4, 3, cnum, self.devs[0], self.wholeranges[0, :],
                     self.wholeranges[1, :], cposinfo)
        cposinfo = '$(q_a, q_c)$' + "=({:.1f}, {:.0f})".format(
                   cposval[0], cposval[2])
        self.plotter(lims_int[1, :], lims_int[3, :], '$q_b(rlu)$', 'E(meV)',
                     data4[self.cpos[0], :, self.cpos[2], :],
                     4, 3, cnum+3, self.devs[1], self.wholeranges[1, :],
                     self.wholeranges[3, :], cposinfo)
        cposinfo = '$(q_b, q_c)$' + "=({:.1f}, {:.0f})".format(
                   cposval[1], cposval[2])
        self.plotter(lims_int[0, :], lims_int[3, :], '$q_a(rlu)$', 'E(meV)',
                     data4[:, self.cpos[1], self.cpos[2], :],
                     4, 3, cnum+6, self.devs[2], self.wholeranges[0, :],
                     self.wholeranges[3, :], cposinfo)
        cposinfo = '$(q_a, q_b)$' + "=({:.1f}, {:.0f})".format(
                   cposval[0], cposval[1])
        self.plotter(lims_int[2, :], lims_int[3, :], '$q_c(rlu)$', 'E(meV)',
                     data4[self.cpos[0], self.cpos[1], :, :],
                     4, 3, cnum+9, self.devs[2], self.wholeranges[2, :],
                     self.wholeranges[3, :], cposinfo)

    def plotter(self, xlim, ylim, lx, ly, data, vn, hn, cn, dev, xr, yr,
                cposinfo):

        cdict = {'red':   ((0.0, 0.3, 0.3),
                           (1.0, 1.0, 1.0)),

                 'green': ((0.0, 0.3, 0.3),
                           (1.0, 1.0, 1.0)),

                 'blue':  ((0.0, 0.3, 0.3),
                           (1.0, 1.0, 1.0))
                 }

        custom_cmap = LinearSegmentedColormap('truncated_gray', cdict)
        custom_cmap.set_under(color='black')
        ax = self.fig.add_subplot(vn,  hn, cn)

        u = np.linspace(xr[0], xr[1], data.shape[0]+1)
        v = np.linspace(yr[0], yr[1], data.shape[1]+1)
        X, Y = np.meshgrid(u, v)

        c = ax.pcolor(X, Y, np.transpose(data), vmin=0,
                      vmax=int(np.max(data[xlim[0]:xlim[1], ylim[0]:ylim[1]])
                               / dev),
                      cmap=custom_cmap)

        #ax.set_xlabel(lx, fontsize=12)
        #ax.set_ylabel(ly, fontsize=12)
        if self.hvlofs:
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
                       right=False, labelleft=True, width=1.5)
        #ax.axis('tight')
        ax.text((xr[1]+xr[0])*.6, yr[1]*1.04, cposinfo)
        self.fig.colorbar(c, ax=ax)

    def create_fig(self, title=None):
        self.fig = plt.figure(figsize=(8.8, 8.8))
        if title is None:
            self.fig.suptitle("crosssections of 4D INS intensity initial " +
                              "aray for bin-widths optimization",
                              fontsize="x-large")
        else:
            self.fig.suptitle(title, fontsize="x-large")


def samplerun():

    outfile = "/home/kazu/desktop/200204/fine/hourbyhour/10h/out_hw_all.hdf5"
    orthotope_lims = np.array([[120, 172], [61, 145], [16, 53], [20, 70]])*1.0
    orgbinwidths = np.array([4, 4, 5, 3])
    binwidths = np.array([1, 1, 1, 1])*1.
    wholeranges = np.array([[-1.65, 4.1], [-2.1, 2.8], [-0.85, 0.9], [0, 40.5]]
                           )
    devs = np.array([1, 1, 1])
    cpos = np.array([27, 36, 7, 0])*orgbinwidths
    cpos = np.array([107, 145, 34, 0])

    pro = CROSS(outfile, orthotope_lims, wholeranges, devs, cpos, hvlofs=False,
                binwidths=binwidths)
    pro.create_fig()
    pro.plot_crosssection(1, xyyx=True)

    pro.outfile = "/home/kazu/desktop/200522/Ei42/veryfineq/14m/" + \
                  "Output4D_00_840.hdf5"
    pro.orthotope_lims = np.array([[114, 200], [69, 122], [11, 20],
                                  [81, 207]])*1.0
    orgbinwidths = np.array([3, 2, 2, 4])
    pro.binwidths = np.array([1, 1, 1, 1])*1.
    pro.wholeranges = np.array([[-0.675, 3.0875], [-0.925, 4.400],
                                [-0.8, 0.60], [-8.0, 36.4]])
    pro.devs = np.array([1, 1, 1])
    pro.cpos = np.array([58, 39, 8, 10])*orgbinwidths

    pro.plot_crosssection(2)

    pro.outfile = "/home/kazu/desktop/200522/Ei24/fineq/26m/" + \
                  "Output4D_00_1560.hdf5"
    pro.orthotope_lims = np.array([[0, 228], [0, 202], [0, 9], [150, 289]])*1.0
    orgbinwidths = np.array([2, 3, 2, 4])
    pro.binwidths = np.array([1, 1, 1, 1])*1.
    pro.wholeranges = np.array([[0.01, 2.31], [-0.67, 1.35], [-0.16, 0.18],
                                [-2.0, 21.20]])
    pro.devs = np.array([1, 1, 1])
    pro.cpos = np.array([75, 39, 2, 6])*orgbinwidths

    pro.plot_crosssection(3)


#samplerun()
#plt.savefig("test.png")
#plt.show()
