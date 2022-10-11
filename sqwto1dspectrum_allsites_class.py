#!/usr/bin/env python
# This script plot 1D map of S(q, w) from the OCLIMAX output csv files
# containing 2D intensity distribution of S(q, w).
# The paths of the csv files are hard coded and u should alter it for proper
# operation.
# 2020/06/23 Kazuyoshi TATSUMI
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
import interpolate_and_subt_class as isc
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)


plt.rcParams['font.size'] = 18
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 2


class sqwto1dspectrum(isc.interpolate_and_subt):
    def __init__(self, qmin=1, qmax=8, sites=None, head=None, tail=None,
                 middles=None):
        self.qmin = qmin
        self.qmax = qmax
        self.cmtomeV = 0.123984  # [meV cm]
        self.sites = sites
        self.head = head
        self.tail = tail
        self.middles = middles

    def get_data(self, infile):
        # The calculated output file of olimax contains header lines of 4, then
        # q, w, intensity are printed line by line.
        # Those lines with the same q values are separated by "\n" line.
        # In the contract of the lines with the same q, w are increased by a
        # constant value line by line.

        # First seek for the number of w.
        f = open(infile)
        for i, line in enumerate(f):
            if line == "\n":
                numomega = i - 4
                break
        f.close
        # Second seek for the omega values.
        omega = np.zeros(numomega)
        intensity = np.zeros(numomega)
        f = open(infile)
        for i, line in enumerate(f):
            if len(line.split(',')) == 3 and i >= 4:
                omega[i - 4] = float(line.split(',')[1])
            if line == "\n":
                break
        f.close
        # Third seek for the intensities wihtin qmin <= q <= qmax for each w
        # and sum them.
        # Finally, we obtain omega and intensity of a 1D spectrum.
        f = open(infile)
        j = 0
        for i, line in enumerate(f):
            if len(line.split(',')) == 3 and i >= 4:
                if float(line.split(',')[0]) >= self.qmin and\
                   float(line.split(',')[0]) <= self.qmax and\
                   float(line.split(',')[2]) > 0:
                    intensity[j] += float(line.split(',')[2])
                    j = j + 1
            elif line == "\n":
                j = 0
        # plt.plot(omega, intensity)
        data = np.zeros((numomega, 2))
        data[:, 0] = omega
        data[:, 1] = intensity
        return(data)

    def create_fig(self, nr=7, figsize=(6, 15), suptitle=None):
        if suptitle is None:
            suptitle = "inc q between "+str(self.qmin)+" and "+str(self.qmax) \
                       + " $\AA^{-1}$"
        self.fig = plt.figure(figsize=figsize)
        #self.fig.suptitle(suptitle)
        self.gs = self.fig.add_gridspec(nr, 1, hspace=0, wspace=0)

    def plotter(self, bottomno=None):
        for isite,  site in enumerate(self.sites):
            #ax = self.fig.add_subplot(len(self.sites)+1, 1, isite + 1)
            ax = self.fig.add_subplot(self.gs[isite, 0])
            for imid, middle in enumerate(self.middles[isite]):
                ax.plot(self.alldata[isite, imid, :, 0],
                        self.alldata[isite, imid, :, 1], label=middle, c='k',
                        lw=2)
            ax.legend(fontsize="xx-small")
            ax.set_ylim(0, np.max(self.alldata[isite, :, :, 1])*1.05)
            ax.set_xlim(0, 250)
            ax.text(10, np.max(self.alldata[isite, :, :, 1])*0.6, site)
            ax.set_yticks([])
            if isite == bottomno:
                ax.set_xlabel('Neutron Energy Loss (meV)')
                ax.tick_params(top=True, right=True, direction='in',
                               which='both', labelbottom=True, width=2)
            else:
                ax.tick_params(top=True, right=True, direction='in',
                               which='both', labelbottom=False, width=2)
            if isite == 3:
                ax.set_ylabel('Intensity')

    def get_alldata(self):
        for isite,  site in enumerate(self.sites):
            for imid, middle in enumerate(self.middles[isite]):
                data = self.get_data(self.head + site + middle + self.tail)
                if imid == 0 and isite == 0:
                    self.alldata = np.zeros((len(self.sites),
                                             len(self.middles[isite]),
                                             data.shape[0], data.shape[1]))
                self.alldata[isite, imid, :, :] = data

    def sumplot(self):
        #ax = self.fig.add_subplot(len(self.sites)+1, 1, len(self.sites) + 1)
        ax = self.fig.add_subplot(self.gs[6, 0])
        xti4 = 0.20
        xv4 = 0.05
        for imid, middle in enumerate(self.middles[0][:]):
            x = self.alldata[1, imid, :, 0]
            y = self.alldata[1, imid, :, 1]*(1-xti4-xv4) +\
                self.alldata[4, 0, :, 1] * xv4 +\
                self.alldata[0, imid, :, 1]*xti4
            ax.plot(x, y, label=middle, c='k', lw=2)
        ax.set_yticks([])
        ax.set_ylim(0, np.max(y)*1.05)
        ax.set_xlim(0, 250)
        ax.set_xlabel('Neutron Energy Loss (meV)')
        ax.text(10, np.max(y)*0.9, "ti3sbh3 * 0.75 + 24k4_6dv1 * 0.05 +\
                ti4sbh2*0.20", fontsize="xx-small")
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=True, width=2)

    def pngplot(self, pngfile, nr=2, extent=(-1.63, 1.63, -1, 1)):
        ax = self.fig.add_subplot(self.gs[nr:, 0])
        ax.imshow(mpimg.imread(pngfile),  extent=extent)
        ax.set_axis_off()

    def polylani5(self, nr=7):
        head = "/home/kazu/WORK/vasp-phonopy/la2ni10h1/"
        coarse_f = head + "INS_polyLaNi5_bg.csv"
        fine_f = head + "INS_polyLaNi5.csv"
        super(sqwto1dspectrum, self).__init__(coarse_f, fine_f)
        super(sqwto1dspectrum, self).interpolate_subt()
        ax = self.fig.add_subplot(self.gs[nr, 0])
        ax.plot(self.fd[:, 0], self.fd[:, 1])
        ax.plot(self.cd[:, 0], self.cd[:, 1])
        ax.set_xlim(0, 250)
        ax.tick_params(top=True, right=True, direction='in',
                       which='both', labelbottom=False, width=2)
        ax = self.fig.add_subplot(self.gs[nr+1, 0])
        ax.plot(self.sbd[:, 0], self.sbd[:, 1])
        ax.set_xlim(0, 250)
        ax.set_xlabel('Neutron Energy Loss (meV)')
        ax.tick_params(top=True, right=True, direction='in', which='both',
                       labelbottom=True, width=2)
        



def samplerun():
    # q range for integration
    qmin = 1  # [Angs-1]
    qmax = 8  # [Angs-1]
    sites = ["ti4sb2h", "ti3sbh3", "ti3sbh3_24k2_6dv1", "ti3sbh3_24k3_6dv1",
             "ti3sbh3_24k4_6dv1", "ti3sbh3_16i"]
    middles = [["/gga/inc-coh/", "/ggasol/inc-coh/",
                "/lda/I4/optagain/inc-coh/"],
               ["/gga/inc-coh/", "/ggasol/inc-coh/", "/lda/inc-coh/"],
               ["/gga/inc-coh/", "/ggasol/inc-coh/", "/lda/inc-coh/"],
               ["/gga/inc-coh/", "/ggasol/inc-coh/", "/lda/inc-coh/"],
               ["/gga/inc-coh/", "/ggasol/inc-coh/", "/lda/inc-coh/"],
               ["/gga/inc-coh/", "/ggasol/inc-coh/", "/lda/inc-coh/"]]
    middles = [["/lda/I4/optagain/inc-coh/"],
               ["/lda/inc-coh/"],
               ["/lda/inc-coh/"],
               ["/lda/inc-coh/"],
               ["/lda/inc-coh/"],
               ["/lda/inc-coh/"]]
    head = "/home/kazu/WORK/vasp-phonopy/"
    tail = "/out_2Dmesh_coh_4K.csv"
    prj = sqwto1dspectrum(qmin=qmin, qmax=qmax, sites=sites, head=head,
                          middles=middles, tail=tail)
    prj.create_fig()
    prj.get_alldata()
    prj.plotter()
    prj.sumplot()
    plt.show()


#samplerun()
