#!/usr/bin/env python
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
import sqwto1dspectrum_allsites_class as sac


class hpdos(sac.sqwto1dspectrum):
    def __init__(self, sites=None, projdir=None, head=None, middle=None):
        self.THztomeV = 4.13567
        self.sites = sites
        self.projdir = projdir
        self.head = head
        self.middle = middle

    def plotter(self):
        for sidx, site in enumerate(self.sites):
            infile = "/home/kazu/WORK/vasp-phonopy/la2ni10h1/lda/" + site +\
                     "/k334/partial_dos.dat_"+self.projdir
            pdos = np.genfromtxt(infile, dtype=float)
            datasize = pdos.shape
            ax = self.fig.add_subplot(self.gs[sidx, 0])
            #ax = self.fig.add_subplot(len(self.sites), 1,  sidx+1)
            ax.plot(pdos[:, 0]*self.THztomeV,
                    pdos[:, datasize[1]-1]/self.THztomeV)
            ax.text(50, 1.4, '%s' % site)
            ax.tick_params(labelbottom=False)
            ax.tick_params(direction="in")
            if sidx == len(self.sites)-1:
                ax.tick_params(labelbottom=True)
                ax.set_xlabel('Energy (meV)')
            if sidx == 2:
                ax.tick_params(labelleft=True)
                ax.set_ylabel('PDOS')
            ax.set_ylim(0, 1.7)
            ax.set_xlim(30, 210)
            ax.set_xticks([50, 100, 150, 200])
            if sidx == 1:
                ax.set_yticks([0, 1])
            else:
                ax.set_yticks([0, 1])

        plt.subplots_adjust(hspace=0.0)
        #plt.savefig("pdos_la2ni10h1_lda_k334_"+self.projdir+".pdf")


def samplerun():
    head = "/home/kazu/WORK/vasp-phonopy/la2ni10h1/lda/"
    middle = "/k334/partial_dos.dat_"
    sites = ["4h", "6m", "12o", "12n"]
    projdir = "c_axis"
    proj = hpdos(sites=sites, projdir=projdir, head=head, middle=middle)
    proj.create_fig(nr=7, suptitle="lda, La2Ni10H1 along "+projdir)
    proj.plotter()


#samplerun()
#plt.show()
