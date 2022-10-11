#!/usr/bin/env python
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
import numpy as np
import sys
sys.path.append("/home/kazu/ktpro")
import sqwto1dspectrum_allsites_class as sac
import pickle
from matplotlib.ticker import EngFormatter

plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 2 

class hpdos(sac.sqwto1dspectrum):
    def __init__(self, infiles):
        self.infiles = infiles

    def plotter(self):
        for iidx, infile in enumerate(self.infiles):
            self.load_pkl(infile)
            ax = self.fig.add_subplot(self.gs[iidx, 0])
            ax.plot(self.dataset['ene'], self.dataset['spec'])
            ax.text(75, 1.4*10**8, '%s' % infile, fontsize=9)
            ax.tick_params(labelbottom=False)
            ax.tick_params(direction="in")
            if iidx == len(self.infiles)-1:
                ax.tick_params(labelbottom=True)
                ax.set_xlabel('Energy (meV)')
            if iidx == 1:
                ax.tick_params(labelleft=True)
                ax.set_ylabel('sqw')
            #ax.set_ylim(0, 1.7)
            ax.set_xlim(0, 400)
            ax.yaxis.set_major_formatter(EngFormatter())
            #ax.set_xticks([50, 100, 150, 200])
            #if iidx == 1:
            #    ax.set_yticks([0, 1])
            #else:
            #    ax.set_yticks([0, 1])

        plt.subplots_adjust(hspace=0.0)
        #plt.savefig("pdos_la2ni10h1_lda_k334_"+self.projdir+".pdf")

    def load_pkl(self, infile):
        #print("loading ", infile)
        with open(infile, 'rb') as f:
            self.dataset = pickle.load(f)

