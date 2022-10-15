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

    def plotter(self, elow=0, ehigh=400, tx=75, ty=1.4*10**8,
                noxlabel=False, numoffolds=2, short=False):
        for iidx, infile in enumerate(self.infiles):
            self.load_pkl(infile)
            ax = self.fig.add_subplot(self.gs[iidx, 0])
            ax.plot(self.dataset['ene'], self.dataset['spec'])
            tmpylim = ax.get_ylim()
            ax.set_ylim(0, tmpylim[1])
            if short:
                words = [infile.split('/')[-1]]
                if iidx == 0:
                    ax.text(elow, tmpylim[1]*1.05, '%s' % infile.split(words[0])[0], fontsize=9)
                ty=tmpylim[1]*0.8
            else:
                words = self.getwords(infile, numoffolds)
            for iw, word in enumerate(words):
                ax.text(tx, ty-iw*ty/12, '%s' % word, fontsize=9)
            ax.tick_params(labelbottom=False)
            ax.tick_params(direction="in")
            if iidx == len(self.infiles)-1 and not noxlabel:
                ax.tick_params(labelbottom=True)
                ax.set_xlabel('Energy (meV)')
            if iidx == 1:
                ax.tick_params(labelleft=True)
                ax.set_ylabel('sqw')
            #ax.set_ylim(0, 1.7)
            ax.set_xlim(elow, ehigh)
            ax.yaxis.set_major_formatter(EngFormatter())
            #ax.set_xticks([50, 100, 150, 200])
            #if iidx == 1:
            #    ax.set_yticks([0, 1])
            #else:
            #    ax.set_yticks([0, 1])

        plt.subplots_adjust(hspace=0.0)
        #plt.savefig("pdos_la2ni10h1_lda_k334_"+self.projdir+".pdf")

    def getwords(self, infile, numoffolds):
        words = infile.split('/')[1:]
        lens = []
        lenss = [0]
        _il = 0
        _words = []
        for i in range(numoffolds):
            if i == numoffolds + 1:
                lens.append(len(words) % numoffolds)
            else:
                lens.append(len(words) // numoffolds)
        for il, _len in enumerate(lens):
            _il += _len
            lenss.append(_il)
        for il in range(len(lenss)):
            _word = ""
            if il == len(lenss)-1:
                for ill in range(lenss[il], len(words)):
                    _word += "/"
                    _word += words[ill]
            else:
                for ill in range(lenss[il], lenss[il+1]):
                    _word += "/"
                    _word += words[ill]
            _words.append(_word)
        return(_words)

    def oclimaxdataplotter(self, oclimaxfile, elow=0, ehigh=400, tx=75,
                           ty=1.4*10**8, numoffolds=2):
        self.qmin = 1
        self.qmax = 8
        data = self.get_data(oclimaxfile)
        words = self.getwords(oclimaxfile, numoffolds)
        ax = self.fig.add_subplot(self.gs[1, 0])
        ax.plot(data[:, 0], data[:, 1])
        ax.set_xlim(elow, ehigh)
        ax.tick_params(direction="in")
        ax.tick_params(labelbottom=True)
        ax.set_xlabel('Energy (meV)')
        for iw, word in enumerate(words):
            ax.text(tx, ty-iw*ty/12, '%s' % word, fontsize=9)

    def load_pkl(self, infile):
        #print("loading ", infile)
        with open(infile, 'rb') as f:
            self.dataset = pickle.load(f)

