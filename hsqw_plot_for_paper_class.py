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
plt.rcParams['axes.linewidth'] = 1.6

class hpdos(sac.sqwto1dspectrum):
    def __init__(self, infiles):
        self.infiles = infiles
        self.THztomeV = 4.13567

    def plotter(self, elow=0, ehigh=400, tx=75, ty=1.4*10**8,
                noxlabel=False, numoffolds=2, short=False, nr=0,
                c='#1f77b4'):
        for iidx, infile in enumerate(self.infiles):
            #print(infile, iidx)
            self.load_pkl(infile)
            ax = self.fig.add_subplot(self.gs[iidx+nr, 0])
            ax.plot(self.dataset['ene'], self.dataset['spec'], c='k')
            tmpylim = ax.get_ylim()
            ax.set_ylim(0, tmpylim[1])
            if short:
                words = [infile.split('/')[-1]]
                if iidx == 0:
                    ax.text(elow, tmpylim[1]*1.05, '%s' %
                            infile.split(words[0])[0], fontsize=9)
                ty = tmpylim[1]*0.8
            else:
                ty = tmpylim[1]*0.9
                words = self.getwords(infile.split('vasp')[-1], numoffolds)
            for iw, word in enumerate(words):
                ax.text(tx, ty-iw*ty/8, '%s' % word, fontsize=8)
            ax.tick_params(labelbottom=False)
            ax.tick_params(direction="in", length=6, top=True, width=1.6)
            if iidx == len(self.infiles)-1 and not noxlabel:
                ax.tick_params(axis='x', direction="inout", length=12, top=True, width=1.6)
                ax.tick_params(axis='y', direction="in", length=6, top=True, width=1.6)
                #ax.tick_params(labelbottom=False)
                #ax.set_xlabel('Energy (meV)')
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

    def pdosplotter(self, elow=0, ehigh=400, tx=75, ty=1.4*10**8,
                    noxlabel=False, numoffolds=2, short=False, nr=3):
        for iidx, infile in enumerate(self.infiles):
            pdos = np.genfromtxt(infile, dtype=float)
            datasize = pdos.shape
            ax = self.fig.add_subplot(self.gs[nr+iidx, 0])
            ax.plot(pdos[:, 0]*self.THztomeV,
                    pdos[:, datasize[1]-1]/self.THztomeV)
            ax.text(tx, ty, '%s' % infile)
            ax.tick_params(labelbottom=False)
            ax.tick_params(direction="in")
            tmpylim = ax.get_ylim()
            ax.set_ylim(0, tmpylim[1])
            if iidx == len(self.infiles)-1:
                ax.tick_params(labelbottom=True)
                ax.set_xlabel('Energy (meV)')
                ax.tick_params(labelleft=True)
                ax.set_ylabel('PDOS')
            #ax.set_ylim(0, 1.7)
            ax.set_xlim(elow, ehigh)

        plt.subplots_adjust(hspace=0.0)

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
                           ty=1.4*10**8, numoffolds=2, nr=0):
        self.qmin = 1
        self.qmax = 8
        data = self.get_data(oclimaxfile)
        words = self.getwords(oclimaxfile.split('vasp-phonopy')[-1],
                              numoffolds)
        ax = self.fig.add_subplot(self.gs[1+nr, 0])
        ax.plot(data[:, 0], data[:, 1], c='k')
        ax.set_xlim(elow, ehigh)
        ax.set_ylim(0, 150)
        #ax.tick_params(direction="in")
        ax.tick_params(axis='x',direction="inout", length=12, top=True, width=1.6)
        ax.tick_params(axis='y',direction="in", length=6, top=True, width=1.6)
        ax.tick_params(labelbottom=True)
        ax.set_xlabel('Energy (meV)')
        for iw, word in enumerate(words):
            ax.text(tx, ty-iw*ty/12, '%s' % word, fontsize=8)

    def load_pkl(self, infile):
        #print("loading ", infile)
        with open(infile, 'rb') as f:
            self.dataset = pickle.load(f)

    def get_alldata(self):
        for ifx, infile in enumerate(self.infiles):
            self.load_pkl(infile)
            if ifx == 0:
                self.all_data = np.zeros((len(self.infiles)+1,
                                          self.dataset['ene'].shape[0]))
                self.all_data[ifx, :] = self.dataset['ene']
            self.all_data[ifx+1, :] = self.dataset['spec']

    def sumplot(self, fracs, elow=0, ehigh=400, tx=75, ty=1.4*10**8,
                numoffolds=2, nr=0, text=None):
        ax = self.fig.add_subplot(self.gs[1+nr, 0])
        ax.plot(self.all_data[0],
                np.matmul(self.all_data[1:len(self.infiles)+1].T,
                          np.array(fracs)), c='k')
        ax.set_xlim(elow, ehigh)
        tmpylim = ax.get_ylim()
        ax.set_ylim(0, tmpylim[1])
        ax.tick_params(axis='x',direction="inout", length=12, top=True, width=1.6)
        ax.tick_params(axis='y',direction="in", length=6, top=True, width=1.6)
        #ax.tick_params(direction="in")
        ax.tick_params(labelbottom=True)
        ax.yaxis.set_major_formatter(EngFormatter())
        if text:
            ty = tmpylim[1]*0.9
            ax.text(tx, ty-ty/12, '%s' % text, fontsize=8)
        ax.set_xlabel('Energy (meV)')
