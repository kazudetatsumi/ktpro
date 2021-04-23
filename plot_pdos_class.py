#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import yaml

plt.rcParams['font.size'] = 18
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 2


class pdosplot:
    def __init__(self):
        PlanckConstant = 4.13566733e-15       # [eV s]
        self.THzTomev = PlanckConstant * 1e15      # [meV]

    def get_pdos(self, infiledir):
        self.infiledir = infiledir
        ene = []
        pdos = []
        for line in open(self.infiledir+"projected_dos.dat"):
            if "#" not in line:
                ene.append(float(line.split()[0]))
                pdos.append([float(_) for _ in line.split()[1:]])
        self.ene = np.array(ene)*self.THzTomev
        self.pdos = np.array(pdos) / self.THzTomev

    def get_atoms(self):
        with open(self.infiledir+"phonopy.yaml", 'r') as f:
            documents = yaml.safe_load(f)
            natoms = len(documents["primitive_cell"]["points"])
            element = []
            for aidx in range(0, natoms):
                element.append(documents["primitive_cell"]["points"][aidx]
                               ["symbol"])
        self.element = np.array(element)
   
    def create_fig(self, title="pdos of Ti4Sb2H and Ti3SbH3"):
        self.fig = plt.figure(figsize=(8, 8))
        self.fig.suptitle(title)

    def plotter(self, tnr, nr, labelbottom=False):
        colordict = {'H': 'green', 'Ti': 'gray', 'Sb': 'k'}
        ax = self.fig.add_subplot(tnr, 1, nr)
        for ele in np.unique(self.element):
            ax.plot(self.ene, np.sum(self.pdos[:, self.element == ele], axis=1
                                     ), label=ele, c=colordict[ele])
            ax.set_ylim(0, np.max(np.sum(self.pdos[:, self.element == ele],
                        axis=1)))
            if ele == 'H':
                ax.fill_between(self.ene, np.sum(self.pdos[:, self.element
                                == ele], axis=1), color='lightgreen')
        ax.set_ylabel('pdos [1/meV]')
        ax.set_xlim(-65, 200)
        text = ""
        for _ in self.infiledir.split('/')[5:]:
            text += "/" + _
        ax.text(50, np.max(np.sum(self.pdos[:, self.element == ele], axis=1))
                * 0.75, text)
        plt.tick_params(top=True, right=True, direction='in', which='both',
                        labelbottom=labelbottom, width=2)
        if labelbottom:
            ax.set_xlabel('phonon energy (meV)')
        plt.legend()
        plt.subplots_adjust(wspace=0.16, hspace=0.0)


def samplerun():
    infiledir1 = "/home/kazu/WORK/vasp-phonopy/ti4sb2h/lda/"
    infiledir2 = "/home/kazu/WORK/vasp-phonopy/ti4sb2h/lda/I4/optagain/"
    infiledir3 = "/home/kazu/WORK/vasp-phonopy/ti3sbh3/lda/"
    prj = pdosplot()
    prj.get_pdos(infiledir1)
    prj.get_atoms()
    prj.create_fig()
    prj.plotter(3, 1)
    prj.get_pdos(infiledir2)
    prj.get_atoms()
    prj.plotter(3, 2)
    prj.get_pdos(infiledir3)
    prj.get_atoms()
    prj.plotter(3, 3, labelbottom=True)
    plt.show()


#samplerun()

