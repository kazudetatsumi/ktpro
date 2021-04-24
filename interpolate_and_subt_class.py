#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

# plt.rcParams['font.family'] = 'Times New Roman'
# plt.rcParams["xtick.top"] = True
# plt.rcParams["xtick.labeltop"] = True


class interpolate_and_subt:
    def __init__(self, coarse_f, fine_f):
        self.coarse_f = coarse_f
        self.fine_f = fine_f

    def getdata(self, f):
        data = np.genfromtxt(f, dtype=float, delimiter=",")
        return data

    def plotter(self):
        fig = plt.figure(figsize=(6, 8))
        ax = fig.add_subplot(2, 1, 1)
        ax.set_title("LaNi5 INS experimental spectrum")
        ax.plot(self.fd[:, 0], self.fd[:, 1])
        ax.plot(self.cd[:, 0], self.cd[:, 1])
        ax.tick_params(direction="in")
        ax.set_ylabel('S(q,$\omega$)')
        ax.text(210, 4, 'S(q, $\omega$)', color="blue")
        ax.text(50, 20, 'bg', color="orange")
        ax = fig.add_subplot(2, 1, 2)
        ax.plot(self.sbd[:, 0], self.sbd[:, 1])
        # ax.set_ylim(0, 0.06)
        plt.subplots_adjust(wspace=0.4, hspace=0.0)
        ax.set_xlabel('Energy (meV)')
        ax.set_ylabel('Intensity')
        ax.text(190, 2, 'S(q, $\omega$) - bg', color="blue")
        ax.tick_params(direction="in")

    def interpolate_subt(self):
        fd = self.getdata(self.fine_f)
        cd = self.getdata(self.coarse_f)
        ind = np.zeros_like(fd[:, 0])
        sbd = np.zeros_like(fd)
        sbd[:, 0] = fd[:, 0]
        for fidx in range(0, fd.shape[0]):
            for cidx in range(0, cd.shape[0]-1):
                if fd[fidx, 0] >= cd[cidx, 0] and fd[fidx, 0] < cd[cidx+1, 0]:
                    ind[fidx] = (cd[cidx+1, 1] - cd[cidx, 1]) /\
                                (cd[cidx+1, 0] - cd[cidx, 0]) *\
                                (fd[fidx, 0] - cd[cidx, 0]) + cd[cidx, 1]
                    sbd[fidx, 1] = fd[fidx, 1] - ind[fidx]
        self.fd = fd
        self.cd = cd
        self.sbd = sbd


def samplerun():
    head = "/home/kazu/WORK/vasp-phonopy/la2ni10h1/"
    coarse_f = head + "INS_polyLaNi5_bg.csv"
    fine_f = head + "INS_polyLaNi5.csv"
    prj = interpolate_and_subt(coarse_f, fine_f)
    prj.interpolate()
    prj.plotter()
    plt.show()


#samplerun()
