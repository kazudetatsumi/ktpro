#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Times New Roman'
#plt.rcParams["xtick.top"] = True
#plt.rcParams["xtick.labeltop"] = True


def getdata(f):
    data = np.genfromtxt(f, dtype=float, delimiter=",")
    return data


def make_mappable(maxvalue):
    from matplotlib.colors import Normalize
    norm = Normalize(vmin=0, vmax=maxvalue)
    from matplotlib.cm import ScalarMappable, get_cmap
    cmap = get_cmap("jet")
    mappable = ScalarMappable(norm=norm, cmap=cmap)
    mappable._A = []
    return mappable


def plotter(fd, cd, fds):
    fig=plt.figure(figsize=(6, 12))
    ax=fig.add_subplot(2, 1, 1)
    ax.set_title("LaNi5 INS experimental spectrum")
    ax.plot(fd[:, 0], fd[:, 1])
    ax.plot(cd[:, 0], cd[:, 1])
    ax.tick_params(direction = "in")
    ax.set_ylabel('S(q,$\omega$)')
    ax.text(210, 4, 'S(q, $\omega$)',color="dodgerblue")
    ax.text(50, 20, 'bg',color="orange")
    ax=fig.add_subplot(2, 1, 2)
    ax.plot(fds[:, 0], fds[:, 1])
    #ax.set_ylim(0, 0.06)
    plt.subplots_adjust(wspace=0.4, hspace=0.0)
    ax.set_xlabel('Energy (meV)')
    ax.set_ylabel('Intensity')
    ax.text(190, 1, 'S(q, $\omega$) - bg', color="dodgerblue")
    ax.tick_params(direction = "in")




def run():
    coarse_f="/home/kazu/INS_polyLaNi5_bg.csv"
    fine_f="/home/kazu/INS_polyLaNi5.csv"
    coarse_data=getdata(coarse_f)
    fine_data=getdata(fine_f)
    cint = np.zeros(fine_data.shape[0])
    interpolated = np.zeros_like(fine_data)
    fine_data_subtracted = np.zeros_like(fine_data)
    for findx in range(0, fine_data.shape[0]):
        interpolated[findx, 0] = fine_data[findx, 0]
        fine_data_subtracted[findx, 0] = fine_data[findx, 0]
        for cindx in range(0, coarse_data.shape[0]-1):
            if fine_data[findx, 0] >= coarse_data[cindx, 0] and fine_data[findx, 0] < coarse_data[cindx + 1, 0]:
                interpolated[findx, 1] = ( coarse_data[cindx+1, 1] - coarse_data[cindx, 1 ]) /( coarse_data[cindx+1, 0] - coarse_data[cindx, 0] ) * ( fine_data[findx, 0] -  coarse_data[cindx, 0] ) + coarse_data[cindx, 1]
                fine_data_subtracted[findx, 1] = fine_data[findx, 1] - interpolated[findx, 1]
    plotter(fine_data, coarse_data, fine_data_subtracted)



run()
plt.show()
#plt.savefig("crosssection.eps")
