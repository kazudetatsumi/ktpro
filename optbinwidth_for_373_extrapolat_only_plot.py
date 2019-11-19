#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'


def run():
    num_txtfiles = 16
    infile = "result.txt.full"
    fig=plt.figure(figsize=(8, 20))
    f = open(infile, 'r')
    extraFound = 0
    line = f.readline()
    xlist_for_each_n = []
    ylist_for_each_n = []
    while extraFound == 0:
        line = f.readline()
        if line.startswith('extraporation'):
            extraFound = 1
        else:
            values = line.split()
            xlist_for_each_n.append(float(values[3]))
            ylist_for_each_n.append(float(values[4]))
    line = f.readline()
    print(xlist_for_each_n)
    print(ylist_for_each_n)
    for Indx_of_txtfile in range(0, num_txtfiles):
        xlist_for_m = []
        ylist_for_m = []
        ForFound = 0
        while ForFound == 0:
            line = f.readline()
            if line.startswith('For'):
                ForFound = 1
            elif not line:
                break
            else:
                values = line.split()
                xlist_for_m.append(float(values[3]))
                ylist_for_m.append(float(values[4]))
        if (Indx_of_txtfile+1 <= num_txtfiles/2):
            ax=fig.add_subplot(num_txtfiles//2, 2, 2*(Indx_of_txtfile+1)-1)
        else:
            ax=fig.add_subplot(num_txtfiles//2, 2, 2*(Indx_of_txtfile+1)-num_txtfiles)
        ax.scatter(xlist_for_each_n[Indx_of_txtfile:],
                ylist_for_each_n[Indx_of_txtfile:], marker='x', clip_on=False,
                s=50, label="each n")
        ax.scatter(xlist_for_m[Indx_of_txtfile:], ylist_for_m[Indx_of_txtfile:],
                marker='+', clip_on=False, s=72, label="prediction")
        ax.set_xlim(0, 0.0001)
        ax.set_ylim(0, 0.06)
        hour = Indx_of_txtfile + 1
        ax.text(0.25, 0.8, 'n at %dh'%hour,
        transform=ax.transAxes, ha="right")

        if (Indx_of_txtfile+1 == 1):
           ax.legend()

        if (Indx_of_txtfile+1 == 1 or Indx_of_txtfile+1 == num_txtfiles//2 + 1):
           ax.set_yticks([0,0.02,0.04,0.06])
        else:
           ax.set_yticks([0,0.02,0.04])

        ax.tick_params(labelbottom=False)
        ax.tick_params(direction = "in")
        if (Indx_of_txtfile+1 ==  num_txtfiles//2  or Indx_of_txtfile+1 == num_txtfiles):
            ax.tick_params(labelbottom=True)
            ax.set_xlabel('1/m or 1/n')

        #plt.gca().ticklabel_format(style="sci", scilimits=(0,0), axis="y")
        plt.gca().ticklabel_format(style="sci", scilimits=(0,0), axis="x")
        if (Indx_of_txtfile == num_txtfiles//4 or Indx_of_txtfile == num_txtfiles//2 + num_txtfiles//4 ):
            ax.set_ylabel('1/(opt_wx*opt_wy)')
    plt.subplots_adjust(wspace=0.4, hspace=0.0)
    plt.show()

run()


