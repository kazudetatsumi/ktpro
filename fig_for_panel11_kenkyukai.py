#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys

argvs = sys.argv
argvs.append("orthotope_opt_20-70meV_ddscs_again2")

measurementtime = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 10, 24, 73, 120, 124, 211]) # minutes.
measurementtime = np.array([0.5,            2.5,    3.5,            10, 24, 73, 120,         ]) # minutes.
plt.rcParams['font.size'] = 14
plt.rcParams['font.family']= 'Arial'


def getdata(infile):
    f = open(infile, 'r')
    extraFound = 0
    line = f.readline()
    nn = []
    nqx = []
    nqy = []
    nqz = []
    nw = []
    counter = 0
    while extraFound == 0:
        line = f.readline()
        if line.startswith('extrapolation'):
            extraFound = 1
        else:
            counter += 1
            if counter == 1 or counter == 5 or counter == 7 or counter == 11 or counter == 12 or counter == 13 or counter == 14:
               values = line.split()
               nn.append(float(values[0]))
               nqx.append(float(values[1]))
               nqy.append(float(values[2]))
               nqz.append(float(values[3]))
               nw.append(float(values[4]))
    data = np.transpose(np.array([nn, nqx, nqy, nqz, nw]))
    #data[:, 0] = measurementtime
    return data


def run():
    filled = "/home/kazu/desktop/200312/for_cu_new/" + argvs[1] + "/filled/result.txt_vec"
    wo_cond = "/home/kazu/desktop/200312/for_cu_new/" + argvs[1] + "/result.txt_vec"
    condparam_09 = "/home/kazu/desktop/200312/for_cu_new/" + argvs[1] + "/condparam09/result.txt_vec"
    #condparam_01 = "/home/kazu/desktop/200312/for_cu_new/" + argvs[1] + "/condparam01/result.txt_vec"
    #filelist = [filled, wo_cond, condparam_09, condparam_01]
    filelist = [filled, wo_cond, condparam_09]
    wlist = ["$\mathregular{q_a}$", "$\mathregular{q_b}$","$\mathregular{q_c}$","$\mathregular{\omega}$",]
    for i, infile in enumerate(filelist):
        data = getdata(infile)
        if i == 0:
            all_data = np.zeros((len(filelist),data.shape[0],data.shape[1]))   # # of files, # of measurement times, # of props (nn, nqx, nqy, nqz, nw)
        all_data[i, :, :] = data
    fig = plt.figure(figsize=(8,8))
    fig.suptitle("optimal  bin-widths of 4D INS simulated orthotope data for #17714 " + argvs[1])
    for widx in range(1, 5):
        ax = fig.add_subplot(4, 1, widx)
        ax.scatter(all_data[0, :, 0], all_data[1, :,  widx], clip_on=False, s=100, edgecolor="black", label="w_mask", marker="x", facecolors="black", linewidth=2)
        #ax.scatter(all_data[0, :, 0], all_data[2, :,  widx], clip_on=False, s=400, edgecolor="black", label=r'$\alpha=0.9$', marker=".", facecolors="white")
        ax.scatter(all_data[0, :, 0], all_data[0, :,  widx], clip_on=False, s=100, edgecolor="black", label="no_mask", marker="o", facecolors="none", linewidth=2)
        ax.text(np.max(all_data[0, :, 0])*0.9, np.max(all_data[0:3, :, widx])*0.9, wlist[widx-1], fontsize=24)
        ax.tick_params(labelbottom=False)
        ax.set_xlim(np.max(all_data[0, :, 0])*(-0.01), np.max(all_data[0, :, 0])*1.05)
        #ax.set_xlim(0, 600000)
        ax.set_ylim(0, np.max(all_data[0:3, :, widx])*1.13)
        ax.set_yticks(range(0, int(np.max(all_data[0:3, :, widx]))+1))
        if widx == 4:
            ax.tick_params(labelbottom=True)
            ax.set_xlabel('total counts', fontsize=18)
        if widx == 2:
            ax.set_ylabel('bin width (1/step size)', fontsize=18)
    plt.subplots_adjust(wspace=0.4, hspace=0.0)
    #plt.scatter(all_data[0, :, 0], (np.sum(diff_data[7, :, 1:], axis=1)), clip_on=False, s=200, edgecolor="black", label="add_mask_lowE_"+r'$\alpha=0.9$', marker=".", facecolors="white")
    #plt.scatter(all_data[0, :, 0], (np.sum(diff_data[6, :, 1:], axis=1)), clip_on=False, s=60, edgecolor="black", label="add_mask_"+r'$\alpha=0.9$', marker=".", facecolors="black")
    #print(all_data[0,:,0])
    #print(np.sum(diff_data[3, :, 1:], axis=1))
    #plt.xlim(min(all_data[0, :, 0]), max(all_data[0, :, 0]))
    #plt.ylim(0, np.max(np.sum(diff_data[1:4, :, 1:], axis=2)))
    #plt.xlim(0, all_data[0, -1, 0])
    plt.gca().ticklabel_format(style="sci", scilimits=(0, 0), axis="x")
    #plt.xticks(position=(0.0, -0.03))
    #plt.yticks( np.arange(0, np.max(np.sum(diff_data[1:4, :, 1:], axis=2))+1) )
    #plt.xlabel('measurement time (arb. units)')
    #plt.ylabel('SUM of absolute differences in bin-width indexes')
    plt.legend(loc='lower left', fancybox=True, framealpha=0.5)
    ##non-ddscs
    #plt.savefig("difference_of_binwidths.pdf")
    ##ddscs
    plt.savefig("fig_binwidths_" + argvs[1] + ".pdf")
    plt.show()





    
run()
