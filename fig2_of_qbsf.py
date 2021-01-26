#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys

#argvs = sys.argv
params = {'mathtext.default': 'regular'}
plt.rcParams.update(params)

plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'


def getdata(infile):
    nqx = []
    nqy = []
    nqz = []
    nw = []
    nn = []
    with open(infile, 'r') as f:
        for line in f:
            values = line.split()
            nqx.append(float(values[0]))
            nqy.append(float(values[1]))
            nqz.append(float(values[2]))
            nw.append(float(values[3]))
            nn.append(float(values[4]))
    data = np.transpose(np.array([nn, nqx, nqy, nqz, nw]))
    return data


def run():
    wo_cond = "/home/kazu/desktop/200312/for_cu_new/old_orthotope" +\
              "/orthotope_ddscs_again2/expt_orthotope_bd/" + \
              "result_only_extrapolate"
    condparam_09 = "/home/kazu/desktop/200312/for_cu_new/old_orthotope" +\
                   "/orthotope_ddscs_again2/expt_orthotope_bd/condparam09/" +\
                   "result_only_extrapolate"
    #wo_cond = "/home/kazu/desktop/200701/orthotope_again_ddscs/" +\
    #          "result_only_extrapolate"
    #condparam_09 = "/home/kazu/desktop/200701/orthotope_again_ddscs/" +\
    #               "condparam09/result_only_extrapolate"
    dqx = 0.025
    dqy = 0.025
    dqz = 0.025
    de = 0.5
    #dqx = 0.0125
    #dqy = 0.025
    #dqz = 0.050
    #de = 0.2
    filelist = [wo_cond, condparam_09]
    wlist = ["$q_x$", "$q_y$", "$q_z$", "$\omega$"]
    deltas = [dqx, dqy, dqz, de]
    for i, infile in enumerate(filelist):
        data = getdata(infile)
        # convert the unist from step to rlu and meV
        data[:, 1:] = data[:, 1:]*deltas
        if i == 0:
            # the shape of the following array:
            # # of files, # of measure  ment times, # of props
            #                                         (nn, nqx, nqy, nqz, nw)
            all_data = np.zeros((len(filelist), data.shape[0], data.shape[1]))
        all_data[i, :, :] = data
    fig = plt.figure(figsize=(6, 12.5))
    fig.suptitle("extrapolated optimal bin-widths of 4D INS orthotope data" +
                 " for phantom of #17714")
    for widx in range(1, 5):
        dp = deltas[widx -1]*2.0
        ax = fig.add_subplot(4, 1, widx)
        #ax.scatter(np.log10(all_data[0, :, 0]), all_data[0, :,  widx],
        #            clip_on=False, s=70, edgecolor="none", label="wo_cond",
        #            marker="|", facecolors="black")
        #ax.scatter(np.log10(all_data[0, :, 0]), all_data[1, :,  widx],
        #           clip_on=False, s=10, edgecolor="none",
        #           label=r'$\alpha=0.9$', marker="|", facecolors="gray")
        ax.plot(np.log10(all_data[0, :, 0]), all_data[0, :,  widx],
                clip_on=False, linestyle="dotted", label="wo_cond",
                marker="x", color='k')
        ax.plot(np.log10(all_data[0, :, 0]), all_data[1, :,  widx],
                clip_on=False, linestyle="dotted", label=r'$\alpha=0.9$',
                marker=".", color='k')
        ax.text(np.max(np.log10(all_data[0, :, 0]))*0.94,
                np.max(all_data[0:3, :, widx])*0.9, wlist[widx-1], size=20)
        ax.tick_params(labelbottom=False, direction='in', top=True, right=True)
        ax.set_ylim(0, (np.max(all_data[0:2, :, widx])//dp+1)*dp)
        ax.set_yticks(np.arange(0,  (np.max(all_data[0:2, :, widx])//dp+2)*dp, dp))
        xmargin = (np.max(np.log10(all_data[0, :, 0])) - np.min(np.log10(all_data[0, :, 0])))*0.01
        ax.set_xlim(np.min(np.log10(all_data[0, :, 0]))-xmargin, np.max(np.log10(all_data[0, :, 0]))+xmargin)
        print('xmargin', xmargin*100.0)
        if widx == 4:
            ax.tick_params(labelbottom=True)
            ax.set_xlabel('log10(total count)')
            ax.set_ylabel('bin width (meV)')
        if widx == 2:
            ax.set_ylabel('bin width (rlu)')
    plt.subplots_adjust(wspace=0.4, hspace=0.0)
    plt.legend(loc="upper center")
    plt.savefig("fig2_of_qbsf.pdf")
    plt.show()

    
run()
