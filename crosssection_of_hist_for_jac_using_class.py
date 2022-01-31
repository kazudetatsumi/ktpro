#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/kazu/ktpro/')
import crosssection_of_4ddata_class as c4c
plt.rcParams.update({'font.size': 20})



def run():
    head = "/home/kazu/desktop/"
    tail = "hist.hdf5"
    outfile = head + "200204/fine/hourbyhour/10h/" + tail

    orthotope_lims = np.array([[120, 172], [61, 145], [16, 53], [20, 70]])*1.0
    binwidths = np.array([4, 4, 5, 3])*1.
    wholeranges = np.array([[-1.65, 4.1], [-2.1, 2.8], [-0.85, 0.9], [0, 40.5]]
                           )
    devs = np.array([1, 1, 1])
    cpos = np.array([27, 36, 7, 0])
    #cpos = np.array([107, 145, 34, 0])

    pro = c4c.CROSS(outfile, orthotope_lims, wholeranges, devs, cpos,
                    hvlofs=False, binwidths=binwidths, showorthob=False, vivid=True)
    pro.create_fig(xdim=24, ydim=24)
    pro.get_data()
    pro.plot_crosssection(1, xyyx=True)

    pro.outfile = head + "200522/Ei42/veryfineq/14m/" + tail
    pro.orthotope_lims = np.array([[114, 200], [69, 122], [11, 20],
                                  [81, 207]])*1.0
    pro.binwidths = np.array([3, 2, 2, 4])*1.
    pro.wholeranges = np.array([[-0.675, 3.0875], [-0.925, 4.400],
                                [-0.8, 0.60], [-8.0, 36.4]])
    pro.cpos = np.array([58, 39, 8, 10])

    pro.get_data()
    pro.plot_crosssection(2)

    pro.outfile = head + "200522/Ei42/veryfineq/195s/" + tail
    pro.orthotope_lims = np.array([[114, 200], [69, 122], [11, 20],
                                  [81, 207]])*1.0
    pro.binwidths = np.array([4, 4, 2, 4])*1.
    pro.wholeranges = np.array([[-0.675, 3.0875], [-0.925, 4.400],
                                [-0.8, 0.60], [-8.0, 36.4]])
    pro.cpos = np.array([43, 19, 8, 10])

    pro.get_data()
    pro.plot_crosssection(3)
    plt.savefig("fig_jac_hist_with_labels.png")
    plt.subplots_adjust(wspace=0.1, hspace=0.37)
    plt.show()


run()
