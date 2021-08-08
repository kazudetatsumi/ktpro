#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/kazu/ktpro/')
import crosssection_of_4ddata_class as c4c
plt.rcParams.update({'font.size': 20})


def run():

    outfile = "/home/kazu/desktop/200312/for_cu_new/old_filled/filled_ddscs_kfki_debye/4h/eliminated_data.hdf5"
    orthotope_lims = np.array([[120, 172], [61, 145], [16, 53], [20, 70]])*1.0
    common_lims = np.array([[85, 88], [143, 146], [30, 39], [26, 32]])*1.0
    orgbinwidths = np.array([4, 4, 5, 3])
    binwidths = np.array([1, 1, 1, 1])*1.
    wholeranges = np.array([[-1.65, 4.1], [-2.1, 2.8], [-0.85, 0.9], [0, 40.5]]
                           )
    devs = np.array([1, 1, 1])
    cpos = np.array([27, 36, 7, 0])*orgbinwidths
    cpos = np.array([107, 145, 34, 0])

    pro = c4c.CROSS(outfile, orthotope_lims, wholeranges, devs, cpos,
                    hvlofs=False, binwidths=binwidths, common_lims=common_lims,
                    vivid=True)
    pro.create_fig(xdim=24, ydim=24)
    pro.get_data()
    pro.plot_crosssection(1, xyyx=True)

    pro.outfile = "/home/kazu/desktop/200701/all/filled_again_ddscs2/4m/eliminated_data.hdf5"
    pro.orthotope_lims = np.array([[114, 200], [69, 122], [11, 20],
                                  [81, 207]])*1.0
    orgbinwidths = np.array([3, 2, 2, 4])
    pro.binwidths = np.array([1, 1, 1, 1])*1.
    pro.wholeranges = np.array([[-0.675, 3.0875], [-0.925, 4.400],
                                [-0.8, 0.60], [-8.0, 36.4]])
    pro.common_lims = np.array([[172, 177], [56, 59], [14, 19], [105, 119]])*1.0
    pro.devs = np.array([0.75, 0.75, 0.75])
    pro.cpos = np.array([58, 39, 8, 10])*orgbinwidths

    pro.get_data()
    pro.plot_crosssection(2)

    pro.outfile = "/home/kazu/desktop/200903/all/filled_ddscs_kfki_debye/12m/eliminated_data.hdf5"
    pro.orthotope_lims = np.array([[0, 228], [0, 202], [0, 9], [150, 289]])*1.0
    orgbinwidths = np.array([2, 3, 2, 4])
    pro.binwidths = np.array([1, 1, 1, 1])*1.
    pro.wholeranges = np.array([[0.01, 2.31], [-0.67, 1.35], [-0.16, 0.18],
                                [-2.0, 21.20]])
    pro.common_lims = np.array([[146, 153], [114, 121], [2, 7], [188, 220]])*1.0
    pro.devs = np.array([1, 1, 1])
    pro.cpos = np.array([75, 39, 2, 6])*orgbinwidths
    pro.hvlofs = True

    pro.get_data()
    pro.plot_crosssection(3)
    plt.savefig("Fig1_advsoft_theo_crosssections.png")
    plt.show()


run()
