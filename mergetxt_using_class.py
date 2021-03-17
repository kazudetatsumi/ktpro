#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
import mergetxt_class as mgc


def run():
    yfile = "/home/kazu/desktop/200522/Ei24/fineq/add_random_mask/maxfrac.txt"
    xfile = "/home/kazu/desktop/200522/Ei24/fineq/nbins.txt"
    proj = mgc.Merge_Txt(xfile, yfile)
    proj.merge()
    proj.create_fig()
    proj.plotter(r'#3($\alpha$=0)', 'blue')

    proj.yfile = "/home/kazu/desktop/200522/Ei24/fineq/add_random_mask/maxfrac_condparam07.txt"
    proj.xfile = "/home/kazu/desktop/200522/Ei24/fineq/condparam07/nbins.txt"
    proj.merge()
    proj.plotter(r'#3($\alpha$=0.7)', 'k')

    proj.yfile = "/home/kazu/desktop/200522/Ei42/veryfineq/add_random_mask/maxfrac.txt"
    proj.xfile = "/home/kazu/desktop/200522/Ei42/veryfineq/nbins.txt"
    proj.merge()
    proj.plotter(r'#2($\alpha$=0)', 'limegreen')

    proj.yfile = "/home/kazu/desktop/200522/Ei42/veryfineq/add_random_mask/maxfrac_condparam09.txt"
    proj.xfile = "/home/kazu/desktop/200522/Ei42/veryfineq/condparam_09/nbins.txt"
    proj.merge()
    proj.plotter(r'#2($\alpha$=0.9)', 'darkgreen')

    proj.yfile = "/home/kazu/desktop/200204/fine/hourbyhour/add_random_mask/maxfrac.txt"
    proj.xfile = "/home/kazu/desktop/200204/fine/hourbyhour/ortho_opt_without_mask/nbins.txt"
    proj.merge()
    proj.plotter(r'#1($\alpha$=0)', 'red')

    proj.yfile = "/home/kazu/desktop/200204/fine/hourbyhour/add_random_mask/maxfrac_condparam09.txt"
    proj.xfile = "/home/kazu/desktop/200204/fine/hourbyhour/ortho_opt_without_mask/condparam09/nbins.txt"
    proj.merge()
    proj.plotter(r'#1($\alpha$=0.9)', 'brown')
    plt.xlabel('number of bins')
    plt.ylabel('upper limit for additional mask volume fraction (%)')
    plt.tick_params(top=True, right=True, direction='in', which='both')

    plt.legend()

    plt.show()


run()
