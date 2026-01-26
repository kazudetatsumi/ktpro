#!/usr/bin/env python
import matplotlib.pyplot as plt
import sys
sys.path.append("/home/kazu/ktpro")
from sf_bi_spike_at_boudary_phantom import compare_images4_2d_phantom,\
        compare_images4_d_phantom
from sf_bi_spike_at_boudary import compare_images4_2d, compare_images4_d


def compare_2d():
    fig = plt.figure(figsize=(11, 10))
    subfigs = fig.subfigures(2, 1, hspace=0.0)
    for sf, plotter in zip(subfigs, [compare_images4_2d_phantom,
                                     compare_images4_2d]):
        plotter(sf)
    plt.savefig('compare_spikes_2d.eps')
    plt.show()


def compare_1d():
    fig = plt.figure(figsize=(4, 6))
    subfigs = fig.subfigures(2, 1, hspace=0.0)
    for sf, plotter in zip(subfigs, [compare_images4_d_phantom,
                                     compare_images4_d]):
        plotter(sf)
    plt.savefig('compare_spikes_1d.eps')
    plt.show()


compare_2d()
