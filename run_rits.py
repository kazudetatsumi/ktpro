#!/usr/bin/env python
import numpy as np
import textfile
import os
import sys
import matplotlib.pyplot as plt
sys.path.append("/home/kazu/ktpro")
from gp_nrca import draw_sample
from rits_fit_kt import get_sim_spectrum
np.random.seed(126)


def run(lims_mean=[3.5, 3.6], lims_scale=[0.0, 0.05]):
    mean = np.random.uniform(low=lims_mean[0], high=lims_mean[1])
    scale = np.random.uniform(low=lims_scale[0], high=lims_scale[1])
    params = draw_sample(mean=mean, numsample=1, scale=scale)

    run_rits(params[0, :], 'AAAAAA')
    #plt.plot(params[0])
    #plt.ylim([0, np.max(params[0])*1.1])
    #plt.show()


def run_rits(params, string):
    for pidx, param in enumerate(params):
        os.system('cp rits_initial.inp.temp rits_initial.inp')
        _param = "{:.4f}".format(param)
        textfile.replace('rits_initial.inp', string, _param)
        x, y = get_sim_spectrum()
        if pidx == 0:
            bi2d_true = np.zeros((params.shape[0], x.shape[0]))
            bi2d_noisy = np.zeros((params.shape[0], x.shape[0]))
        bi2d_true[pidx] = y
        timescale = 2000
        #sensitivity = 0.01
        #_y = np.random.poisson(y*timescale*sensitivity*x)/(sensitivity*x)
        _y = np.random.poisson(y*timescale*x)/x
        bi2d_noisy[pidx] = _y
        #plt.plot(x,y)
    #os.system('python ./rits_fit_kt.py > /dev/null')
    plt.imshow(bi2d_true)
    #plt.plot(x,y*timescale)
    #plt.plot(x,_y)
    plt.show()


run()
