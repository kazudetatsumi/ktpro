#!/usr/bin/env python
import numpy as np


def alpha(dataname, expt=True):
    vlims = np.array([[5.155200e+04, 5.189990e+05],
                      [1.140130e+05, 5.906150e+05],
                      [7.821000e+03, 9.685300e+04]])
    tcountcommon_expt = np.array([123.0, 483.0, 1885.0])
    tcount_theo = np.array([564026.0, 702056.0, 103368.0])
    tcountcommon_theo = np.array([343.0, 1424.0, 3208.0])
    alphas_expt = tcountcommon_expt/vlims[:, 1]
    alphas_theo = tcountcommon_theo/tcount_theo
    datanametoindex = {'17714': 0, 'Ei42': 1, 'Ei24': 2}
    if expt:
        alphas = alphas_expt
    else:
        alphas = alphas_theo
    return alphas[datanametoindex[dataname]]

