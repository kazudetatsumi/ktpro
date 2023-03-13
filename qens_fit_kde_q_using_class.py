#!/usr/bin/env python
import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append("/home/kazu/ktpro")
from qens_fit_class import qens_fit as qf


def fit(devf, tf, variables, elim):
    prj = qf(devf, tf, elim, showplot=True, leastsq=False)
    prj.icorr()
    prj.preprocess(doicorr=True)
    prj.bg = 0.
    prj.optimize(variables=variables, figname='test.png')


def run():
    prefix = "./"
    elim = [-0.03, 0.07]
    variables = [0.4, 0.026, 0.5, 0.0075, 0.22, 0.01]
    variables = [0.4, 0.026, 0.22, 0.01]
    devf = prefix + "qens_kde_o_divided_by_i_6204.pkl"
    tf = prefix + "qens_kde_o_divided_by_i_6202.pkl"
    fit(devf, tf, variables, elim)


run()
