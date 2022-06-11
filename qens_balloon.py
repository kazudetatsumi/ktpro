#!/usr/bin/env python
# This script tries to form a densit esimation of qens
# using the bin width distribution of the adaptive KDE
# and the usual qens histogram wiht several corections.
import os
import sys
sys.path.append("/home/kazu/ktpro")
from qens_kde_results_odata_divided_by_idata_class\
    import odata_divided_by_idata as odbi
pwd = os.getcwd()


def getbandwidth(kf):
    proj = odbi(kf, 'dummy')
    return proj.read_pkl(kf)['y_ssvk']


def testrun():
    qrange = pwd.split('q')[1]
    qmin = float(qrange.split('-')[0])
    qmax = float(qrange.split('-')[1])
    sprefix = "/home/kazu/desktop/210108/Tatsumi/srlz/0000025/"
    allqf = sprefix + "run6204s.pkl"
    kf = "./qens_run6204united_kde_results_on_data_qsel.pkl"
    ky = getbandwidth(kf)
    print(ky[0].shape)


testrun()
