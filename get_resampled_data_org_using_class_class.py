#!/usr/bin/env python
# This script read the raw neutron count data of DNA without any corrections and resample the count data
# and apply the necessary corrections to deduce resampled QENS data corresponding to dobule differential 
# cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
import sys
sys.path.append("/home/kazu/ktpro")
#from get_resampled_data_mpi_class import Sget_qlist as sgq
from get_resampled_data_org_class_class import SSget_qlist as ssgq


#class SSget_qlist(sgq):
#    def __init__(self, runNo, TimeParam, qmin, qmax):
#        super().__init__(pklfile="run" + str(runNo) + "spectrab.pkl")
#        self.get_org_data("0.000025", runNo, TimeParam=TimeParam)
#        self.get_org_spectra(qmin, qmax)


def samplerun():
    qmin = 0.9
    qmax = 1.0
    #runNo = 6204
    #TimeParam = "10225.0, 12445.0"
    runNo = 6202
    #TimeParam = "8764.0, 10225.0"
    TimeParam = "-1.0/-1.0"
    print("TimeParam=", TimeParam)
    prj = ssgq(runNo, TimeParam, qmin, qmax)
    prj.save_pkl()


samplerun()
