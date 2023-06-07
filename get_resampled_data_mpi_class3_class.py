#!/usr/bin/env python
# This script read the raw neutron count data of DNA without any corrections and resample the count data
# and apply the necessary corrections to deduce resampled QENS data corresponding to dobule differential
# cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
# get_resampled_data_using_class_w_uncorr_mpi_run6202.py is altered as a class over get_resampled_data_mpi_class.
# Kazuyoshi TATSUMI 2023/05/30
import sys
sys.path.append("/home/kazu/ktpro")
from get_resampled_data_mpi_class3 import Sget_qlist as sgq
import numpy as np
#from mpi4py import MPI
#import datetime


class SSget_qlist(sgq):
    def __init__(self, runNo, nbs, TimeParam, qmin, qmax):
        super().__init__(pklfile="run" + str(runNo) + "spectrab.pkl")
        self.get_org_data("0.000025", runNo, TimeParam=TimeParam)
        self.get_all_data()
        self.orgintensityshape = self.dataset['intensity'].shape
        self.orgintensity1d = np.array(self.dataset['intensity'].flatten()).astype('uint32')
        self.get_boot_strap_sampled_spectra(nbs, qmin, qmax, wnocorr=True)


def sample_run(runNo, nbs, TimeParam, qmin, qmax):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    prj = SSget_qlist(runNo, nbs, TimeParam, qmin, qmax)
    if rank == 0:
        print(datetime.datetime.now(), 'boot_strap_sampling ended')
        prj.save_pkl()
        print(datetime.datetime.now(), 'boot_strap_sampled_spectra saved')


#sample_run(6202, 2,  "-1.0/-1.0", 0.7, 0.8)