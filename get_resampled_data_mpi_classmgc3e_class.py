#!/usr/bin/env python
# This script read the raw neutron count data of DNA without any corrections and resample the count data
# and apply the necessary corrections to deduce resampled QENS data corresponding to dobule differential
# cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
# get_resampled_data_using_class_w_uncorr_mpi_run6202.py is altered as a class over get_resampled_data_mpi_class.
# Kazuyoshi TATSUMI 2023/05/30
import sys
import numpy as np
sys.path.append("/home/kazu/ktpro")
from get_resampled_data_mpi_classmgc2e import Sget_qlist as sgq


class SSget_qlist(sgq):
    def __init__(self, runNo, nbs, TimeParam, qmin, qmax,
                 maskfile="maskTY140218ForAfterRun52.txt"):
        super().__init__(pklfile="run" + str(runNo) + "spectrab.pkl",
                         randfile='randomstates.pkl.' + str(runNo) + '.30000')
        self.get_org_data("0.000025", runNo, TimeParam=TimeParam)
        self.get_all_data3()
        self.orgintensityshape = self.dataset['intensity'].shape
        self.orgintensity1d = np.array(self.dataset['intensity'].flatten())
        self.get_boot_strap_sampled_spectra(nbs, qmin, qmax, wnocorr=True,
                                            maskfile=maskfile)


def sample_run(runNo, nbs, TimeParam, qmin, qmax):
    from mpi4py import MPI
    import datetime
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    prj = SSget_qlist(runNo, nbs, TimeParam, qmin, qmax)
    if rank == 0:
        print(datetime.datetime.now(), 'boot_strap_sampling ended')
        prj.save_pkl()
        print(datetime.datetime.now(), 'boot_strap_sampled_spectra saved')


#sample_run(6202, 2,  "-1.0/-1.0", 0.7, 0.8)
