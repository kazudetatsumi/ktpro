#!/usr/bin/env python
# This script read the raw neutron count data of DNA without any corrections and resample the count data
# and apply the necessary corrections to deduce resampled QENS data corresponding to dobule differential 
# cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
from get_resampled_data_mpi_class import Sget_qlist as sgq
import datetime
import numpy as np
from mpi4py import MPI


def run_org(runNo, TimeParam):
    prj = sgq(pklfile="run" + str(runNo) + "spectraorg.pkl")
    prj.get_org_data("0.000025", runNo, TimeParam=TimeParam)
    prj.get_qemap()
    prj.get_all_sdata()
    prj.spect(0.55, 0.70, prj.dataset, isplot=True)
    prj.spectrab = np.zeros((1, 2, prj.ene.shape[0]))
    prj.spectrab[0, 0, :] = prj.ene
    prj.spectrab[0, 1, :] = prj.spectra
    prj.get_all_data()
    prj.spect(0.55, 0.70, prj.dataset, isplot=True)
    prj.spectrab = np.concatenate((prj.spectrab, prj.spectra.reshape(1, 1, -1)), axis=1)
    prj.save_pkl()


runNo = 6202
TimeParam = "8764.0, 10225.0"
print("TimeParam=", TimeParam)
run_org(runNo, TimeParam)