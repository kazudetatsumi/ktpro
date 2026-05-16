#!/usr/bin/env python
# This script read the raw neutron count data of DNA without any corrections and resample the count data
# and apply the necessary corrections to deduce resampled QENS data corresponding to dobule differential 
# cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
from get_resampled_data_mpi_class import Sget_qlist as sgq
import datetime
import numpy as np
from mpi4py import MPI

def run():
    runNo = 6202
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    prj = sgq(pklfile="run" + str(runNo) + "spectrab.pkl")
    #prj.get_org_data("0.000025")
    prj.get_org_data("0.000025", runNo, TimeParam="0.0, 1460.7")
    print(datetime.datetime.now(), 'org_data ended')
    prj.get_org_intensity_array()
    print(datetime.datetime.now(), 'org_intensity_array ended')
    nbs = 96
    #nbs = 8
    qmin = 0.55
    qmax = 0.70
    prj.get_boot_strap_sampled_spectra(nbs, qmin, qmax, restart=False)
    print(datetime.datetime.now(), 'boot_strap_sampled_spectra ended', rank)
    if rank == 0:
        prj.save_pkl()
    prj.save_pkl()
    #intensities = np.sort(np.unique(prj.intensity.flatten()))
    #print(intensities[0:5])
    #print(intensities[-5:-1])


def check():
    prj = sgq(pklfile="run6202spectraorg.pkl")
    prj.load_pkl()
    for inb in range(prj.spectrab.shape[0]):
        plt.plot(prj.spectrab[inb, 0, :], prj.spectrab[inb, 1, :], lw=0.2)
    prj = Sget_qlist(pklfile="run6202spectrab.pkl")
    prj.load_pkl()
    for inb in range(prj.spectrab.shape[0]):
        plt.plot(prj.spectrab[inb, 0, :], prj.spectrab[inb, 1, :], lw=0.2)
    #plt.xlim(-0.025, 0.075)
    plt.ylim(.0, 8.)
    plt.show()


run()
#check()
