#!/usr/bin/env python
# This script read the raw neutron count data of DNA without any corrections and resample the count data
# and apply the necessary corrections to deduce resampled QENS data corresponding to dobule differential 
# cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
from get_resampled_data_mpi_class import Sget_qlist as sgq
import datetime
import numpy as np
from mpi4py import MPI
import copy

def run():
    runNo = 6204
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    prj = sgq(pklfile="run" + str(runNo) + "spectrab.pkl")
    #prj.get_org_data("0.000025")
    #prj.get_org_data("0.000025", runNo, TimeParam="0.0, 1460.7")
    prj.get_org_data("0.000025", runNo)
    if rank == 0:
        print(datetime.datetime.now(), 'org_data ended')
    #prj.get_org_intensity_array()
    prj.get_all_data()
    if rank == 0:
        print("chk dataset['omega'] shape from ecm of self.DAT by prj.get_org_data:", prj.dataset['omega'].shape)
        print(prj.dataset['omega'][0,0,0:10])
    # The object name of prj.intensity instead of prj.dataset['intenity'] is needed for  prj.get_boot_strap_sampled_spectra.
    prj.intensity = np.array(prj.dataset['intensity'])
    # Saving the present prj.dataset as prj.datasetnocorr to use it in the latter half of prj.get_boot_strap_sampled_spectra.
    prj.datasetnocorr = copy.deepcopy(prj.dataset)
    print(datetime.datetime.now(), 'org_intensity_array ended')
    nbs = 4
    qmin = 0.55
    qmax = 0.70
    prj.get_boot_strap_sampled_spectra(nbs, qmin, qmax, restart=False, wnocorr=True)
    if rank == 0:
        print(datetime.datetime.now(), 'boot_strap_sampling ended')
        prj.save_pkl()
        print(datetime.datetime.now(), 'boot_strap_sampled_spectra saved')
    #intensities = np.sort(np.unique(prj.intensity.flatten()))
    #print(intensities[0:5])
    #print(intensities[-5:-1])


def run_org():
    runNo = 6202
    prj = sgq(pklfile="run" + str(runNo) + "spectrab.pkl")
    prj.get_org_data("0.000025", runNo)
    prj.get_qemap()
    prj.get_all_sdata()
    prj.spect(0.55, 0.70, isplot=True)
    prj.spectrab = np.zeros((1, 2, prj.ene.shape[0]))
    prj.spectrab[0, 0, :] = prj.ene
    prj.spectrab[0, 1, :] = prj.spectra
    prj.save_pkl()


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
#run_org()
