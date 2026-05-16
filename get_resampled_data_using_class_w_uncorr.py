#!/usr/bin/env python
# This script read the raw neutron count data of DNA without any corrections and resample the count data
# and apply the necessary corrections to deduce resampled QENS data corresponding to dobule differential 
# cross-sections.
# Kazuyoshi TATSUMI 2023/02/15
import sys
sys.path.append("/home/kazu/ktpro")
from get_resampled_data_class import Sget_qlist as sgq
import datetime
import numpy as np
import copy

def run(runNo, nbs, restart, TimeParam, qmin, qmax):
    prj = sgq(pklfile="run" + str(runNo) + "spectrab.pkl")
    prj.get_org_data("0.0025", runNo, TimeParam=TimeParam)
    print(datetime.datetime.now(), 'org_data ended')
    prj.get_all_data()
    print("chk dataset['omega'] shape from ecm of self.DAT by\
          prj.get_org_data:", prj.dataset['omega'].shape)
    print(prj.dataset['omega'][0, 0, 0:10])
    # The object name of prj.intensity instead of prj.dataset['intenity'] is
    # needed for prj.get_boot_strap_sampled_spectra.
    prj.intensity = np.array(prj.dataset['intensity'])
    # Saving the present prj.dataset as prj.datasetnocorr to use it in the
    # latter half of prj.get_boot_strap_sampled_spectra.
    prj.datasetnocorr = copy.deepcopy(prj.dataset)
    print(datetime.datetime.now(), 'org_intensity_array ended')
    prj.get_boot_strap_sampled_spectra(nbs, qmin, qmax, restart=restart,
                                       wnocorr=True)
    print(datetime.datetime.now(), 'boot_strap_sampling ended')
    prj.save_pkl()
    print(datetime.datetime.now(), 'boot_strap_sampled_spectra saved')
    #intensities = np.sort(np.unique(prj.intensity.flatten()))
    #print(intensities[0:5])
    #print(intensities[-5:-1])


def run_org(runNo, TimeParam, tidx):
    prj = sgq(pklfile="run" + str(runNo) + "spectraorg" + str(tidx) + ".pkl")
    prj.get_org_data("0.000025", runNo, TimeParam=TimeParam)
    prj.get_qemap()
    prj.get_all_sdata()
    prj.spect(0.55, 0.70, isplot=True)
    prj.spectrab = np.zeros((1, 2, prj.ene.shape[0]))
    prj.spectrab[0, 0, :] = prj.ene
    prj.spectrab[0, 1, :] = prj.spectra
    prj.save_pkl()


def run_orgs(runNo, TimeParams):
    prj = sgq(pklfile="run" + str(runNo) + "spectraorgs.pkl")
    for tidx, TimeParam in enumerate(TimeParams):
        prj.get_org_data("0.000025", runNo, TimeParam=TimeParam)
        prj.get_qemap()
        prj.get_all_sdata()
        prj.spect(0.55, 0.70, isplot=True)
        if tidx == 0:
            prj.spectrab = np.zeros((len(TimeParams), 2, prj.ene.shape[0]))
        prj.spectrab[tidx, 0, :] = prj.ene
        prj.spectrab[tidx, 1, :] = prj.spectra
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


#run(6204, 35, True,  "10225.0, 12445.0")
run(6204, 1, False,  "-1.0/-1.0", 0.9, 1.0)
#run(6202, 35,  True, "8764.0, 10225.0")
#check()
#run_orgs(6204, ["8881.0, 11101.0", "11101.0, 13321.0", "13321.0, 15542.0",
#                "15542.0, 17762.0"])
#run_orgs(6202, ["5842.0, 7303.0", "7303.0, 8764.0", "8764.0, 10225.0",
#                "10225.0, 11686.0"])
