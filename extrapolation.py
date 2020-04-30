#!/usr/bin/env python
# This script gathers the Cn.hdf5 files generated by optbinwidth4D_wholefort.py  for several INS data sets with different measurement times
# and summarizes the statistics, plot the 1/m - 1/delta (extrapolation)  and 1/n - 1/delta (optimization for the data).
# 2020/02/15 Kazuyoshi TATSUMI
# a bug in the extrapolation routine was fixed.
# 2020/04/26
# usage:  extrapolation.py num_of_hdfs
import numpy as np
import h5py
import sys


def get_optid(m, n, Cn, kave, delta):
    ex = np.zeros((delta.shape[0], delta.shape[1], delta.shape[2], delta.shape[3]))
    # This is a bag, originated from an old purely python script of the extrapolation method.
    #ex[1:, 1:, 1:, 1:] = (1/m - 1/n) * kave[1:, 1:, 1:, 1:] / (delta[1:, 1:, 1:, 1:]**2*n)
    #ex[0, :, :, :] = 0.0
    #ex[:, 0, :, :] = 0.0
    #ex[:, :, 0, :] = 0.0
    #ex[:, :, :, 0] = 0.0
    ex[:, :, :, :] = (1/m - 1/n) * kave[:, :, :, :] / (delta[:, :, :, :]**2*n)
    Cm = ex + Cn
    optid = np.unravel_index(np.argmin(Cm, axis=None), Cm.shape)
    return optid


def unite_hdf5(num_hdfs, head):
    for i in range(0, num_hdfs):
        outfile = head + str(i+1) + "h/Cn.hdf5"
        f = h5py.File(outfile, 'r')
        if i == 0:
            tCn = np.zeros((f["Cn"].shape[0], f["Cn"].shape[1], f["Cn"].shape[2], f["Cn"].shape[3], num_hdfs))
            tdelta = np.zeros((f["delta"].shape[0], f["delta"].shape[1], f["delta"].shape[2], f["delta"].shape[3], num_hdfs))
            tkave = np.zeros((f["kave"].shape[0], f["kave"].shape[1], f["kave"].shape[2], f["kave"].shape[3], num_hdfs))
        tCn[:, :, :, :, i] = f["Cn"]
        tdelta[:, :, :, :, i] = f["delta"]
        tkave[:, :, :, :, i] = f["kave"]
        datafile = head + str(i+1) + "h/eliminated_data.hdf5"
        g = h5py.File(datafile, 'r')
        if i == 0:
            tTotInt = np.zeros((num_hdfs))
        tTotInt[i] = np.sum(g["data4"])
    outfile = head + "Cn_all.hdf5"
    with h5py.File(outfile, 'w') as hf:
        hf.create_dataset('tCn', data=tCn)
        hf.create_dataset('tdelta', data=tdelta)
        hf.create_dataset('tkave', data=tkave)
        hf.create_dataset('tTotInt', data=tTotInt)


def get_stat_via_unitehdf5(head):
    outfile = head + "Cn_all.hdf5"
    f = h5py.File(outfile, 'r')
    tCn = f["tCn"]
    tdelta = f["tdelta"]
    tkave = f["tkave"]
    tTotInt = f["tTotInt"]
    return tCn, tdelta, tkave, tTotInt


def run():
    num_hdf = int(sys.argv[1])
    head = "./"
    #num_hdf = 10
    unite_hdf5(num_hdf, head)
    tCn, tdelta, tkave, tTotInt = get_stat_via_unitehdf5(head)
    outfile = "result.txt_vec"
    result = []

    for hid in range(1, num_hdf + 1):
        print("Extrapolation with n = ", tTotInt[hid - 1])
        l_result = []
        Cn = tCn[:, :, :, :, hid - 1]
        kave = tkave[:, :, :, :, hid - 1]
        delta = tdelta[:, :, :, :, hid - 1]
        n = tTotInt[hid - 1]
        Cn = Cn / (n**2)   # This is according to the Cn in NeCo(2007)
        for m in tTotInt:  # Here we extrapolate the optimization from the total intensity of the data, n, toward m taken as each of total intensities acquired with different measurement times.
            optid = get_optid(m, n, Cn, kave, delta)
            l_result.append([optid[0]+1, optid[1]+1, optid[2]+1, optid[3]+1, m]) # Because numpy index starts from 0, here we add 1 for treating the actual amounts of the width.
            if m == tTotInt[-1]:
                max_optid = optid
        for mn_ratio in np.arange(1.1, 20.1, 0.1):
            m = mn_ratio * tTotInt[hid - 1]
            optid = get_optid(m, n, Cn, kave, delta)
            if 1/(np.prod(np.array(optid)+1)*1.0) > 1/(np.prod(np.array(max_optid)+1)*1.0):
                max_optid = optid
                l_result.append([optid[0]+1, optid[1]+1, optid[2]+1, optid[3]+1, m]) # Because numpy index starts from 0, here we add 1 for treating the actual amounts of the width.
        result.append(l_result)
    with open(outfile, mode='w') as f:
        f.write("optimization results for 17714 data: n, optid_x, optid_y, optid_z, optid_w, 1/n, 1/(optidnx_x*optid_y*optid_z*optid_w) \n")
        for hid in range(0, num_hdf):
            f.write("%e %d %d %d %d %e %e\n" %
             (
              result[hid][hid][4], result[hid][hid][0], result[hid][hid][1], result[hid][hid][2], result[hid][hid][3], 1/(result[hid][hid][4]*1.0),
              1/(np.prod(np.array(result[hid][hid][0:4]))*1.0)
             )
            )
        f.write("extrapolation results for 17714 data: m, optid_x, optid_y, optid_z, optid_w, 1/m, 1/(optidnx_x*optid_y*optid_z*optid_w) \n")
        for hid in range(0, num_hdf):
            f.write("For n = %e \n" % tTotInt[hid])
            for mid in range(0, len(result[hid])):
                f.write("%e %d %d %d %d %e %e\n" %
                 (
                  result[hid][mid][4], result[hid][mid][0], result[hid][mid][1], result[hid][mid][2], result[hid][mid][3],
                  1/(result[hid][mid][4]*1.0), 1/(np.prod(np.array(result[hid][mid][0:4]))*1.0)
                 )
                )


run()
