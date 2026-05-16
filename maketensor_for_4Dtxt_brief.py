#!/usr/bin/env python
import numpy as np
import h5py
import sys
import os


def get4ddata(txtfile, nxyzw):
    #   genformtxt is too slow and requires too large memory for the purpose of
    #   extracting only one column. So, I comment out and use a for-loop without readlines().
    # data = np.genfromtxt(f, dtype=float, comments='#', delimiter=',')
    # intensity = data[:, 4]  # data[:, 4] are mask if the values are 1e*30
    intensity = []
    for line in open(txtfile):
            intensity.append(float(line[:-1].split(',')[-2]))
    intensity = np.array(intensity[4:])
    nx = nxyzw[0]
    ny = nxyzw[1]
    nz = nxyzw[2]
    nw = nxyzw[3]

    if intensity.shape[0] != nx*ny*nz*nw:
        print("number of data elements is different from your expected")
        sys.exit()
    karr2 = np.reshape(intensity, (nw, nx, ny, nz))
    karr2 = np.transpose(karr2, (1, 2, 3, 0))
    if abs(np.sum(karr2) - np.sum(intensity)) > 1:
        print("total intensity is not preserved!!",
              np.sum(karr2), np.sum(intensity))
        sys.exit()
    condition = karr2 < 0.9e+30 # mask if the values are 1e*30
    #karr2 = karr2*condition/12.5   # 12.5 means one neutron count in the txtfile of Ei24
    #karr2 = karr2*condition/2.0   # 2.0 means one neutron count in the txtfile of Ei24 with 0.025 0.025 0.025 0.5
    karr2 = karr2*condition  
    return karr2, condition


def gen_hdf5(num_txtfiles, head):
    for i in range(0, num_txtfiles):
        txtfile = head + "Output4D_00_" + str((i+1)*60) + ".txt"
        print(txtfile)
        if os.path.exists(txtfile):
           outfile = head + "Output4D_00_" + str((i+1)*60) + ".hdf5"
           hfile = "head4line"
           os.system("head --line=4 " + txtfile + " > " + hfile)
           nxyzw = []
           for line in open('head4line', "r").readlines():
              nxyzw.append(int(float(line[:-1].split(',')[-1])))
           os.system("rm " + hfile)
           data4, condition = get4ddata(txtfile, nxyzw)
           #data4, condition = get4ddata_2(txtfile, nxyzw)
           with h5py.File(outfile, 'w') as hf:
              hf.create_dataset('data4', data=data4)
              hf.create_dataset('condition', data=condition)


def gen_hdf5_seconds(seconds, head):
    for second in seconds:
        txtfile = head + "Output4D_00_" + second + ".txt"
        print(txtfile)
        if os.path.exists(txtfile):
           outfile = head + "Output4D_00_" + second + ".hdf5"
           hfile = "head4line"
           os.system("head --line=4 " + txtfile + " > " + hfile)
           nxyzw = []
           for line in open('head4line', "r").readlines():
              nxyzw.append(int(float(line[:-1].split(',')[-1])))
           os.system("rm " + hfile)
           data4, condition = get4ddata(txtfile, nxyzw)
           #data4, condition = get4ddata_2(txtfile, nxyzw)
           with h5py.File(outfile, 'w') as hf:
              hf.create_dataset('data4', data=data4)
              hf.create_dataset('condition', data=condition)


def run():
    num_txtfiles = 26
    num_txtfiles = 1
    head = "./"
    seconds = ["000", "120", "360"]
    seconds = ["000"]
    #gen_hdf5(num_txtfiles, head)
    gen_hdf5_seconds(seconds, head)


run()
