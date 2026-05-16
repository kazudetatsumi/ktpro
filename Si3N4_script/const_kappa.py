#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from matplotlib import rc
homedir = "/home/kazu/"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
Temp = 300
max_freq = 35
c = cdir + "noiso/kappa-m8810.hdf5"
s = sdir + "noiso/kappa-m8820.hdf5"
g = gdir + "noiso/kappa-m121212.hdf5"

def parse_gamma(filename,temp):
    f = h5py.File(filename,'r')
    temperature = f["temperature"].value
    i=0
    for t in temperature:
        if t == temp:
            tindex=i
        i += 1
    freq=f["frequency"]
    gamma=f["gamma"][tindex,]
    hcapa=f["heat_capacity"][tindex,]
    kappa=f["kappa"][tindex,]
    w=f["weight"]
    gvgv=f["gv_by_gv"]
    m=f["mesh"]
    kuc=f["kappa_unit_conversion"]
    return(freq,gamma,hcapa,w,gvgv,m,kuc,kappa)

def run():
    freqs,gammas,hcapas,ws,gvgvs,ms,kucs,kappas=parse_gamma(s,Temp)
    print freqs.shape
    print gammas.shape
    print hcapas.shape
    print ws.shape
    print gvgvs.shape


    pro1 = gvgvs[:,:,2] * gammas
    pro2 = pro1 * hcapas
    print pro2.shape
    pro3 = np.dot(ws, pro2) 
    print pro3.shape
    print np.sum(pro2)
    print kappas[2]


run()
