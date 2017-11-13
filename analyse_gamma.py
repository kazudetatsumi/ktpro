#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
homedir = "/home/kazu/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
Temp = 300
nbins = 300
#y_max = 0.12
max_freq = 15
fs = 9
s = sdir + "noiso/kappa-m141432.noiso.hdf5"
sv = sdir + "qpoints.hdf5"


def parse_gamma(filename,temp):
    f = h5py.File(filename,'r')
    temperature = f["temperature"].value
    i=0
    for t in temperature:
        if t == temp:
            tindex=i
        i += 1
    gamma_k=f["gamma"][tindex,]
    omega_k=f["frequency"][:,:]
    qpoint_k=f["qpoint"][:]
    #print np.array(omega_k).shape
    print np.array(gamma_k).shape
    #print np.array(qpoint_k).shape
    #print omega_k[0,0:32]

    return(omega_k,gamma_k,qpoint_k)

def parse_eigenvec(filename):
    #freqs = []
    #mode_prop = []
    f = h5py.File(filename,'r')
    omega_q=f["frequency"][:,:]
    qpoint_q=f["qpoint"][:]
    eigenvec_q=f["eigenvector"][:]

    #for freq, e in zip(omega, eigenvec):
    #    condition = freq < max_freq
    #    _freq = np.extract(condition, freq)
    #    _e = np.extract(condition, e)
    #    freqs += list(_freq)
    #    mode_prop += list(_e)
    #gamma1=np.array(mode_prop).ravel()
    #omega1=np.array(freqs).ravel()
    #return(omega1,gamma1,qpoint)
    #omega=np.array(omega)
    #mode_prop=np.array(eigenvec)
    #print np.array(omega_q).shape
    print np.array(eigenvec_q).shape
    #print np.array(qpoint_q).shape
    #print omega_q[0,0:32]

    return(omega_q,eigenvec_q,qpoint_q)

def check(data1,data2):
    data1_1d = np.array(data1).ravel()
    data2_1d = np.array(data2).ravel()
    for d1, d2 in zip(data1_1d,data2_1d):
        if abs(d1 - d2) > 0.1:
           print "Too large difference!"
           print d1, d2
    print "check finished"

def project_eigenvec(edata):
    datasize=edata.shape
    ampdata = np.zeros((datasize[0],datasize[2]))
    for i in range(0,datasize[0]):
        for j in range(0,datasize[2]):
            amp = 0
            for k in range(0,datasize[1]):
                if k % 3 == 0 or k % 3 == 1:
                    amp += abs(edata[i,k,j])**2
            ampdata[i,j]=amp
    return(ampdata)

def select_mode(omega,gamma,sqamp):
    for f, g, s in zip(omega,gamma,sqamp):
        condition = f < max_freq
        _f = np.extract(condition,f)
        _g = np.extract(condition,g)
        _s = np.extract(condition,s)
        f1 = np.array(_f).ravel()
        g1 = np.array(_g).ravel()
        s1 = np.array(_s).ravel()
    return(f1,g1,s1)

def run():
    omegak,gammak,qpointk=parse_gamma(s,Temp)
    omegaq,eigenvecq,qpointq=parse_eigenvec(sv)
    check(omegak,omegaq)
    check(qpointk,qpointq)
    sqamp=project_eigenvec(eigenvecq) 
    #print sqamp[0,:]
    #print sqamp.shape
    print omegak.shape
    print gammak.shape
    print sqamp.shape
    omega1d,gamma1d,sqamp1d=select_mode(omegak,gammak,sqamp)
    plt.figure(figsize=(12,12))
    plt.scatter(sqamp1d,gamma1d)
run()
plt.show()
