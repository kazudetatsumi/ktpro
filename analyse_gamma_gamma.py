#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import math
homedir = "/home/kazu/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
Temp = 300
max_freq = 30
fs = 9
g = gdir + "noiso/kappa-m121212.hdf5"
gv = gdir + "qpoints.hdf5"


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
    #print np.array(gamma_k).shape
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
    #print np.array(eigenvec_q).shape
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
            tmpamp = []
            tmpamp.append(get_amp(edata[i,:,j],1,1,1))
            tmpamp.append(get_amp(edata[i,:,j],-1,1,1))
            tmpamp.append(get_amp(edata[i,:,j],1,-1,1))
            tmpamp.append(get_amp(edata[i,:,j],1,1,-1))
            ampdata[i,j]=max(tmpamp)
            print max(tmpamp)
    return(ampdata)

def get_amp(data,l,m,n):
    datasize = data.shape
    amp = 0
    for k in range(0,datasize[0]):
        if k % 3 == 0:
           amp += (abs(data[k]*l + data[k+1]*m + data[k+2]*n))**2 / (l**2 + m**2 + n**2)
    return(1 - amp)

    

def select_mode(omega,gamma,sqamp):
    freqs = []
    gammas = []
    sqamps = []
    for f, g, s in zip(omega,gamma,sqamp):
        condition = f < max_freq 
        _f = np.extract(condition,f)
        _g = np.extract(condition,g)
        _s = np.extract(condition,s)
        freqs += list(_f)
        gammas += list(_g)
        sqamps += list(_s)
    f1 = np.array(freqs).ravel()
    g1 = np.array(gammas).ravel()
    s1 = np.array(sqamps).ravel()
    f1 = np.append(f1,0)
    g1 = np.append(g1,0)
    s1 = np.append(s1,0)
    return(f1,g1,s1)

def caserun(casefile,vcasefile,n,phase):
    omegak,gammak,qpointk=parse_gamma(casefile,Temp)
    omegaq,eigenvecq,qpointq=parse_eigenvec(vcasefile)
    check(omegak,omegaq)
    check(qpointk,qpointq)
    sqamp=project_eigenvec(eigenvecq) 
    omega1d,gamma1d,sqamp1d=select_mode(omegak,gammak,sqamp)
    plt.scatter(omega1d,gamma1d,c=sqamp1d,linewidth=0.01,s=15, label=phase)
    plt.yscale("log")
    plt.ylim(0.0005,0.350)
    plt.xlabel('omega / THz')
    plt.ylabel('gamma')
    plt.legend()
    plt.colorbar(label='sum of squares of eigenvector comp on {111}')


def run():
    plt.figure(figsize=(12,12))
    #caserun(c,cv,1,"alpha")
    #caserun(s,sv,2,"beta")
    caserun(g,gv,1,"gamma")

run()
plt.show()
