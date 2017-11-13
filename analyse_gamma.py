#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
homedir = "/home/kazu/"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
Temp = 300
max_freq = 15
fs = 9
c = cdir + "noiso/kappa-m141416.noiso.hdf5"
s = sdir + "noiso/kappa-m141432.noiso.hdf5"
cv = cdir + "qpoints.hdf5"
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
    return(f1,g1,s1)

def caserun(casefile,vcasefile,n,phase):
    omegak,gammak,qpointk=parse_gamma(casefile,Temp)
    omegaq,eigenvecq,qpointq=parse_eigenvec(vcasefile)
    check(omegak,omegaq)
    check(qpointk,qpointq)
    sqamp=project_eigenvec(eigenvecq) 
    omega1d,gamma1d,sqamp1d=select_mode(omegak,gammak,sqamp)
    plt.subplot(2,1,n)
    plt.scatter(omega1d,gamma1d,c=sqamp1d,linewidth=0.01,s=2, label=phase)
    plt.yscale("log")
    plt.ylim(0.0005,0.10)
    plt.xlabel('omega / THz')
    plt.ylabel('gamma')
    plt.legend()
    plt.colorbar(label='sum of squares of eigenvector_x_y_components')

def run():
    plt.figure(figsize=(12,12))
    caserun(c,cv,1,"alpha")
    caserun(s,sv,2,"beta")

run()
plt.show()
