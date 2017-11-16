#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import math
homedir = "/home/kazu/"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
Temp = 300
max_freq = 5
fs = 9
c = cdir + "noiso/kappa-m141416.noiso.hdf5"
s = sdir + "noiso/kappa-m141432.noiso.hdf5"
g = gdir + "noiso/kappa-m121212.hdf5"
cv = cdir + "qpoints.hdf5"
sv = sdir + "qpoints.hdf5"
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
            for k in range(0,datasize[1]):
                if k % 3 == 0 or k % 3 == 1:
                    amp += abs(edata[i,k,j])**2
            ampdata[i,j]=amp
    return(ampdata)

def trans(qpoint,edata):
    datasize=edata.shape
    longidata = np.zeros((datasize[0],datasize[2]))
    for i in range(0,datasize[0]):
        norm = (qpoint[i,0]**2 + qpoint[i,1]**2 + qpoint[i,2]**2)**0.5
        if norm > 0:
           for j in range(0,datasize[2]):
               longi = 0
               for k in range(0,datasize[1]):
                   if k % 3 == 0:
                      longi += abs(edata[i,k,j]*qpoint[i,0]/norm)**2
                   if k % 3 == 1:
                      longi += abs(edata[i,k,j]*qpoint[i,1]/norm)**2
                   if k % 3 == 2:
                      longi += abs(edata[i,k,j]*qpoint[i,2]/norm)**2
               longidata[i,j]=longi
    return(longidata)

    

def select_mode(omega,gamma,sqamp,lamp):
    freqs = []
    gammas = []
    sqamps = []
    lamps = []
    for f, g, s, l in zip(omega,gamma,sqamp,lamp):
        condition = f < max_freq 
        _f = np.extract(condition,f)
        _g = np.extract(condition,g)
        _s = np.extract(condition,s)
        _l = np.extract(condition,l)
        freqs += list(_f)
        gammas += list(_g)
        sqamps += list(_s)
        lamps += list(_l)
    f1 = np.array(freqs).ravel()
    g1 = np.array(gammas).ravel()
    s1 = np.array(sqamps).ravel()
    l1 = np.array(lamps).ravel()
    return(f1,g1,s1,l1)

def caserun(casefile,vcasefile,n,phase):
    omegak,gammak,qpointk=parse_gamma(casefile,Temp)
    omegaq,eigenvecq,qpointq=parse_eigenvec(vcasefile)
    check(omegak,omegaq)
    check(qpointk,qpointq)
    sqamp=project_eigenvec(eigenvecq) 
    longiamp=trans(qpointq,eigenvecq)
    omega1d,gamma1d,sqamp1d,longiamp1d=select_mode(omegak,gammak,sqamp,longiamp)
    plt.subplot(3,2,n)
    plt.scatter(omega1d,gamma1d,c=longiamp1d,linewidth=0.01,s=5, label=phase)
    plt.yscale("log")
    plt.ylim(0.0005,0.0050)
    plt.xlabel('omega / THz')
    plt.ylabel('gamma')
    plt.legend()
    plt.colorbar(label='sum of squares of eigenvector component along q')
    plt.subplot(3,2,n+2)
    plt.scatter(omega1d,gamma1d,c=sqamp1d,linewidth=0.01,s=5, label=phase)
    if phase == "beta":
        x = 0.1*np.arange(0,50)
        y = 0.001*np.exp(0.16*(x-0.5)) 
        plt.plot(x,y)
        n_up = 0
        n_dn = 1
        for omega,gamma in zip(omega1d,gamma1d):
            if gamma >=  0.001*math.exp(0.16*(omega-0.5)):
                    n_up += 1
            else:
                    n_dn += 1
        print "num of up states:",n_up
        print "num of dn states:",n_dn

    plt.yscale("log")
    plt.ylim(0.0005,0.0050)
    plt.xlabel('omega / THz')
    plt.ylabel('gamma')
    plt.legend()
    plt.colorbar(label='sum of squares of eigenvector x_y_component')


def run():
    plt.figure(figsize=(12,12))
    caserun(c,cv,1,"alpha")
    caserun(s,sv,2,"beta")
    caserun(g,gv,3,"gamma")

run()
plt.show()
