#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import re
import h5py
from math import exp


h = 6.626070040 * 1e-34
kb = 1.3806488 * 1e-23
Temp = 300
max_freq = 35
width = 9
homedir = "/home/kazu/"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
fc = cdir +  'noiso/gvaccum.dat'
gc = cdir + 'noiso/kaccum.dat'
fs = sdir + 'noiso/gvaccum.dat'
gs = sdir + 'noiso/kaccum.dat'
fg = gdir + 'noiso/gvaccum.dat'
gg = gdir + 'noiso/kaccum.dat'
c = cdir + "noiso/kappa-m8810.hdf5"
s = sdir + "noiso/kappa-m8820.hdf5"
g = gdir + "noiso/kappa-m121212.hdf5"


vc = 298.78277783767
vs = 148.59906968339
vg = 118.0350288

my_tau = 20.0

def parse_gaccum(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line.strip() == "":
                continue
            data.append([float(x) for x in line.split()])
    return np.array(data)

def parse_kaccum(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                values = line.split()
                if int(values[1]) == Temp:
                    break

        for line in f:
            if line[0] == '#':
                break
            if line.strip() == "":
                continue
            data.append([float(x) for x in line.split()])
    return np.array(data)

def modelkaccum2(file_gaccum, file_kaccum, file_kappa):
    g_data = parse_gaccum(file_gaccum)
    k_data = parse_kaccum(file_kaccum)

    omegak = g_data[:,0]
    gvaccum = g_data[:,1:4]
    dgvaccum = g_data[:,7:10]
    kaccum = k_data[:,1:4]

    efac = h * omegak * 10**12 / kb / Temp
    cv2 = kb * efac**2 * np.exp(efac) / (np.exp(efac) - 1)**2
    #Here gvaccum is multiplied with a constant mode heat capacity and a constant lifetime

    modelkaccum2 = np.zeros_like(dgvaccum)
    #max_omega_index = 98
    max_omega_index = dgvaccum.shape[0]-2
    print dgvaccum.shape
    for j in range(1, max_omega_index + 2):
        d_omegak = omegak[j] - omegak[j-1]
        delta = cv2[j] * dgvaccum[j,:] * d_omegak * 1e34 * 1e-12 * my_tau
        modelkaccum2[j,:] = modelkaccum2[j-1,:] + delta

    return (omegak, gvaccum, kaccum,  modelkaccum2)

def modelkaccum3(file_gaccum, file_kappa):
    gv_data = parse_gaccum(file_gaccum)
    gm_data = parse_gamma(file_kappa,Temp,1)
    gm_datasma = sma(gm_data)

    omegak = gv_data[:,0]
    gvaccum = gv_data[:,1:4]
    dgvaccum = gv_data[:,7:10]


    efac = h * omegak * 10**12 / kb / Temp
    cv2 = kb * efac**2 * np.exp(efac) / (np.exp(efac) - 1)**2

    modelkaccum3 = np.zeros_like(dgvaccum)
    max_omega_index = dgvaccum.shape[0]-2
    print dgvaccum.shape
    for j in range(1, max_omega_index + 2): 
        gamma = getNearestOmegaGamma(gm_datasma,omegak[j-1])
        d_omegak = omegak[j] - omegak[j-1]
        delta = cv2[j] * dgvaccum[j,:] * d_omegak * 1e34 * 1e-12 / (2*gamma*2*np.pi)
        modelkaccum3[j,:] = modelkaccum3[j-1,:] + delta

    return (omegak, modelkaccum3)


def getNearestOmegaGamma(data, num):
    idx = np.abs(np.asarray(data[:,0]) - num).argmin()
    return data[idx,1] 

def parse_gamma(filename,temp,frag):
    freqs = []
    mode_prop = []
    f = h5py.File(filename,'r')
    temperature = f["temperature"].value
    i=0
    for t in temperature:
        if t == temp:
            tindex=i
        i += 1
    gamma=f["gamma"][tindex,]
    omega=f["frequency"][:,:]

    for freq, g in zip(omega, gamma):
        condition = freq < max_freq
        _freq = np.extract(condition, freq)
        _g = np.extract(condition, g)
        freqs += list(_freq)
        mode_prop += list(_g)
    gamma1=np.array(mode_prop).ravel()
    omega1=np.array(freqs).ravel()
    if frag == 1:
        omeganz1,gammanz1=remove_zero(omega1,gamma1)
    else:
        omeganz1,gammanz1=omega1,gamma1
    data = np.c_[omeganz1,gammanz1]
    data_sorted = np.array(sorted(data, key=lambda omeganz1:omeganz1[0]))
    #return(omeganz1,gammanz1)
    return(data_sorted)


def remove_zero(omega,gamma):
    gs = []
    freqs = []
    for freq, g in zip(omega, gamma):
        condition = g > 0
        _g = np.extract(condition,g)
        _freq = np.extract(condition,freq)
        gs += list(_g)
        freqs += list(_freq)
    gsnz = np.array(gs)
    fsnz = np.array(freqs)
    return(fsnz,gsnz)


def sma(data):
    b = np.ones(width)/float(width)
    a = data[:,1]
    data_ave = np.c_[data[:,0],np.convolve(a,b,mode='same')]
    return data_ave

def getgamma(x,y,om):
    data = np.c_[x,y]
    datas = np.array(sorted(data, key=lambda x:x[0]))

def run():
    print "alpha"
    omegac,gvaccumc,kaccumc,modelkaccum2c=modelkaccum2(fc,gc,c)
    print omegac.shape,kaccumc.shape
    print "beta"
    omegas,gvaccums,kaccums,modelkaccum2s=modelkaccum2(fs,gs,s)
    print "gamma"
    omegag,gvaccumg,kaccumg,modelkaccum2g=modelkaccum2(fg,gg,g)
    #omegac1,gammac1=parse_gamma(c,Temp,1)
    #print gammac1
    #omegas1,gammas1=parse_gamma(s,Temp,1)
    #omegag1,gammag1=parse_gamma(g,Temp,1)
    #omegac1s,gammac1s=sortomega(omegac1,gammac1)
    omega3c,modelkaccum3c=modelkaccum3(fc,g)
    omega3s,modelkaccum3s=modelkaccum3(fs,g)
    omega3g,modelkaccum3g=modelkaccum3(fg,g)
    print modelkaccum3c[modelkaccum3c.shape[0]-1,0]
    print modelkaccum3c[modelkaccum3c.shape[0]-1,2]
    print modelkaccum3s[modelkaccum3s.shape[0]-1,0]
    print modelkaccum3s[modelkaccum3s.shape[0]-1,2]
    print modelkaccum3g[modelkaccum3g.shape[0]-1,0]



    plt.figure(figsize=(6,15))
    plt.subplot(3,1,1)
    plt.plot(omegac,kaccumc[:,0],"b",label="alpha_kaccumxx")
    plt.plot(omegac,kaccumc[:,2],"b",label="alpha_kaccumzz")
    plt.plot(omegas,kaccums[:,0],"g",label="beta_kaccumxx")
    plt.plot(omegas,kaccums[:,2],"g",label="beta_kaccumzz")
    plt.plot(omegag,kaccumg[:,0],"r",label="gamma_kaccumxx")
    plt.xlim(0,40)
    plt.ylim(0,300)
    plt.title("kaccum")
    plt.subplot(3,1,2)
    plt.plot(omegac,modelkaccum2c[:,0],"b",label="alpha_mkaccumxx")
    plt.plot(omegac,modelkaccum2c[:,2],"b",label="alpha_mkaccumzz")
    plt.plot(omegas,modelkaccum2s[:,0],"g",label="beta_mkaccumxx")
    plt.plot(omegas,modelkaccum2s[:,2],"g",label="beta_mkaccumzz")
    plt.plot(omegag,modelkaccum2g[:,0],"r",label="gamma_mkaccumxx")
    plt.xlim(0,40)
    plt.ylim(0,300)
    plt.title("gv*gv*clammda*tau, tau is the fitting param")
    plt.subplot(3,1,3)
    plt.plot(omega3c,modelkaccum3c[:,0],"b",label="alpha_mkaccumxx")
    plt.plot(omega3c,modelkaccum3c[:,2],"b",label="alpha_mkaccumzz")
    plt.plot(omega3s,modelkaccum3s[:,0],"g",label="beta_mkaccumxx")
    plt.plot(omega3s,modelkaccum3s[:,2],"g",label="beta_mkaccumzz")
    plt.plot(omega3g,modelkaccum3g[:,0],"r",label="gamma_mkaccumxx")
    plt.xlim(0,40)
    plt.ylim(0,300)
    plt.title("gv*gv*clammda*tau, tau is from kaccum file")
    #plt.plot(omegac1s,gammac1s)

run()
plt.show()
#plt.savefig("kaccum_modelkaccum5_gvaccum_plot.eps")

# print(parse_gaccum(fb))

