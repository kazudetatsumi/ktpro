#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from matplotlib import rc

Temp = 300
nbins = 100
y_max = 0.12
max_freq = 15
numr = 3
c = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/noiso/kappa-m8810.hdf5"
s = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/noiso/kappa-m8820.hdf5"
g = "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/noiso/kappa-m121212.hdf5"
apc= "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/gpjob_m8810_fullpp/kappa-m8810.hdf5"
aps= "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/gpjob_m8820_fullpp/kappa-m8820.hdf5"
apg= "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/gpjob_m121212_fullpp/kappa-m121212.hdf5"
cjc= "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/kappa-m8810.const_ave1.hdf5"
cjs= "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/kappa-m8820.const_ave1.hdf5"
cjg= "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/kappa-m121212.const_ave1.hdf5"



 

def parse_gamma(filename,temp):
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
    return(omega1,gamma1) 

def parse_avepp(filename):
    freqs = []
    mode_prop = []
    f = h5py.File(filename,'r')
    avepp=f["ave_pp"][:,:]
    omega=f["frequency"][:,:]
    for freq, ap in zip(omega, avepp):
        condition = freq < max_freq
        _freq = np.extract(condition, freq)
        _ap = np.extract(condition, ap)
        freqs += list(_freq)
        mode_prop += list(_ap)
    avepp1=np.array(mode_prop).ravel()
    omega1=np.array(freqs).ravel()
    return(omega1,avepp1) 

def remove_zero(gamma, mode_prop):
    gs = []
    ps = []
    for g, p in zip(gamma, mode_prop):
        condition = g > 0
        _g = np.extract(condition,g)
        _p = np.extract(condition,p)
        gs += list(_g)
        ps += list(_p)
    gsnz = np.array(gs)
    psnz = np.array(ps)
    return(gsnz,psnz)

def eachplot11(sn,phase,gamma,avepp):
   plt.subplot(numr,3,sn)
   plt.title("gamma-avepp_for_" + phase)
   plt.scatter(gamma,avepp,s=0.1)

def eachplot12(sn,phase,gamma,mode_prop,xmin,xmax,ymin,ymax,title):
   plt.subplot(numr,3,sn)
   plt.title('gamma-' + title + '_for_' + phase)
   gnz, pnz = remove_zero(gamma,mode_prop) 
   pnz2, gnz2 = remove_zero(pnz,gnz) 
   print np.corrcoef(gnz2,pnz2)[0,1],phase,title
   #plt.scatter(gamma,mode_prop,s=3)
   plt.scatter(gnz2,pnz2,s=0.1)
   plt.xscale("log")
   plt.yscale("log")
   plt.xlabel("gamma")
   plt.ylabel(title)
   string = "corr. coef = %3.2f " % (np.corrcoef(gnz2,pnz2)[0,1])
   plt.text(np.max(gnz2)*0.1,np.min(pnz2)*1.1,string)
   #print min(gnz),phase,title
   #print min(pnz),phase,title
   plt.xlim(np.min(gnz2)*0.9,np.max(gnz2)*1.1)
   plt.ylim(np.min(pnz2)*0.9,np.max(pnz2)*1.1)
   #plt.xlim(0.00000000000,0.0500000000)


def run():
   omegac1,gammac1=parse_gamma(c,Temp)
   omegas1,gammas1=parse_gamma(s,Temp)
   omegag1,gammag1=parse_gamma(g,Temp)
   omegaapc1,apc1=parse_avepp(apc)
   omegaaps1,aps1=parse_avepp(aps)
   omegaapg1,apg1=parse_avepp(apg)
   omegacjc1,gammacjc1=parse_gamma(cjc,Temp)
   omegacjs1,gammacjs1=parse_gamma(cjs,Temp)
   omegacjg1,gammacjg1=parse_gamma(cjg,Temp)

   #gas=np.c_[gammas1,aps1]
   #print gas.shape



   plt.figure(figsize=(16,16))
   eachplot12(1,"alpha",gammac1,apc1,0.001,0.03,0.00000000005,0.0000000010,"avepp")
   eachplot12(2,"beta",gammas1,aps1,0.001,np.max(gammas1),0.00000000005,np.max(aps1),"avepp")
   eachplot12(3,"gamma",gammag1,apg1,0.001,np.max(gammag1),0.00000000005,np.max(apg1),"avepp")
   eachplot12(4,"alpha",gammac1,gammacjc1/4,0.001,0.03,3.0*10**7,2.5*10**8,"wjdos")
   eachplot12(5,"beta",gammas1,gammacjs1,0.001,0.03,3.0*10**7,2.5*10**8,"wjdos")
   eachplot12(6,"gamma",gammag1,gammacjg1,0.001,0.03,1.6*10**7,0.75*10**8,"wjdos")
   eachplot12(7,"alpha",gammac1,apc1*gammacjc1,0.001,0.03,0.001,0.25,"wjdos*avepp")
   eachplot12(8,"beta",gammas1,aps1*gammacjs1,0.001,0.03,0.001,0.25,"wjdos*avepp")
   eachplot12(9,"gamma",gammag1,apg1*gammacjg1,0.001,0.03,0.001,0.25,"wjdos*avepp")
   plt.tight_layout()
   #plt.savefig("tst_plot.pdf")

run()
plt.show()
