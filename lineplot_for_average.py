#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from matplotlib import rc

Temp = 300
nbins = 100
y_max = 0.12
numr = 3
homedir = "/home/kazu/"
cdir = homedir + "asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
#c = cdir + "noiso/kappa-m8810.hdf5"
#s = sdir + "noiso/kappa-m8820.hdf5"
#g = gdir + "noiso/kappa-m121212.hdf5"
#apc= cdir + "gpjob_m8810_fullpp/kappa-m8810.hdf5"
apc= cdir + "gpjob_m101014_fullpp/kappa-m101014.hdf5"
#aps= sdir + "gpjob_m8820_fullpp/kappa-m8820.hdf5"
aps= sdir + "gpjob_m101026_fullpp/kappa-m101026.hdf5"
apg= gdir + "gpjob_m121212_fullpp/kappa-m121212.hdf5"
#cjc= cdir + "kappa-m8810.const_ave1.hdf5"
#cjs= sdir + "kappa-m8820.const_ave1.hdf5"
#cjg= gdir + "kappa-m121212.const_ave1.hdf5"
cz=4
sz=2
gz=2




 

def parse_gamma(filename,temp,max_freq):
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

def parse_avepp(filename,max_freq):
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


   

def run():
   x = np.array([])
   cc = np.array([])
   cs = np.array([])
   cg = np.array([])
   agc = np.array([])
   ags = np.array([])
   agg = np.array([])
   aac = np.array([])
   aas = np.array([])
   aag = np.array([])
   aagc = np.array([])
   aags = np.array([])
   aagg = np.array([])
   maxf = np.arange(2.0,35,0.1)
   for max_freq in maxf:
      #omegac1,gammac1=parse_gamma(c,Temp,max_freq)
      #omegas1,gammas1=parse_gamma(s,Temp,max_freq)
      #omegag1,gammag1=parse_gamma(g,Temp,max_freq)
      omegaapc1,apc1=parse_avepp(apc,max_freq)
      omegaaps1,aps1=parse_avepp(aps,max_freq)
      omegaapg1,apg1=parse_avepp(apg,max_freq)
      #omegacjc1,gammacjc1=parse_gamma(cjc,Temp,max_freq)
      #omegacjs1,gammacjs1=parse_gamma(cjs,Temp,max_freq)
      #omegacjg1,gammacjg1=parse_gamma(cjg,Temp,max_freq)
      #cc = np.append(cc, np.corrcoef(gammac1,apc1*gammacjc1)[0,1])
      #cs = np.append(cs,np.corrcoef(gammas1,aps1*gammacjs1)[0,1])
      #cg = np.append(cg,np.corrcoef(gammag1,apg1*gammacjg1)[0,1])
      #agc = np.append(agc,np.average(gammac1))
      #ags = np.append(ags,np.average(gammas1))
      #agg = np.append(agg,np.average(gammag1))
      aac = np.append(aac,np.average(apc1))
      aas = np.append(aas,np.average(aps1))
      aag = np.append(aag,np.average(apg1))
      #aagc = np.append(aagc,np.average(apc1*gammacjc1))
      #aags = np.append(aags,np.average(aps1*gammacjs1))
      #aagg = np.append(aagg,np.average(apg1*gammacjg1))
      x = np.append(x,max_freq)
   #gas=np.c_[gammas1,aps1]
   #print gas.shape
   plt.figure(figsize=(16,16))
   #plt.subplot(3,1,1)
   #plt.plot(x,cc,label="cc_alpha")
   #plt.plot(x,cs,label="cc_beta")
   #plt.plot(x,cg,label="cc_gamma")
   plt.legend(loc="lower right")
   #plt.subplot(3,1,2)
   plt.plot(x,aac*cz**2,label="ap_alpha")
   plt.plot(x,aas*sz**2,label="ap_beta")
   plt.plot(x,aag*gz**2,label="ap_gamma")
   print x[130],aac[130]*cz**2,aas[130]*sz**2,aag[130]*gz**2
   print x[329],aac[329]*cz**2,aas[329]*sz**2,aag[329]*gz**2
   print x[100],aac[100]*cz**2,aas[100]*sz**2,aag[100]*gz**2
   print x[280],aac[280]*cz**2,aas[280]*sz**2,aag[280]*gz**2
   plt.legend(loc="upper right")
   #eachplot12(1,"alpha",x,cc,0.001,0.03,0.00000000005,0.0000000010,"avepp")
#   eachplot12(2,"beta",gammas1,aps1,0.001,np.max(gammas1),0.00000000005,np.max(aps1),"avepp")
#   eachplot12(3,"gamma",gammag1,apg1,0.001,np.max(gammag1),0.00000000005,np.max(apg1),"avepp")
#   eachplot12(4,"alpha",gammac1,gammacjc1/4,0.001,0.03,3.0*10**7,2.5*10**8,"wjdos")
#   eachplot12(5,"beta",gammas1,gammacjs1,0.001,0.03,3.0*10**7,2.5*10**8,"wjdos")
#   eachplot12(6,"gamma",gammag1,gammacjg1,0.001,0.03,1.6*10**7,0.75*10**8,"wjdos")
#   eachplot12(7,"alpha",gammac1,apc1*gammacjc1,0.001,0.03,0.001,0.25,"wjdos*avepp")
#   eachplot12(8,"beta",gammas1,aps1*gammacjs1,0.001,0.03,0.001,0.25,"wjdos*avepp")
#   eachplot12(9,"gamma",gammag1,apg1*gammacjg1,0.001,0.03,0.001,0.25,"wjdos*avepp")
#   plt.tight_layout()
   #plt.savefig("tst_plot.pdf")

#maxf = np.arange(2.0,35,1.0)
#for max_freq in maxf:
#    print max_freq
run()
plt.show()
