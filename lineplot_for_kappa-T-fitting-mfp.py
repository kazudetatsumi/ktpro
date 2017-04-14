#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import csv
import scipy.optimize 

t = np.arange(250,1610,10)
initialCoeff = 0.5
texp_single=[300,300]
kexp_single=[180,69]
alpha_hdf5 = "/home/kazu/asi3n4/phono3py_112_fc2_222_sym_monk_shift/kappa-m141416.hdf5"
beta_hdf5 = "/home/kazu/bsi3n4_m/phono3py_113_fc2_224_sym_monk_shift/gpjob_m141432_1micron/kappa-m141432.hdf5"
gamma_hdf5 = "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/kappa-m222222.hdf5"
aexp_file = "/home/kazu/asi3n4/aSi3N4_kappa_expt_cvd.csv"
bexp_file = "/home/kazu/bsi3n4_m/bSi3N4_kappa_expt_hirosaki.csv"
fixed_vols_t=[300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600]

def parse_hdf5(filename):
    print filename
    f = h5py.File(filename)
    kdata = f["kappa"]
    tdata = f["temperature"] 
    return(tdata,kdata)


def get_txxdata(temp):
    filename = "/home/kazu/bsi3n4_m/fixed_vols/t"+str(temp)+"/phono3py_113_fc2_224_sym_monk/kappa-m141432.hdf5"
    f = h5py.File(filename)
    temperature = f["temperature"]
    i=0
    for t in temperature:
        if t == temp:
            tindex=i 
        i += 1
    print  tindex,f["temperature"][tindex],f["kappa"][tindex,0]
    return(f["kappa"][tindex,0],f["kappa"][tindex,2])

def fixed_vols_kappa(tlist):
    kdata = np.zeros((len(tlist),2))
    j=0
    for i in tlist:
        kdata[j,0],kdata[j,1]=get_txxdata(i)
        j += 1
    return kdata

def invTfunc(x,a,b):
	return a*x**(-1.0)+b

def invT_fit(tdata,data):
    p, c = scipy.optimize.curve_fit(invTfunc,data[:,0],data[:,1])
    fity = invTfunc(tdata,p[0],p[1])
    return(fity)

def get_nearestpoints(texp,tdata,kdata):
    tmpsize = texp.shape
    tmpsize2 = kdata.shape
    tnear = np.zeros((tmpsize[0]))
    knear = np.zeros((tmpsize[0],tmpsize2[1]))
    i=0
    for t in texp:
        err = np.abs(np.ones_like(tdata)*t-tdata)
        index = np.argmin(err)
        tnear[i]=tdata[index]
        knear[i,:]=kdata[index,:]
        i += 1
    return tnear,knear 

def objFunc(coef,texp,kexp,tdata,kdata):
    tnear,knear=get_nearestpoints(texp,tdata,kdata)
    tmpsize=texp.shape
    base = np.zeros((tmpsize[0],2)) 
    base[:,0]=knear[:,0]
    base[:,1]=knear[:,2]
    r = kexp - theoFunc(coef,base)
    return r

def theoFunc(coef,base):
    f = coef[0]*base[:,0] + (1-coef[0])*base[:,1]
    return f 

def get_sigma(coef,texp,kexp,tdata,kdata):
    rs = objFunc(coef,texp,kexp,tdata,kdata)
    numele = len(rs)
    rsdot = (np.dot(rs,rs))**(0.5)/numele
    return rsdot
def get_relerr(coef,texp,kexp,tdata,kdata):
    rs = objFunc(coef,texp,kexp,tdata,kdata)
    rele = rs/kexp
    return rele


def fit_proc(coef,texp,kexp,tdata,kdata):
    coeffID = scipy.optimize.leastsq(objFunc, coef, args=(texp,kexp,tdata,kdata))
    print "coef:",coeffID
    base = np.c_[kdata[:,0],kdata[:,2]]
    k_fit = theoFunc(coeffID[0],base)
    sigma = get_sigma(coeffID[0],texp,kexp,tdata,kdata)
    print "sigma:",sigma
    rele = get_relerr(coeffID[0],texp,kexp,tdata,kdata)
    print "rele:",rele
    return k_fit

def plotset():
   plt.xlim(298,1600)
   plt.ylim(0,200)
   plt.xlabel("Temp. [K]")
   plt.ylabel("Kappa [W/K.m]")
   plt.legend(loc='upper right')

def run():    
   t_a,k_a = parse_hdf5(alpha_hdf5)
   t_b,k_b = parse_hdf5(beta_hdf5)
   t_g,k_g = parse_hdf5(gamma_hdf5)

   aexpdata = np.loadtxt(aexp_file,delimiter=",",skiprows=1)
   bexpdata = np.loadtxt(bexp_file,delimiter=",",skiprows=1)
   ykexpa = invT_fit(t,aexpdata)
   ykexpb = invT_fit(t,bexpdata)

   k_a_fit = fit_proc(initialCoeff, aexpdata[:,0],aexpdata[:,1],t_a,k_a)
   k_b_fit = fit_proc(initialCoeff, bexpdata[:,0],bexpdata[:,1],t_b,k_b)

   fixed_vols_k=fixed_vols_kappa(fixed_vols_t) 

   plt.figure(figsize=(10,16))
   plt.subplot(2,1,1)
   plt.title("kappa vs temp for alpha")
   plt.plot(t_a,k_a[:,0],"r",label="kxx_alpha")
   plt.plot(t_a,k_a[:,2],"r",label="kzz_alpha")
   plt.plot(t_a,k_a_fit,"r",label="k_fit",linewidth=2)
   plt.plot(t,ykexpa,"r--",label="kexp_a/T+b")
   plt.plot(aexpdata[:,0],aexpdata[:,1],"ro",label="kexp_alpha",markersize=10)
   #plt.plot(t_g,k_g[:,0],"g",label="kxx_gamma")
   plotset()

   plt.subplot(2,1,2)
   plt.title("kappa vs temp for beta")
   plt.plot(t_b,k_b[:,0],"b",label="kxx_beta")
   plt.plot(t_b,k_b[:,2],"b",label="kzz_beta")
   plt.plot(t_b,k_b_fit,"b",label="k_fit",linewidth=2)
   plt.plot(bexpdata[:,0],bexpdata[:,1],'o',label="kexp_betaH",markersize=10)
   plt.plot(fixed_vols_t,fixed_vols_k[:,0],'x',label="fixed_vols_kxx",markersize=10)
   plt.plot(fixed_vols_t,fixed_vols_k[:,1],'x',label="fixed_vols_kzz",markersize=10)
   plt.plot(t,ykexpb,"b--",label="kexpH_a/T+b")
   plt.plot(texp_single,kexp_single,'bd',label="kexp_signle",markersize=10)
   #plt.plot(t_g,k_g[:,0],"g",label="kxx_gamma")
   plotset()


   plt.show()
  # plt.savefig("kappa-temp-fitting-1600K.eps")

run()
