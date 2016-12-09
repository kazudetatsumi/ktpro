#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
import csv


f = h5py.File("kappa-m8820.hdf5")
kxx = f['kappa'][:,0]
kzz = f['kappa'][:,2]
t = f['temperature'][:]

g = h5py.File("/home/kazu/waln/phono3py_332_fc2_443_sym_monk/kappa-m212111.hdf5")
kxx_aln = g['kappa'][:,0]
kzz_aln = g['kappa'][:,2]

h = h5py.File("/home/kazu/asi3n4/phono3py_112_fc2_222_sym_monk_shift/kappa-m8810.hdf5")
kxx_a = h['kappa'][:,0]
kzz_a = h['kappa'][:,2]

i = h5py.File("/home/kazu/mgsin2/phono3py_222_sym_monk/kappa-m11911.hdf5")
kxx_mg = i['kappa'][:,0]
kyy_mg = i['kappa'][:,1]
kzz_mg = i['kappa'][:,2]

j = h5py.File("/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/kappa-m121212.hdf5")
kxx_g = j['kappa'][:,0]


#11111csv obj = csv.render(open("/home/kazu/asi3n4/aSi3N4_kappa_expt_cvd.csv",'r'))
aexpdata = np.loadtxt("/home/kazu/asi3n4/aSi3N4_kappa_expt_cvd.csv",delimiter=",",skiprows=1)
bexpdatah = np.loadtxt("/home/kazu/bsi3n4_m/bSi3N4_kappa_expt_hirosaki.csv",delimiter=",",skiprows=1)


texp=[300,600,900]
kexp=[107,76,45]
kexpaln=[256,96,55]
kexpmg=[23,15,11]
texp_single=[300,300]
kexp_single=[180,69]

taniso=[300]
k_beta_xx=[71]
k_beta_zz=[193]
k_alpha_xx=[65.9]
k_alpha_zz=[95.6]
k_aln_xx=[217]
k_aln_zz=[190]
k_mgsin2_xx=[35.4]
k_mgsin2_yy=[40.3]
k_mgsin2_zz=[38.9]

plt.figure(figsize=(8,6))
plt.title("kappa vs temp")
plt.plot(t,kxx,"b",label="kxx_beta")
plt.plot(t,kzz,"b",label="kzz_beta")
#plt.plot(t,(kxx*2+kzz)/3,"b",label="k_beta")
##plt.plot(t,(kzz+kxx*2)/3,label="kam_beta")
##plt.plot(t,3*(1/kzz+1/kxx*2)**(-1),label="khm_beta")
plt.plot(t,kxx_aln,"k",label="kxx_aln")
plt.plot(t,kzz_aln,"k",label="kzz_aln")
plt.plot(t,(kxx_aln*2+kzz_aln)/3,"k",label="k_aln")
##plt.plot(t,(kxx_aln*2+kzz_aln)/3,label="kam_aln")
##plt.plot(t,3*(1/kzz_aln+1/kxx_aln*2)**(-1),label="khm_aln")
##plt.plot(t,(kzz_aln+kxx_aln*2)/3,label="kave_aln")
plt.plot(t,kxx_a,"r",label="kxx_alpha")
plt.plot(t,kzz_a,"r",label="kzz_alpha")
#plt.plot(t,(kxx_a*2+kzz_a)/3,"r",label="k_alpha")
##plt.plot(t,(kzz_a+kxx_a*2)/3,label="kam_alpha")
##plt.plot(t,3*(1/kzz_a+1/kxx_a*2)**(-1),label="kam_alpha")

plt.plot(t,kxx_g,"violet",label="k_gamma")
plt.plot(t,kxx_mg,"g",label="kxx_mgsin2")
plt.plot(t,kyy_mg,"g",label="kyy_mgsin2")
#plt.plot(t,kzz_mg,"g",label="kzz_mgsin2")
#plt.plot(t,(kxx_mg+kyy_mg+kzz_mg)/3,"g",label="k_mgsin2")
#plt.plot(texp,kexp,"b:s",label="kexp_beta",markersize=6)
#plt.plot(texp,kexpaln,'k--o',label="kexp_aln",markersize=6)
#plt.plot(texp,kexpmg,'g--o',label="kexp_mgsin2",markersize=6)
plt.plot(texp_single,kexp_single,'o',label="kexp_signle",markersize=6)

#plt.plot(aexpdata[:,0],aexpdata[:,1],"r--o",label="kexp_alpha",markersize=6)
#plt.plot(bexpdatah[:,0],bexpdatah[:,1],'b--o',label="kexp_betaH",markersize=6)

#plt.plot(taniso,k_beta_xx,"co",label="k_beta_xx",markersize=6)
#plt.plot(taniso,k_beta_zz,"co",label="k_beta_zz",markersize=6)
#plt.plot(taniso,k_alpha_xx,"ro",label="k_alpha_xx",markersize=6)
#plt.plot(taniso,k_alpha_zz,"ro",label="k_alpha_zz",markersize=6)
#plt.plot(taniso,k_aln_xx,"k_",label="k_aln_xx",markersize=6)
#plt.plot(taniso,k_aln_zz,"k_",label="k_aln_zz",markersize=6)
#plt.plot(taniso,k_mgsin2_xx,"g_",label="k_mgsin2_xx",markersize=6)
#plt.plot(taniso,k_mgsin2_zz,"g_",label="k_mgsin2_zz",markersize=6)
#plt.plot(taniso,k_mgsin2_yy,"g_",label="k_mgsin2_yy",markersize=6)

#plt.legend()
plt.xlim(250,1000)
#plt.ylim(0,300)
plt.ylim(9,300)
plt.xlabel("Temp. [K]")
plt.ylabel("Kappa [W/K.m]")
plt.yscale("log")
#plt.yticks([5,10,20,40,80,160])
#plt.xscale("log")
#plt.show()
plt.savefig("kappa-temp_withaniso.eps")
