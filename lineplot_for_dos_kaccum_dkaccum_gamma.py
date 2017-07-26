#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from matplotlib import rc

Temp = 300
nbins = 100
y_max = 0.12
numr = 7
jdosc1 = np.loadtxt('/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g0-t300.dat',comments='#',dtype='float')
jdoss1 = np.loadtxt('/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g0-t300.dat',comments='#',dtype='float')
jdosg1 = np.loadtxt('/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/jdos-m222222-g0-t300.dat',comments='#',dtype='float')
jdosc2 = np.loadtxt('/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g1673-t300.dat',comments='#',dtype='float')
jdoss2 = np.loadtxt('/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g3241-t300.dat',comments='#',dtype='float')
jdosg2 = np.loadtxt('/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/jdos-m222222-g5577-t300.dat',comments='#',dtype='float')
dosc  = np.loadtxt('/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/total_dos_m292935.dat',comments='#',dtype='float')
doss  = np.loadtxt('/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/total_dos_m292967.dat',comments='#',dtype='float')
dosg  = np.loadtxt('/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/total_dos.dat',comments='#',dtype='float')
gc = '/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/noiso/kaccum.dat'
gs = '/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/noiso/kaccum.dat'
gg = '/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/noiso/kaccum.dat'
ggc = '/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/noiso/gvaccum.dat'
ggs = '/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/noiso/gvaccum.dat'
ggg = '/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/noiso/gvaccum.dat'
c = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/kappa-m8810.hdf5"
s = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/noiso/kappa-m8820.hdf5"
g = "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/noiso/kappa-m121212.hdf5"
grc = "/home/kazu/asi3n4/gruneisen/gruneisen.hdf5"
grs = "/home/kazu/bsi3n4_m/gruneisen/gruneisen.hdf5"
grg = "/home/kazu/gamma-si3n4-unit/gruneisen/gruneisen.hdf5"
nuc = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/kappa-m8810.hdf5"
nus = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/noiso/gpjob_m8820_nu/kappa-m8820.hdf5"
nug = "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/noiso/gpjob_m121212_nu/kappa-m121212.hdf5"
apc= "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift//gpjob_m8810_fullpp/kappa-m8810.hdf5"
aps= "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/gpjob_m8820_fullpp/kappa-m8820.hdf5"
apg= "/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/gpjob_m121212_fullpp/kappa-m121212.hdf5"


def parse_kaccum(filename):
    data = []
    count = 0 
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                values = line.split()
                if float(values[1]) == Temp:
                    break
                count += 1

        for line in f:
            if line[0] == '#':
                break
            if line.strip() == "": 
                continue
            data.append([float(x) for x in line.split()])
    kaccum = np.array(data)[:,1:4]
    dkaccum = np.array(data)[:,7:10]
    omega = np.array(data)[:,0]
    return (omega,kaccum,dkaccum)

def run_KDE(x, y, nbins, y_max):
    x_min = 0
    x_max = np.rint(x.max() * 1.1)
    y_min = 0
    _y_max = y_max
    values = np.vstack([x.ravel(), y.ravel()])
    kernel = stats.gaussian_kde(values)

    xi, yi = np.mgrid[x_min:x_max:nbins*1j, y_min:_y_max:nbins*1j]
    positions = np.vstack([xi.ravel(), yi.ravel()])
    zi = np.reshape(kernel(positions).T, xi.shape)

    return xi, yi, zi



def parse_gvaccum(filename):
    data = []
    count = 0 
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                break
            if line.strip() == "": 
                continue
            data.append([float(x) for x in line.split()])
    gvaccum = np.array(data)[:,1:4]
    dgvaccum = np.array(data)[:,7:10]
    omega = np.array(data)[:,0]
    return (omega,gvaccum,dgvaccum)    
 
def parse_gamma(filename,temp):
    f = h5py.File(filename,'r')
    temperature = f["temperature"].value
    i=0 
    for t in temperature:
        if t == temp:
            tindex=i 
        i += 1
    gamma=f["gamma"][tindex,]
    omega=f["frequency"][:,:]
    weights=f["weight"][:]
    gammashape=gamma.shape
    gamma1=gamma.reshape( (gammashape[0]*gammashape[1]) )
    omega1=omega.reshape( (gammashape[0]*gammashape[1]) )
    dgamma=[]
    domega=[]
    j=0 
    for i in omega1:
        if i <= 12: 
            dgamma.append(gamma1[j])
            domega.append(i)
        j=j+1
    print np.average(dgamma),np.std(dgamma)

    freqs = []
    mode_prop = []
    mode_weights = []
    for w, freq, g in zip(weights, omega, gamma):
        tau = 1.0 / np.where(g > 0, g, -1) / (2 * 2 * np.pi) * 0.001
        
        condition = tau > 0
        _tau = np.extract(condition, tau)
        _freq = np.extract(condition, freq)

        freqs += list(_freq) * w
        mode_prop += list(_tau) * w
    x = np.array(freqs)
    y = np.array(mode_prop)
    return(omega1,gamma1,x,y) 

def parse_gammanu(filename,temp):
    f = h5py.File(filename,'r')
    print filename
    temperature = f["temperature"].value
    i=0 
    for t in temperature:
        if t == temp:
            tindex=i 
        i += 1
    gammau=f["gamma_U"][tindex,]
    gamman=f["gamma_N"][tindex,]
    omega=f["frequency"][:,:]
    weights=f["weight"][:]
    gammashape=gammau.shape
    gammau1=gammau.reshape( (gammashape[0]*gammashape[1]) )
    gamman1=gamman.reshape( (gammashape[0]*gammashape[1]) )
    omega1=omega.reshape( (gammashape[0]*gammashape[1]) )
    return(omega1,gammau1,gamman1) 

def parse_avepp(filename):
    f = h5py.File(filename,'r')
    avepp=f["ave_pp"][:,:]
    omega=f["frequency"][:,:]
    aveppshape=avepp.shape
    print avepp[10,3]
    avepp1=avepp.reshape( (aveppshape[0]*aveppshape[1]) )
    omega1=omega.reshape( (aveppshape[0]*aveppshape[1]) )
    return(omega1,avepp1) 

def parse_gruneisen(filename):
    f = h5py.File(filename,'r')
    gruneisen=f["gruneisen"][:,:]
    omega=f["frequency"][:,:]
    weights=f["weight"][:]

    freqs = []
    mode_prop = []
    mode_weights = []
    for w, freq, g in zip(weights, omega, gruneisen):
        freqs += list(freq) * w 
        mode_prop += list(g) * w 
    x = np.array(freqs)
    y = np.array(mode_prop)

    shapex=x.shape
    print shapex
    shapey=y.shape
    print shapey

    return(x,y)


def eachplot(sn,phase,omega,kaccum,dkaccum):
   plt.subplot(numr,3,sn)
   plt.title("kaccum_for_" + phase)
   plt.plot(omega,kaccum[:,0],label=phase + "_kxx")
   plt.plot(omega,kaccum[:,2],label=phase + "_kzz")
   plt.plot(omega,dkaccum[:,0]*10,label=phase + "_dkxx")
   plt.plot(omega,dkaccum[:,2]*10,label=phase + "_dkzz")
   plt.ylim(0,255)
   plt.yticks([0,100,200])
   plt.xlim(0,15)
   #plt.xscale("log")
   #plt.xlabel("mfp [micro-meter]")
   #plt.ylabel("kaccum [W/m.K]")
   #plt.legend(loc='upper left')

def eachplot2(sn,phase,omega,dkaccum):
   plt.subplot(numr,3,sn)
   plt.title("gv_for_" + phase)
   plt.plot(omega,dkaccum[:,0]*10,label=phase + "_dkxx")
   plt.plot(omega,dkaccum[:,2]*10,label=phase + "_dkzz")
   plt.ylim(0,105)
   plt.yticks([0,50,100])
   plt.xlim(0,15)

def eachplot3(sn,phase,omega,dos):
   plt.subplot(numr,3,sn)
   plt.title("dos_for_" + phase)
   plt.plot(omega,dos,label=phase + "_dos")
   plt.ylim(0,3)
   plt.yticks([0,1,2,3])
   plt.xlim(0,15)

def eachplot4(sn,phase,omega,gamma):
   plt.subplot(numr,3,sn)
   plt.title("lifetime_for_" + phase)
   plt.scatter(omega,1/(2*gamma*2*np.pi),s=0.1,label=phase +"_tau")
   #plt.scatter(omega,gamma,s=0.1,label=phase +"_scattering_rate")
   plt.yscale("log")
   plt.ylim(2,100)
   #plt.yticks([0,0.04,0.08,0.12])
   plt.xlim(0,15)

def eachplot5(xi, yi, zi, nbins, sn,phase,omega,gamma):
   plt.subplot(numr,3,sn)
   plt.title("tau_for_" + phase)
   #plt.pcolormesh(xi[:,:nbins], yi[:,:nbins], zi[:,:nbins],cmap='OrRd')
   #plt.pcolormesh(xi[:,:nbins], yi[:,:nbins], zi[:,:nbins])
   #plt.colorbar()
   plt.scatter(omega,1/(2*gamma*2*np.pi)*0.001,s=0.1,label=phase +"_tau")
   plt.ylim(0,0.12)
   #plt.yticks([0,50,100])
   plt.xlim(0,15)

def eachplot6(sn,phase,omega1,dos1,omega2,dos2):
   plt.subplot(numr,3,sn)
   plt.title("jdos_for_" + phase)
   plt.plot(omega1,dos1,label=phase + "_jdos1")
   plt.plot(omega2,dos2,label=phase + "_jdos2")
   plt.ylim(0,10)
   #plt.yticks([0,1,2,3])
   plt.xlim(0,15)

def eachplot7(sn,phase,omega,gruneisen):
   plt.subplot(numr,3,sn)
   plt.title("gruneisen_for_" + phase)
   plt.scatter(omega,gruneisen,s=0.1,label=phase +"_gruneisen")
   #plt.yscale("log")
   plt.ylim(-2.5,2.5)
   #plt.yticks([0,0.04,0.08,0.12])
   plt.xlim(0,15)

def eachplot8(sn,phase,omega,gammau,gamman):
   plt.subplot(numr,3,sn)
   plt.title("unscattering_rate_for_" + phase)
   plt.scatter(omega,gammau,s=0.1,color="red",label=phase +"_uscattering_rate")
   plt.scatter(omega,gamman,s=0.1,color="blue",label=phase +"_nscattering_rate")
   plt.ylim(0.00,0.03)
   plt.xlim(0,15)
   plt.legend(loc="upper left")

def eachplot9(sn,phase,omega,avepp):
   plt.subplot(numr,3,sn)
   plt.title("avepp_for_" + phase)
   plt.scatter(omega,avepp,s=0.1,label=phase +"_avepp_rate")
   plt.yscale("log")
   plt.ylim(0.0000000020,0.00000000004)
   #plt.yticks([0,0.04,0.08,0.12])
   plt.xlim(0,15)

def run():
   omegac,kaccumc,dkaccumc=parse_kaccum(gc)
   omegas,kaccums,dkaccums=parse_kaccum(gs)
   omegag,kaccumg,dkaccumg=parse_kaccum(gg)
   omegagc,gvaccumc,dgvaccumc=parse_gvaccum(ggc)
   omegags,gvaccums,dgvaccums=parse_gvaccum(ggs)
   omegagg,gvaccumg,dgvaccumg=parse_gvaccum(ggg)
   omegac1,gammac1,xc,yc=parse_gamma(c,Temp)
   omegas1,gammas1,xs,ys=parse_gamma(s,Temp)
   omegag1,gammag1,xg,yg=parse_gamma(g,Temp)
   omeganus1,gammasu1,gammasn1=parse_gammanu(nus,Temp)
   omeganug1,gammagu1,gammagn1=parse_gammanu(nug,Temp)
   omegaaps1,aps1=parse_avepp(aps)
   omegaapg1,apg1=parse_avepp(apg)
   ##xci, yci, zci = run_KDE(xc, yc, nbins, y_max)
   ##xsi, ysi, zsi = run_KDE(xs, ys, nbins, y_max)
   ##xgi, ygi, zgi = run_KDE(xg, yg, nbins, y_max)
   xc,yc=parse_gruneisen(grc)
   xs,ys=parse_gruneisen(grs)
   xg,yg=parse_gruneisen(grg)



   plt.figure(figsize=(16,20))
#   rc('text', usetex=True)
   rc('font', family='serif')
   rc('font', serif='Times New Roman')
   plt.rcParams['pdf.fonttype'] = 42

   eachplot3(1,"alpha",dosc[:,0],dosc[:,1]/4)
   eachplot3(2,"beta",doss[:,0],doss[:,1]/2)
   eachplot3(3,"gamma",dosg[:,0],dosg[:,1]/2)
   eachplot(4,"alpha",omegac,kaccumc,dkaccumc)
   eachplot(5,"beta",omegas,kaccums,dkaccums)
   eachplot(6,"gamma",omegag,kaccumg,dkaccumg)
   eachplot4(7,"alpha",omegac1,gammac1)
   eachplot4(8,"beta",omegas1,gammas1)
   eachplot4(9,"gamma",omegag1,gammag1)
   #eachplot4(10,"alpha",omegac1,gammac1)
   #eachplot4(11,"betau",omeganus1,gammasu1)
   #eachplot4(12,"gamman",omeganug1,gammagu1)
   #eachplot4(13,"alpha",omegac1,gammac1)
   #eachplot4(14,"betan",omeganus1,gammasn1)
   #eachplot4(15,"gamman",omeganug1,gammagn1)
   #eachplot8(11,"beta",omeganus1,gammasu1,gammasn1)
   #eachplot8(12,"gamma",omeganug1,gammagu1,gammagn1)
   #eachplot5(xgi,ygi,zgi,nbins,9,"gamma",omegag1,gammag1)
   eachplot2(10,"alpha",omegagc,dgvaccumc)
   eachplot2(11,"beta",omegags,dgvaccums)
   eachplot2(12,"gamma",omegagg,dgvaccumg)
   eachplot6(13,"alpha",jdosc1[:,0],(jdosc1[:,1]+jdosc1[:,2])/16,jdosc2[:,0],(jdosc2[:,1]+jdosc2[:,2])/16)
   eachplot6(14,"beta",jdoss1[:,0],(jdoss1[:,1]+jdoss1[:,2])/4,jdoss2[:,0],(jdoss2[:,1]+jdoss2[:,2])/4)
   eachplot6(15,"gamma",jdosg1[:,0],(jdosg1[:,1]+jdosg1[:,2])/4,jdosg2[:,0],(jdosg2[:,1]+jdosg2[:,2])/4)
   eachplot7(16,"alpha",xc,yc)
   eachplot7(17,"beta",xs,ys)
   eachplot7(18,"gamma",xg,yg)
   eachplot9(19,"beta",omegaaps1,aps1)
   eachplot9(20,"beta",omegaaps1,aps1)
   eachplot9(21,"gamma",omegaapg1,apg1)
   plt.tight_layout()
   #plt.savefig("tst_plot.pdf")

run()
plt.show()
