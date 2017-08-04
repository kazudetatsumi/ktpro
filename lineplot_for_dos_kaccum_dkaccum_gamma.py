#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from matplotlib import rc
homedir = "/home/kazu/stern"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
Temp = 300
nbins = 300
#y_max = 0.12
numr = 7
max_freq = 15
jdosc1 = np.loadtxt(cdir + 'jdos-m141416-g0-t300.dat',comments='#',dtype='float')
jdoss1 = np.loadtxt(sdir + 'jdos-m141432-g0-t300.dat',comments='#',dtype='float')
jdosg1 = np.loadtxt(gdir + 'jdos-m222222-g0-t300.dat',comments='#',dtype='float') 
jdosc2 = np.loadtxt(cdir + 'jdos-m141416-g1673-t300.dat',comments='#',dtype='float')
jdoss2 = np.loadtxt(sdir + 'jdos-m141432-g3241-t300.dat',comments='#',dtype='float')
jdosg2 = np.loadtxt(gdir + 'jdos-m222222-g5577-t300.dat',comments='#',dtype='float')
dosc  = np.loadtxt(cdir + 'total_dos_m292935.dat',comments='#',dtype='float')
doss  = np.loadtxt(sdir + 'total_dos_m292967.dat',comments='#',dtype='float')
dosg  = np.loadtxt(gdir + 'total_dos.dat',comments='#',dtype='float')
gc = cdir + 'noiso/kaccum.dat'
gs = sdir + 'noiso/kaccum.dat'
gg = gdir + 'noiso/kaccum.dat'
ggc = cdir +  'noiso/gvaccum.dat'
ggs = sdir + 'noiso/gvaccum.dat'
ggg = gdir + 'noiso/gvaccum.dat'
c = cdir + "noiso/kappa-m8810.hdf5"
s = sdir + "noiso/kappa-m8820.hdf5"
g = gdir + "noiso/kappa-m121212.hdf5"
grc = homedir + "/asi3n4/gruneisen/gruneisen.hdf5"
grs = homedir + "/bsi3n4_m/gruneisen/gruneisen.hdf5"
grg = homedir + "/gamma-si3n4-unit/gruneisen/gruneisen.hdf5"
nuc = cdir + "kappa-m8810.hdf5"
nus = sdir + "noiso/gpjob_m8820_nu/kappa-m8820.hdf5"
nug = gdir + "noiso/gpjob_m121212_nu/kappa-m121212.hdf5"
apc= cdir + "gpjob_m8810_fullpp/kappa-m8810.hdf5"
aps= sdir + "gpjob_m8820_fullpp/kappa-m8820.hdf5"
apg= gdir + "gpjob_m121212_fullpp/kappa-m121212.hdf5"
cjc= cdir + "kappa-m8810.const_ave1.hdf5"
cjs= sdir + "kappa-m8820.const_ave1.hdf5"
cjg= gdir + "kappa-m121212.const_ave1.hdf5"


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
            if line[0] > max_freq:
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


def run_KDE2(x, y, nbins, y_max=None):
    x_min = 0
    x_max = np.rint(x.max() * 1.0)
    y_min = 0
    if y_max is None:
        y_max = np.rint(y.max() * 1.0)
    values = np.vstack([x, y])
    kernel = stats.gaussian_kde(values)

    xi, yi = np.mgrid[x_min:x_max:nbins*1j, y_min:y_max:nbins*1j]
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
            if line[0] > max_freq:
                data.append([float(x) for x in line.split()])
    gvaccum = np.array(data)[:,1:4]
    dgvaccum = np.array(data)[:,7:10]
    omega = np.array(data)[:,0]
    return (omega,gvaccum,dgvaccum)    
 

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
    return(omeganz1,gammanz1)

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
    gammau1=gammau.ravel()
    gamman1=gamman.ravel()
    omega1=omega.ravel()
    return(omega1,gammau1,gamman1) 

def parse_avepp(filename,frag):
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
    if frag == 1:
        omeganz1,aveppnz1=remove_zero(omega1,avepp1)
    else:
        omeganz1,aveppnz1=omega1,avepp1
    return(omeganz1,aveppnz1)


def parse_gruneisen(filename):
    f = h5py.File(filename,'r')
    gruneisen=f["gruneisen"][:,:]
    omega=f["frequency"][:,:]
    weights=f["weight"][:]

    freqs = []
    mode_prop = []
    mode_weights = []
    for w, freq, g in zip(weights, omega, gruneisen):
        condition = freq < max_freq
        _freq = np.extract(condition, freq)
        _g = np.extract(condition, g)
        freqs += list(_freq) * w 
        mode_prop += list(_g) * w 
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
   plt.xlim(0,max_freq)
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
   plt.xlim(0,max_freq)

def eachplot3(sn,phase,omega,dos):
   plt.subplot(numr,3,sn)
   plt.title("dos_for_" + phase)
   plt.plot(omega,dos,label=phase + "_dos")
   plt.ylim(0,3)
   plt.yticks([0,1,2,3])
   plt.xlim(0,max_freq)


def eachplot5(xi, yi, zi, nbins, sn,phase,omega,gamma,xmin,xmax,ymin,ymax,title):
   plt.subplot(numr,3,sn)
   plt.title(title + "_" + phase)
   plt.pcolormesh(xi[:,:nbins], yi[:,:nbins], zi[:,:nbins],cmap='OrRd')
   #plt.pcolormesh(xi[:,:nbins], yi[:,:nbins], zi[:,:nbins])
   #plt.colorbar()
   plt.scatter(omega,gamma,s=0.1,label=phase + "_" + title)
   plt.yscale("log")
   plt.ylim(ymin,ymax)
   #plt.yticks([0,50,100])
   plt.xlim(xmin,xmax)

def eachplot6(sn,phase,omega1,dos1,omega2,dos2):
   plt.subplot(numr,3,sn)
   plt.title("wjdos_for_" + phase)
   plt.plot(omega1,dos1,label=phase + "_jdos1")
   plt.plot(omega2,dos2,label=phase + "_jdos2")
   plt.ylim(0,10)
   #plt.yticks([0,1,2,3])
   plt.xlim(0,15)


def eachplot8(sn,phase,omega,gammau,gamman):
   plt.subplot(numr,3,sn)
   plt.title("unscattering_rate_for_" + phase)
   plt.scatter(omega,gammau,s=0.1,color="red",label=phase +"_uscattering_rate")
   plt.scatter(omega,gamman,s=0.1,color="blue",label=phase +"_nscattering_rate")
   plt.ylim(0.00,0.03)
   plt.xlim(0,15)
   plt.legend(loc="upper left")

def eachplot12(sn,phase,omega,gamma,xmin,xmax,ymin,ymax,title):
   plt.subplot(numr,3,sn)
   plt.title( title + "_for_" +  phase)
   plt.scatter(omega,gamma,c="None", s=0.1,label=phase +"_" + title)
   plt.yscale("log")
   plt.ylim(ymin,ymax)
   plt.xlim(xmin,xmax)

def run():
   omegac,kaccumc,dkaccumc=parse_kaccum(gc)
   omegas,kaccums,dkaccums=parse_kaccum(gs)
   omegag,kaccumg,dkaccumg=parse_kaccum(gg)
   omegagc,gvaccumc,dgvaccumc=parse_gvaccum(ggc)
   omegags,gvaccums,dgvaccums=parse_gvaccum(ggs)
   omegagg,gvaccumg,dgvaccumg=parse_gvaccum(ggg)
   omegac1,gammac1=parse_gamma(c,Temp,1)
   omegas1,gammas1=parse_gamma(s,Temp,1)
   omegag1,gammag1=parse_gamma(g,Temp,1)
   omeganus1,gammasu1,gammasn1=parse_gammanu(nus,Temp)
   omeganug1,gammagu1,gammagn1=parse_gammanu(nug,Temp)
   omegaapc1,apc1=parse_avepp(apc,0)
   omegaaps1,aps1=parse_avepp(aps,0)
   omegaapg1,apg1=parse_avepp(apg,0)
   #xci, yci, zci = run_KDE2(omegac1,1/(2*gammac1*2*np.pi), nbins)
   #xsi, ysi, zsi = run_KDE2(omegas1,1/(2*gammas1*2*np.pi), nbins)
   #xgi, ygi, zgi = run_KDE2(omegag1,1/(2*gammag1*2*np.pi), nbins)
   #xapci, yapci, zapci = run_KDE2(omegaapc1,apc1*4*10**(11), nbins, 1.5*10**(2))
   #xapsi, yapsi, zapsi = run_KDE2(omegaaps1,aps1*10**(11), nbins, 1.5*10**(2))
   #xapgi, yapgi, zapgi = run_KDE2(omegaapg1,apg1*10**(11), nbins, 1.5*10**(2))
   #xc,yc=parse_gruneisen(grc)
   #xs,ys=parse_gruneisen(grs)
   #xg,yg=parse_gruneisen(grg)
   omegacjc1,gammacjc1=parse_gamma(cjc,Temp,0)
   omegacjs1,gammacjs1=parse_gamma(cjs,Temp,0)
   omegacjg1,gammacjg1=parse_gamma(cjg,Temp,0)
   #xcjci, ycjci, zcjci = run_KDE2(omegacjc1,gammacjc1/4, nbins, 1.6*10**(8))
   #xcjsi, ycjsi, zcjsi = run_KDE2(omegacjs1,gammacjs1, nbins, 1.6*10**(8))
   #xcjgi, ycjgi, zcjgi = run_KDE2(omegacjg1,gammacjg1, nbins, 1.6*10**(8))
   #xcjapci, ycjapci, zcjapci = run_KDE2(omegacjc1,gammacjc1*apc1*10**3, nbins, 1.0*10**2)
   #xcjapsi, ycjapsi, zcjapsi = run_KDE2(omegacjs1,gammacjs1*aps1*10**3, nbins, 1.0*10**2)
   #xcjapgi, ycjapgi, zcjapgi = run_KDE2(omegacjg1,gammacjg1*apg1*10**3, nbins, 1.0*10**2)



   plt.figure(figsize=(12,16))
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
   #eachplot4(10,"alpha",omegac1,gammac1)
   #eachplot4(11,"betau",omeganus1,gammasu1)
   #eachplot4(12,"gamman",omeganug1,gammagu1)
   #eachplot4(13,"alpha",omegac1,gammac1)
   #eachplot4(14,"betan",omeganus1,gammasn1)
   #eachplot4(15,"gamman",omeganug1,gammagn1)
   #eachplot8(11,"beta",omeganus1,gammasu1,gammasn1)
   #eachplot8(12,"gamma",omeganug1,gammagu1,gammagn1)
   #eachplot5(xci,yci,zci,nbins,7,"alpha",omegac1,1/(2*gammac1*2*np.pi),0,15,2,100,"tau")
   #eachplot5(xsi,ysi,zsi,nbins,8,"beta",omegas1,1/(2*gammas1*2*np.pi),0,15,2,100,"tau")
   #eachplot5(xgi,ygi,zgi,nbins,9,"gamma",omegag1,1/(2*gammag1*2*np.pi),0,15,2,100,"tau")
   eachplot2(7,"alpha",omegagc,dgvaccumc)
   eachplot2(8,"beta",omegags,dgvaccums)
   eachplot2(9,"gamma",omegagg,dgvaccumg)
   #eachplot12(10,"alpha",omegac1,1/(2*gammac1*2*np.pi),0,max_freq,2,100,"lifetime")
   #eachplot12(11,"beta",omegas1,1/(2*gammas1*2*np.pi),0,max_freq,2,100,"lifetime")
   #eachplot12(12,"gamma",omegag1,1/(2*gammag1*2*np.pi),0,max_freq,2,100,"lifetime")
   eachplot12(10,"alpha",omegac1,gammac1,0,max_freq,0.0005,0.08,"gamma")
   eachplot12(11,"beta",omegas1,gammas1,0,max_freq,0.0005,0.08,"gamma")
   eachplot12(12,"gamma",omegag1,gammag1,0,max_freq,0.0005,0.08,"gamma")
   #eachplot6(13,"alpha",jdosc1[:,0],(jdosc1[:,1]+jdosc1[:,2])/16,jdosc2[:,0],(jdosc2[:,1]+jdosc2[:,2])/16)
   #eachplot6(14,"beta",jdoss1[:,0],(jdoss1[:,1]+jdoss1[:,2])/4,jdoss2[:,0],(jdoss2[:,1]+jdoss2[:,2])/4)
   #eachplot6(15,"gamma",jdosg1[:,0],(jdosg1[:,1]+jdosg1[:,2])/4,jdosg2[:,0],(jdosg2[:,1]+jdosg2[:,2])/4)
#   eachplot7(16,"alpha",xc,yc)
#   eachplot12(16,"alpha",xc,yc,0,15,-2.5,2.5,"gruneisen")
#   eachplot12(17,"beta",xs,ys,0,15,-2.5,2.5,"gruneisen")
#   eachplot12(18,"gamma",xg,yg,0,15,-2.5,2.5,"gruneisen")
   #eachplot9(13,"alpha*4",omegaapc1,apc1*4)
   #eachplot5(xapci,yapci,zapci,nbins,13,"alpha*4*10**(11)",omegaapc1,apc1*4*10**11,0,15,4*10**(0),1.5*10**(2),"avepp")
   #eachplot5(xapsi,yapsi,zapsi,nbins,14,"beta*10**(11)",omegaaps1,aps1*10**11,0,15,4*10**(0),1.5*10**(2),"avepp")
   #eachplot5(xapgi,yapgi,zapgi,nbins,15,"gamma*10**(11)",omegaapg1,apg1*10**11,0,15,4*10**(0),1.5*10**(2),"avepp")
   eachplot12(16,"alpha*4",omegaapc1,apc1*4,0,max_freq,4*10**(-11),1.5*10**(-9),"avepp")
   eachplot12(17,"beta",omegaaps1,aps1,0,max_freq,4*10**(-11),1.5*10**(-9),"avepp")
   eachplot12(18,"gamma",omegaapg1,apg1,0,max_freq,4*10**(-11),1.5*10**(-9),"avepp")
   #eachplot5(xcjci,ycjci,zcjci,nbins,16,"alpha/4",omegacjc1,gammacjc1/4,0,15,10**(7),1.6*10**(8),"wjdos")
   #eachplot5(xcjsi,ycjsi,zcjsi,nbins,17,"beta",omegacjs1,gammacjs1,0,15,10**(7),1.6*10**(8),"wjdos")
   #eachplot5(xcjgi,ycjgi,zcjgi,nbins,18,"gamma",omegacjg1,gammacjg1,0,15,10**(7),1.6*10**(8),"wjdos")
   eachplot12(13,"alpha/4",omegacjc1,gammacjc1/4,0,max_freq,10**7,1.6*10**8,"wjdos")
   eachplot12(14,"beta",omegacjs1,gammacjs1,0,max_freq,10**7,1.6*10**8,"wjdos")
   eachplot12(15,"gamma",omegacjg1,gammacjg1,0,max_freq,10**7,1.6*10**8,"wjdos")
   #eachplot5(xcjapci,ycjapci,zcjapci,nbins,19,"alpha*10**3",omegacjc1,gammacjc1*apc1*10**3,0,15,5*10**(-1),1.0*10**(2),"wjdos*avepp")
   #eachplot5(xcjapsi,ycjapsi,zcjapsi,nbins,20,"beta*10**3",omegacjs1,gammacjs1*aps1*10**3,0,15,5*10**(-1),1.0*10**(2),"wjdos*avepp")
   #eachplot5(xcjapgi,ycjapgi,zcjapgi,nbins,21,"gamma*10**3",omegacjg1,gammacjg1*apg1*10**3,0,15,5*10**(-1),1.0*10**(2),"wjdos*avepp")
   eachplot12(19,"alpha",omegacjc1,gammacjc1*apc1,0,max_freq,0.0005,0.1,"wjdos*avepp")
   eachplot12(20,"beta",omegacjs1,gammacjs1*aps1,0,max_freq,0.0005,0.1,"wjdos*avepp")
   eachplot12(21,"gamma",omegacjg1,gammacjg1*apg1,0,max_freq,0.0005,0.1,"wjdos*avepp")
   plt.tight_layout()
   plt.savefig("tst_plot.eps")

run()
#plt.show()
