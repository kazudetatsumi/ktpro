#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from matplotlib import rc
#plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'
homedir = "/home/kazu/"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
Temp = 300
nbins = 300
#y_max = 0.12
numr = 7
max_freq = 33
fs = 9
dosc  = np.loadtxt(cdir + 'total_dos_m292935.dat',comments='#',dtype='float')
doss  = np.loadtxt(sdir + 'total_dos_m292967.dat',comments='#',dtype='float')
dosg  = np.loadtxt(gdir + 'total_dos.dat',comments='#',dtype='float')
gc = cdir + 'noiso/kaccum_m101014.dat'
gs = sdir + 'gpjob_m101026_v1129/kaccum_m101026.dat'
gg = gdir + 'noiso/kaccum.dat'
ggc = cdir +  'noiso/gvaccum_m101014.dat'
ggs = sdir + 'gpjob_m101026_v1129/gvaccum_m101026.dat'
ggg = gdir + 'noiso/gvaccum.dat'
c = cdir + "noiso/kappa-m101014.noiso.hdf5"
s = sdir + "gpjob_m101026_v1129/kappa-m101026.hdf5"
g = gdir + "noiso/kappa-m121212.hdf5"
cv=298.78
sv=143.78
gv=472.14


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
        freq_s = sorted(freq)
        #condition = freq_s[0] == freq
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

def eachplot2(sn,phase,omega,dkaccum):
   plt.subplot(numr,3,sn)
   plt.title("gv_for_" + phase)
   plt.plot(omega,dkaccum[:,0],label=phase + "_dkxx")
   plt.plot(omega,dkaccum[:,2],label=phase + "_dkzz")
   plt.ylim(0,10.5)
   plt.yticks([0,5,10])
   plt.xlim(0,max_freq)

def eachplot3(sn,phase,omega,dos):
   plt.subplot(numr,3,sn)
   plt.title("dos_for_" + phase)
   plt.plot(omega,dos,label=phase + "_dos")
   plt.ylim(0,0.10)
   plt.yticks([0,0.05,0.10])
   plt.xlim(0,max_freq)





def eachplot12(sn,phase,omega,gamma,xmin,xmax,ymin,ymax,title):
       plt.subplot(numr,3,sn)
       plt.title( title + "_for_" +  phase)
       plt.scatter(omega,gamma,c="None", s=0.1,label=phase +"_" + title)
       omega_a,gamma_a=sortomega(omega,gamma)
       print omega_a.shape,gamma_a.shape
       plt.yscale("log")
       plt.ylim(ymin,ymax)
       plt.xlim(xmin,xmax)


def sortomega(x,y):
    data = np.c_[x,y]
    datas = np.array(sorted(data, key=lambda x:x[0]))
    b = np.ones(fs)/float(fs)
    a = datas[:,1]
    data_ave = np.convolve(a,b,'same')
    ds = datas.shape
    print ds
    print data_ave.shape
    cs = (fs-1)/2
    return datas[:,0],data_ave 

def run():
   omegac,kaccumc,dkaccumc=parse_kaccum(gc)
   omegas,kaccums,dkaccums=parse_kaccum(gs)
   omegag,kaccumg,dkaccumg=parse_kaccum(gg)
   omegagc,gvaccumc,dgvaccumc=parse_gvaccum(ggc)
   omegags,gvaccums,dgvaccums=parse_gvaccum(ggs)
   omegagg,gvaccumg,dgvaccumg=parse_gvaccum(ggg)
   omegac1,gammac1=parse_gamma(c,Temp,0)
   omegas1,gammas1=parse_gamma(s,Temp,0)
   omegag1,gammag1=parse_gamma(g,Temp,0)


   plt.figure(figsize=(11,15.375))
   plt.rcParams['pdf.fonttype'] = 42

   eachplot3(1,"alpha",dosc[:,0],dosc[:,1]/cv)
   eachplot3(2,"beta",doss[:,0],doss[:,1]/sv)
   eachplot3(3,"gamma",dosg[:,0],dosg[:,1]/gv*4)
   eachplot(4,"alpha",omegac,kaccumc,dkaccumc)
   eachplot(5,"beta",omegas,kaccums,dkaccums)
   eachplot(6,"gamma",omegag,kaccumg,dkaccumg)
   eachplot2(7,"alpha",omegagc,dgvaccumc)
   eachplot2(8,"beta",omegags,dgvaccums)
   eachplot2(9,"gamma",omegagg,dgvaccumg)
   eachplot12(10,"alpha",omegac1,gammac1,0,max_freq,0.0005,0.35,"gamma")
   eachplot12(11,"beta",omegas1,gammas1,0,max_freq,0.0005,0.35,"gamma")
   eachplot12(12,"gamma",omegag1,gammag1,0,max_freq,0.0005,0.35,"gamma")
   #plt.savefig("tst_plot.eps")



run()
plt.show()
