#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import stats
from matplotlib import rc
plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'
homedir = "/home/kazu/"
cdir = homedir + "/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
sdir = homedir + "/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/"
gdir = homedir + "/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/"
Temp = 300
nbins = 300
numr = 7
max_freq = 33
fs = 9
gc = cdir + 'noiso/kaccum_m101014.dat'
gs = cdir + 'noiso/kaccum_v1129.dat'
#gs = sdir + 'noiso/kaccum_m101026.dat'
gg = gdir + 'noiso/kaccum.dat'
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




def run():
   omegac,kaccumc,dkaccumc=parse_kaccum(gc)
   omegas,kaccums,dkaccums=parse_kaccum(gs)
   omegag,kaccumg,dkaccumg=parse_kaccum(gg)
   plt.figure(figsize=(11,10.050))
   plt.rcParams['pdf.fonttype'] = 42

   eachplot(4,"alpha",omegac,kaccumc,dkaccumc)
   eachplot(4,"alphav1129",omegas,kaccums*10*10*14,dkaccums*10*10*14)
   eachplot(6,"gamma",omegag,kaccumg,dkaccumg)
   plt.savefig("tst_plot.eps")



run()
#plt.show()
