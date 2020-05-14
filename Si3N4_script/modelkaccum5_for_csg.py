#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import re
from math import exp


h = 6.626070040 * 1e-34
kb = 1.3806488 * 1e-23
Temp = 300
fs = '/home/kazu/bsi3n4_m/phono3py_113_fc2_224_sym_monk_shift/noiso_for_csige/gvaccum.dat'
gs = '/home/kazu/bsi3n4_m/phono3py_113_fc2_224_sym_monk_shift/noiso_for_csige/kaccum.dat'
fc = '/home/kazu/bc3n4_m/phono3py_113_fc2_224_sym/noiso/gvaccum.dat'
gc = '/home/kazu/bc3n4_m/phono3py_113_fc2_224_sym/noiso/kaccum.dat'
fg = '/home/kazu/bge3n4_m/phono3py_113_fc2_224_sym/noiso/gvaccum.dat'
gg = '/home/kazu/bge3n4_m/phono3py_113_fc2_224_sym/noiso/kaccum.dat'

my_tau = 11.5

def parse_gaccum(filename):
    data = []
    with open(filename) as f:
        for line in f:
            if line.strip() == "":
                continue
            data.append([float(x) for x in line.split()])
    return np.array(data)

def parse_kaccum(filename, index=30):
    data = []
    count = 0
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                if count == index:
                    break
                count += 1

        for line in f:
            if line[0] == '#':
                break
            if line.strip() == "":
                continue
            data.append([float(x) for x in line.split()])
    return np.array(data)

def modelkaccum2(file_gaccum, file_kaccum):
    g_data = parse_gaccum(file_gaccum)
    k_data = parse_kaccum(file_kaccum)

    omegak = g_data[:,0]
    gvaccum = g_data[:,1:4]
    dgvaccum = g_data[:,7:10]
    kaccum = k_data[:,1:4]

    efac = h * omegak * 10**12 / kb / Temp
    cv2 = kb * efac**2 * np.exp(efac) / (np.exp(efac) - 1)**2
    cv = kb
    #Here gvaccum is multiplied with a constant mode heat capacity and a constant lifetime
    gvaccum *= cv * 1e34 * 1e-12 * my_tau

    modelkaccum2 = np.zeros_like(dgvaccum)
    max_omega_index = 98
    for j in range(1, max_omega_index + 2):
        d_omegak = omegak[j] - omegak[j-1]
        delta = cv2[j] * dgvaccum[j,:] * d_omegak * 1e34 * 1e-12
        modelkaccum2[j,:] = modelkaccum2[j-1,:] + delta

    limj = max_omega_index
    sumgh1 = sum(kaccum[limj,:] * gvaccum[limj,:])
    sumhh1 = sum(gvaccum[limj,:] * gvaccum[limj,:])
    a1 = sumgh1 / sumhh1
    print("a1:", a1)

    da1 = np.sqrt(sum((kaccum[limj,:] / gvaccum[limj,:] - a1)**2) / 3)
    print("da1:", da1 / a1)

    res1 = np.sqrt(sum((kaccum[limj,:] - a1 * gvaccum[limj,:])**2)/3)
    print("res1:", res1)

    sumgh2 = sum(kaccum[limj,:] * modelkaccum2[limj,:])
    sumhh2 = sum(modelkaccum2[limj,:]**2)
    a2 = sumgh2 / sumhh2
    print("a2:", a2)

    da2 = np.sqrt(sum((kaccum[limj,:] / modelkaccum2[limj,:] - a2)**2) / 3)
    print("da2:", da2 / a2)

    res2 = np.sqrt(sum((kaccum[limj,:] - a2 * modelkaccum2[limj,:])**2) / 3)
    print("res2:", res2)
    
    return (omegak, gvaccum, kaccum, a2 * modelkaccum2)


def run():
    print "C3N4"
    omegac,gvaccumc,kaccumc,modelkaccum2c=modelkaccum2(fc,gc)
    print "Si3N4"
    omegas,gvaccums,kaccums,modelkaccum2s=modelkaccum2(fs,gs)
    print "Ge3N4"
    omegag,gvaccumg,kaccumg,modelkaccum2g=modelkaccum2(fg,gg)
    plt.figure(figsize=(6,15))
    plt.subplot(3,1,1)
    plt.plot(omegas,kaccums[:,0],"b",label="Si_kaccumxx")
    plt.plot(omegas,kaccums[:,2],"b",label="Si_kaccumzz")
    plt.plot(omegag,kaccumg[:,0],"g",label="Ge_kaccumxx")
    plt.plot(omegag,kaccumg[:,2],"g",label="Ge_kaccumzz")
    plt.plot(omegac,kaccumc[:,0],"r",label="C_kaccumxx")
    plt.plot(omegac,kaccumc[:,2],"r",label="C_kaccumzz")
    plt.xlim(0,40)
    plt.ylim(0,300)
    plt.title("kaccum")
    plt.subplot(3,1,2)
    plt.plot(omegas,gvaccums[:,0],"b",label="Si_gvaccumxx")
    plt.plot(omegas,gvaccums[:,2],"b",label="Si_gvaccumzz")
    plt.plot(omegag,gvaccumg[:,0],"g",label="Ge_gvaccumxx")
    plt.plot(omegag,gvaccumg[:,2],"g",label="Ge_gvaccumzz")
    plt.plot(omegac,gvaccumc[:,0],"r",label="C_gvaccumxx")
    plt.plot(omegac,gvaccumc[:,2],"r",label="C_gvaccumzz")
    plt.xlim(0,40)
    plt.ylim(0,300)
    plt.title("gvaccum*c*tau, c and tau are constants")
    plt.subplot(3,1,3)
    plt.plot(omegas,modelkaccum2s[:,0],"b",label="Si_mkaccumxx")
    plt.plot(omegas,modelkaccum2s[:,2],"b",label="Si_mkaccumzz")
    plt.plot(omegag,modelkaccum2g[:,0],"g",label="Ge_mkaccumxx")
    plt.plot(omegag,modelkaccum2g[:,2],"g",label="Ge_mkaccumzz")
    plt.plot(omegac,modelkaccum2c[:,0],"r",label="C_mkaccumxx")
    plt.plot(omegac,modelkaccum2c[:,2],"r",label="C_mkaccumzz")
    plt.xlim(0,40)
    plt.ylim(0,300)
    plt.title("gv*gv*clammda*tau, tau is the fitting param")

run()
#plt.show()
plt.savefig("kaccum_modelkaccum5_gvaccum_plot.eps")

# print(parse_gaccum(fb))

