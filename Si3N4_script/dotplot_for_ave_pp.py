#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import h5py


sf = "/home/kazu/bsi3n4_m/phono3py_113_fc2_224_sym_monk_shift/noiso_ave_pp/kappa-m8820.ave_pp.hdf5"
gf = "/home/kazu/bge3n4_m/phono3py_113_fc2_224_sym/noiso_ave_pp/kappa-m8820.ave_pp.hdf5"
cf = "/home/kazu/bc3n4_m/phono3py_113_fc2_224_sym/noiso_ave_pp/kappa-m8820.ave_pp.hdf5"
limomega = 20


def parse_data(filename):
    f=h5py.File(filename)
    ap=f['ave_pp'][:,:]
    om=f['frequency'][:,:]
    aps=ap.shape
    ap1=ap.reshape( (aps[0]*aps[1]) )
    om1=om.reshape( (aps[0]*aps[1]) )
    dap=[]
    dom=[]
    j=0
    for i in om1:
        if i <= limomega:
            dap.append(ap1[j])
            dom.append(i)
        j=j+1
    return(dap,dom)

def plotter(x,y,i,clab,label):
    plt.subplot(1,1,i)
    plt.title("ave_pp vs omega")
    plt.plot(x,y,'.',color=clab,label=label,markersize=1,fillstyle='full') 
    plt.legend(loc='upper left')
    plt.xlim(0,limomega)
    #plt.ylim(0,1.5)
    #plt.yticks([0,0.5,1.0,1.5])
    plt.xlabel("Omega. [THz]")
    plt.ylabel("ave_pp [xx]")



def run():
    sap,som=parse_data(sf)
    gap,gom=parse_data(gf)
    cap,com=parse_data(cf)
    plt.figure(figsize=(6,6))
    #plt.plot(som,sap,'.',markersize=3,fillstyle='full') 
    plotter(som,sap,1,'blue','Si3N4')
    plotter(gom,gap,1,'green','Ge3N4')
    plotter(com,cap,1,'red','C3N4')
    #plt.show()
    plt.savefig("avepp-omega.eps")

run()

