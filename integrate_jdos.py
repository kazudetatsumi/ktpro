#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt 
import h5py,sys

dirname = "/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/jdos_t300/tmp/"
irfile = dirname + "ir_grid_points.yaml"
kappafile = dirname + "kappa-m8820.hdf5"

def get_freqs(kf,n):
    f = h5py.File(kf)
    freqs = f['frequency']
    qpts = f['qpoint']
    sumw = sum(f['weight'])
    if n != sumw:
       print "n and sumw are different!"
       sys.exit()
    return(freqs,qpts,sumw)
    

def get_irdata(irf):
    g=[]
    w=[]
    qs=np.empty((0,3),float)
    rvec=np.zeros((3,3))
    f = open(irfile,'r')
    for l in f:
        if "mesh:" in l:
            values = l.split()
            n= int(values[2][:-1])*int(values[3][:-1])*int(values[4])
        if "grid_point:" in l:
            values = l.split()
            g.append(int(values[2]))
        if "weight:" in l:
            values = l.split()
            w.append(int(values[1]))
        if "q-point:" in l:
            values = l.split()
            qx= float(values[2][:-1])
            qy= float(values[3][:-1])
            qz= float(values[4])
            qs = np.append(qs, np.array([[qx,qy,qz]]), axis=0)
        if "a*" in l:
            values = l.split()
            rvec[0,0]= float(values[2][:-1])
            rvec[0,1]= float(values[3][:-1])
            rvec[0,2]= float(values[4][:-1])
        if "b*" in l:
            values = l.split()
            rvec[1,0]= float(values[2][:-1])
            rvec[1,1]= float(values[3][:-1])
            rvec[1,2]= float(values[4][:-1])
        if "c*" in l:
            values = l.split()
            rvec[2,0]= float(values[2][:-1])
            rvec[2,1]= float(values[3][:-1])
            rvec[2,2]= float(values[4][:-1])
            vol = 1/np.dot(rvec[0,:],np.cross(rvec[1,:],rvec[2,:]))
    
    f.close()
    g=np.array(g)
    w=np.array(w)
    return(g,w,qs,vol,n)

def jdos_wj(g,omega):
    jdata=np.loadtxt(dirname + "jdos-m8820-g" + str(g) + "-t300.dat",dtype='float')
    totx=0
    for omj in omega:
        mindiff = abs(jdata[0,0] - omj)
        mink = 0
        k = 0
        for om in jdata[:,0]:
            if abs(om - omj) < mindiff:
               mindiff = abs(om - omj)
               mink = k
            k += 1
        if mink < jdata.shape[0]-1 and mink > 0:
           if jdata[mink,0] < omj and omj < jdata[mink+1,0]:
              x = (jdata[mink+1,1:3] - jdata[mink,1:3]) / (jdata[mink+1,0] - jdata[mink,0]) * (omj - jdata[mink,0]) + jdata[mink,1:3]  
              totx += x
           if jdata[mink-1,0] < omj and omj < jdata[mink,0]:
              x = (jdata[mink+1,1:3] - jdata[mink,1:3]) / (jdata[mink+1,0] - jdata[mink,0]) * (omj - jdata[mink,0]) + jdata[mink,1:3]  
              totx += x
    return(np.sum(totx))

def get_ps(gp,wt,qp,omega,qs):
    tmp_jdata=np.loadtxt(dirname + "jdos-m8820-g" + str(gp[0]) + "-t300.dat",dtype='float')
    totjdata=np.zeros_like(tmp_jdata)
    totjdata[:,0]=tmp_jdata[:,0]
    i=0
    ps=0
    for g, w in zip(gp,wt):
        if np.linalg.norm(qp[i,:] - qs[i,:]) > 0.01:
           print "qp and qs are much different! \n"
           sys.exit()
         
        ps += jdos_wj(g,omega[i,:])*w
        i += 1

    return(ps)
        
    


def run():
    gp,wt,qp,vol,n=get_irdata(irfile)
    omega,qs,sumw=get_freqs(kappafile,n)
    print "phase space=", get_ps(gp,wt,qp,omega,qs) / float(sumw) / vol**2
   
            
run() 
