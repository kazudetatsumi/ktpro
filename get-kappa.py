#!/user/bin/env python
import numpy as np
import h5py,sys


subjs = []
dirs = []
meshes =[]
subjs.append("bc3n4")
dirs.append("/home/kazu/bc3n4_m/phono3py_113_fc2_338_sym/jdos_t300/tmp/")
meshes.append("101026")
subjs.append("bsi3n4")
dirs.append("/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym/jdos_t300/tmp/")
meshes.append("101026")
subjs.append("bge3n4")
dirs.append("/home/kazu/bge3n4_m/phono3py_113_fc2_338_sym/jdos_t300/tmp/")
meshes.append("101026")
subjs.append("ac3n4")
dirs.append("/home/kazu/ac3n4/phono3py_112_fc2_334_sym_monk_shift/jdos_t300/tmp/")
meshes.append("101014")
subjs.append("asi3n4")
dirs.append("/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos_t300/tmp/")
meshes.append("101014")
subjs.append("age3n4")
dirs.append("/home/kazu/age3n4/phono3py_112_fc2_334_sym_monk_shift/jdos_t300/tmp/")
meshes.append("101014")

temp = 300

def generate_filenames(dirs, meshes):
    kappafiles = []
    irfiles = []

    for d, m in zip(dirs, meshes):
        kappafiles.append(d + "kappa-m" + m + ".hdf5")
        irfiles.append(d + "ir_grid_points.yaml")
    return(kappafiles,irfiles)
 

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
    f = open(irf,'r')
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
            rvlc[2,2]= float(values[4][:-1])
            vol = 1/np.dot(rvec[0,:],np.cross(rvec[1,:],rvec[2,:]))

    f.close()
    g=np.array(g)
    return(g,w,qs,vol,n)

def jdos_wj(g,omega,dirname,mesh):
    jdata=np.loadtxt(dirname + "jdos-m"+mesh+"-g" + str(g) + "-t300.dat",dtype='float')
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

def get_ps(gp,wt,qp,omega,qs,dirname,mesh):
    tmp_jdata=np.loadtxt(dirname + "jdos-m"+mesh+"-g" + str(gp[0]) + "-t300.dat",dtype='float')
    totjdata=np.zeros_like(tmp_jdata)
    totjdata[:,0]=tmp_jdata[:,0]
    i=0
    ps=0
    for g, w in zip(gp,wt):
        if np.linalg.norm(qp[i,:] - qs[i,:]) > 0.01:
           print "qp and qs are much different! \n"
           sys.exit()

        ps += jdos_wj(g,omega[i,:],dirname,mesh)*w
        i += 1
    return(ps)

def get_kappa(hdf):
    f = h5py.File(hdf)
    temperature = f["temperature"].value
    i=0
    for t in temperature:
        if t == temp:
            tindex=i
        i += 1
    kappa=f["kappa"][tindex,] 
    return(kappa)

def run():
    kappafiles,irfiles = generate_filenames(dirs,meshes)
    # get kappa
    for f,s in zip(kappafiles,subjs):
       print s
       kappa=get_kappa(f)
       print kappa[0],kappa[2]
    # get phase space
    for ds,ms,kappafile, irfile in zip(dirs,meshes,kappafiles,irfiles):
       gp,wt,qp,vol,n = get_irdata(irfile)
       omega,qs,sumw = get_freqs(kappafile,n)
       print get_ps(gp,wt,qp,omega,qs,ds,ms) / float(sumw) / vol**2

run()
    
