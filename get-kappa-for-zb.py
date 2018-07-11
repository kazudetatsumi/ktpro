#!/user/bin/env python
import numpy as np
import h5py,sys
import matplotlib.pyplot as plt


subjs = []
dirs = []
meshes =[]
subjs.append("gap")
dirs.append("/home/kazu/zb/gap/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("gasb")
dirs.append("/home/kazu/zb/gasb/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("ge")
dirs.append("/home/kazu/zb/ge/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("alsb")
dirs.append("/home/kazu/zb/alsb/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("gaas")
dirs.append("/home/kazu/zb/gaas/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("si")
dirs.append("/home/kazu/zb/si/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("cdte")
dirs.append("/home/kazu/zb/cdte/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("inp")
dirs.append("/home/kazu/zb/inp/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("insb")
dirs.append("/home/kazu/zb/insb/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")
subjs.append("znte")
dirs.append("/home/kazu/zb/znte/phono3py_222_fc2_444_sym_monk/")
meshes.append("191919")

temp = 300

def generate_filenames(dirs, meshes):
    kappafiles = []
    kappa_constavp_files = []
    irfiles = []

    for d, m in zip(dirs, meshes):
        kappafiles.append(d + "kappa-m" + m + ".hdf5")
        kappa_constavp_files.append(d + "kappa-m" + m + ".const_ave-pp_noiso.hdf5")
        irfiles.append(d + "/jdos/tmp/ir_grid_points.yaml")
    return(kappafiles,kappa_constavp_files,irfiles)
 

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
            rvec[2,2]= float(values[4][:-1])
            vol = 1/np.dot(rvec[0,:],np.cross(rvec[1,:],rvec[2,:]))

    f.close()
    g=np.array(g)
    return(g,w,qs,vol,n)

def jdos_wj(g,omega,dirname,mesh):
    #jdata=np.loadtxt(dirname + "jdos_t300/tmp/jdos-m"+mesh+"-g" + str(g) + "-t300.dat",dtype='float')
    jdata=np.loadtxt(dirname + "jdos/tmp/jdos-m"+mesh+"-g" + str(g) + ".dat",dtype='float')
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
    #tmp_jdata=np.loadtxt(dirname + "jdos_t300/tmp/jdos-m"+mesh+"-g" + str(gp[0]) + "-t300.dat",dtype='float')
    #totjdata=np.zeros_like(tmp_jdata)
    #totjdata[:,0]=tmp_jdata[:,0]
    i=0
    ps=0
    for g, w in zip(gp,wt):
        j=0
        for qpe,qse in zip(qp[i,:],qs[i,:]):
            if qse < 0:
               qs[i,j] +=1 
            if qpe < 0:
               qp[i,j] +=1 
            j += 1
        if np.linalg.norm(qp[i,:] - qs[i,:]) > 0.01:
           print "qp and qs are much different! \n"
           print np.linalg.norm(qp[i,:] - qs[i,:])
           print qp[i,:]
           print qs[i,:]
           print i
           sys.exit()

        ps += jdos_wj(g,omega[i,:],dirname,mesh)*w
        i += 1
    return(ps)

def get_kappa(hdf):
    print hdf
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
    kappafiles,kappa_constavp_files,irfiles = generate_filenames(dirs,meshes)
    # get kappa
    kappas=np.empty((0,6),float)
    kappacs=np.empty((0,6),float)
    for f,fc,s in zip(kappafiles,kappa_constavp_files,subjs):
       print s
       kappa=get_kappa(f)
       kappas = np.append(kappas,([kappa]), axis=0)
       kappac=get_kappa(fc)
       kappacs = np.append(kappacs,([kappac]), axis=0)

    print "kappa:\n",kappas
    print "kappa_const_avepp:\n",kappacs
    # get phase space
    pss=np.empty((0),float)
    for ds,ms,kappafile, irfile in zip(dirs,meshes,kappafiles,irfiles):
       gp,wt,qp,vol,n = get_irdata(irfile)
       omega,qs,sumw = get_freqs(kappafile,n)
       ps=get_ps(gp,wt,qp,omega,qs,ds,ms) / float(sumw) / vol**2
       pss=np.append(pss,np.array([ps]))
    print "phase spaece:\n",pss

    plt.figure(figsize=(10,20))
    plt.subplot(2,1,1)
    #plt.scatter(pss,kappas[:,0],label="kxx")
    plt.scatter(kappas[:,2],pss,label="kzz")
    #plt.scatter(pss,np.average(kappas[:,0:3], axis=1),label="kave")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa")
    plt.ylabel("phase space")
    plt.subplot(2,1,2)
    #plt.scatter(kappacs[:,0],kappas[:,0],label="kxx")
    plt.scatter(kappacs[:,2],kappas[:,2],label="kzz")
    #plt.scatter(np.average(kappacs[:,0:3],axis=1),np.average(kappas[:,0:3], axis=1),label="kave")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa const avepp")
    plt.ylabel("kappa")
    plt.legend()
run()

plt.show()
    
