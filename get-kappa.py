#!/user/bin/env python
import numpy as np
import h5py,sys
import matplotlib.pyplot as plt


subjs = []
dirs = []
meshes =[]
masss =[]
zs =[]
#subjs.append("bc3n4")
#dirs.append("/home/kazu/bc3n4_m/phono3py_113_fc2_338_sym/")
#meshes.append("101026")
#masss.append(184.109)
#zs.append(2)
subjs.append("bsi3n4")
dirs.append("/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym/")
meshes.append("101026")
masss.append(280.5554)
zs.append(2)
subjs.append("bge3n4")
dirs.append("/home/kazu/bge3n4_m/phono3py_113_fc2_338_sym/")
meshes.append("101026")
masss.append(547.8314)
zs.append(2)
#subjs.append("ac3n4")
#dirs.append("/home/kazu/ac3n4/phono3py_112_fc2_334_sym_monk_shift/")
#meshes.append("101014")
#masss.append(368.218)
subjs.append("asi3n4")
dirs.append("/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/")
meshes.append("101014")
masss.append(561.1109)
zs.append(4)
#subjs.append("age3n4")
#dirs.append("/home/kazu/age3n4/phono3py_112_fc2_334_sym_monk_shift/")
#meshes.append("101014")
#masss.append(1095.663)
#zs.append(4)
#subjs.append("be2sio4")
#dirs.append("/home/kazu/be2sio4/phono3py_111_fc2_222/")
#meshes.append("121212")
#masss.append(660.624)
#zs.append(6)
#subjs.append("gsi3n4")
#dirs.append("/home/kazu//gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/")
#meshes.append("121212")
#masss.append(280.5554)
#zs.append(2)
#subjs.append("wAlN")
#dirs.append("/home/kazu/waln/phono3py_332_fc2_443_sym_monk/")
#meshes.append("212111")

temp = 300

zs=np.array(zs)
print zs
def generate_filenames(dirs, meshes):
    kappafiles = []
    kappa_constavp_files = []
    irfiles = []
    kappa_avp_files = []

    for d, m in zip(dirs, meshes):
        kappafiles.append(d + "kappa-m" + m + ".hdf5")
        kappa_constavp_files.append(d + "kappa-m" + m + ".const_ave-pp_noiso.hdf5")
        kappa_avp_files.append(d + "/gpjob_m" + m + "_fullpp/kappa-m" + m + ".ave_pp.hdf5")
        irfiles.append(d + "/jdos_t300/tmp/ir_grid_points.yaml")
    return(kappafiles,kappa_constavp_files,irfiles,kappa_avp_files)
 

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
            vol = np.dot(rvec[0,:],np.cross(rvec[1,:],rvec[2,:]))

    f.close()
    g=np.array(g)
    return(g,w,qs,vol,n)

def jdos_wj(g,omega,dirname,mesh):
    jdata=np.loadtxt(dirname + "jdos_t300/tmp/jdos-m"+mesh+"-g" + str(g) + "-t300.dat",dtype='float')
    #jdata=np.loadtxt(dirname + "jdos/tmp/jdos-m"+mesh+"-g" + str(g) + ".dat",dtype='float')
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
              x = (jdata[mink,1:3] - jdata[mink-1,1:3]) / (jdata[mink,0] - jdata[mink-1,0]) * (omj - jdata[mink-1,0]) + jdata[mink-1,1:3]
              totx += x
    return(np.sum(totx))

def get_ps(gp,wt,qp,omega,qs,dirname,mesh):
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
    kappafiles,kappa_constavp_files,irfiles,kappa_avp_files = generate_filenames(dirs,meshes)
    # get kappa
    kappas=np.empty((0,6),float)
    kappacs=np.empty((0,6),float)
    kappaavps=np.empty((0,6),float)
    for f,fc,favp,s in zip(kappafiles,kappa_constavp_files,kappa_avp_files,subjs):
       print s
       kappa=get_kappa(f)
       kappas = np.append(kappas,([kappa]), axis=0)
       kappac=get_kappa(fc)
       kappacs = np.append(kappacs,([kappac]), axis=0)
       kappaavp=get_kappa(favp)
       kappaavps = np.append(kappaavps,([kappaavp]), axis=0)

    print "material:\n",subjs
    print "kappa:\n",kappas
    print "kappa_ave\n",np.average(kappas[:,0:3], axis=1)
    print "kappa_const_avepp:\n",kappacs
    print "kappa_avepp:\n",kappaavps
    # get phase space
    pss=np.empty((0),float)
    rhos=np.empty((0),float)
    for ds,ms,kappafile,irfile,mass in zip(dirs,meshes,kappafiles,irfiles,masss):
       gp,wt,qp,vol,n = get_irdata(irfile)
       omega,qs,sumw = get_freqs(kappafile,n)
       nj = omega.shape[1]
       ps=get_ps(gp,wt,qp,omega,qs,ds,ms) / float(sumw) / nj**3 * 2 / 3
       pss=np.append(pss,np.array([ps]))
       rho=mass*vol
       rhos=np.append(rhos,np.array([rho]))
    print "phase spaece:\n",pss

    plt.figure(figsize=(13,13))

    plt.subplot(3,3,1)
    plt.scatter(np.average(kappas[:,0:3], axis=1),pss,label="kave")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa")
    plt.ylabel("phase space")
    plt.legend()

    plt.subplot(3,3,2)
    plt.scatter(kappas[:,0],pss,label="kxx")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa")
    plt.ylabel("phase space")
    plt.legend()

    plt.subplot(3,3,3)
    plt.scatter(kappas[:,2],pss,label="kzz")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa")
    plt.ylabel("phase space")
    plt.legend()

    plt.subplot(3,3,4)
    plt.scatter(np.average(kappas[:,0:3],axis=1),np.average(kappacs[:,0:3], axis=1),label="kave")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa")
    plt.ylabel("kappa const avepp")
    plt.legend()

    plt.subplot(3,3,5)
    plt.scatter(kappas[:,0],kappacs[:,0],label="kxx")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa")
    plt.ylabel("kappa const avepp")
    plt.legend()

    plt.subplot(3,3,6)
    plt.scatter(kappas[:,2],kappacs[:,2],label="kzz")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa")
    plt.ylabel("kappa const avepp")
    plt.legend()

    plt.subplot(3,3,7)
    plt.scatter(np.average(kappas[:,0:3],axis=1),np.average(kappaavps[:,0:3],axis=1),label="kavep")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.ylabel("kappa avepp")
    plt.legend()

    plt.subplot(3,3,8)
    plt.scatter(kappas[:,0],kappaavps[:,0],label="kxx")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xlabel("kappa")
    plt.ylabel("kappa avepp")
    plt.legend()

    plt.subplot(3,3,9)
    plt.scatter(kappas[:,2],kappaavps[:,2],label="kavep_zz")
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.ylabel("kappa avepp")
    plt.legend()

run()

plt.show()
#plt.savefig("ph-kappa-wjdos-si3n4.eps")
    
