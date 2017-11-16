#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt 
import math

max_length = 20
tolerance = 0.000000001
nbins = 50
sigma = 0.10
asigma = 0.01

ccar = "/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/POSCAR"
scar = "/home/kazu//bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/POSCAR"


def parse_poscar(filename):
    with open(filename) as f:
        latvec = []
        ele = []
        num = []
        pos = []
        elepos = []
        atomid = []
        j = 0
        for line in f:
            if j >= 2 and j <= 4:
                latvec.append([float(x) for x in line.split()]) 
            if j == 5:
                ele = [str(x) for x in line.split()]
            if j == 6:
                num = [int(x) for x in line.split()]
            if j >= 8 and j <= 7 + sum(num):
                pos.append([float(x) for x in line.split()]) 
                atomid.append(j)
            j += 1
        limatom = 0
        for numa,elea in zip(num,ele):
            for k in range(limatom,limatom+numa):
                elepos.append(elea)
            limatom += numa 
        return np.array(latvec),np.array(ele),np.array(num),np.array(pos),np.array(elepos),np.array(atomid)

def expand_lattice(posa, elea, atomid,nlat):
    expos = []
    exele = []
    exid = []
    for ix in range(-1*nlat[0], nlat[0]):
       for iy in range(-1*nlat[1], nlat[1]):
          for iz in range(-1*nlat[2], nlat[2]):
             for pos, ele, aid in zip(posa,elea, atomid):
                expos.append([ ix + pos[0], iy + pos[1], iz + pos[2] ])
                exele.append(ele)
                exid.append(aid)
    return np.array(expos),np.array(exele), np.array(exid)


def count_bondlength(expos, exele, exid, pos, ele, atomid, lat, el):
    data = []
    cartdata = []
    for e1 in el:
        for e2 in el:
            #tmpdata = [e1,e2]
            print e1, e2
            tmpdata = []
            tmpcartdispos = []
            for pos1, ele1, id1 in zip(pos, ele, atomid):
                for pos2, ele2, id2 in zip(expos, exele, exid):
                    if ele1 == e1 and ele2 == e2:
                        dispos =  pos2 - pos1
                        cartdispos = np.dot(dispos,lat)
                        length = np.linalg.norm(cartdispos)
                        if length <= max_length and length > 0:
                           if ele1 == ele2 and pos2[0] >= 0 and pos2[0] <= 1  and pos2[1] >= 0 and pos2[1] <= 1 and pos2[2] >= 0 and pos2[2] <= 1:
                              if id2 > id1:
                                 tmpdata.append(length)
                                 tmpcartdispos.append(cartdispos)
                           else:
                              tmpdata.append(length) 
                              tmpcartdispos.append(cartdispos)
            data.append(tmpdata)
            cartdata.append(tmpcartdispos)
    return data,cartdata

def rdf(data):
    f = np.zeros(200)
    for dist in data:
        for i in range(1,200):
            x = i * 0.1
            expa = (x - dist)**2 / (2 * sigma**2)
            f[i] += 1/(2*sigma)**0.5 * math.exp(-1.0*expa) / x**2 
    return f 

def stereo(cartposdata,rmin,rmax):
    stereodata=[]
    for cartpos in cartposdata:
        tmpvec = np.array(cartpos)
        length = np.linalg.norm(tmpvec)
        if length >= rmin and length <= rmax:
           tmpvec = tmpvec / length 
           if tmpvec[2] < 0:
              tmpvec = tmpvec * (-1.0)
           tmpvec[2] = tmpvec[2] + 1.0
           tmpvec = tmpvec * 2 / tmpvec[2]
           stereodata.append(tmpvec)
    return np.array(stereodata)

def rdfangle(data):
    f = np.zeros(400)
    for angle,icount in zip(data[:,0],data[:,1]):
        for i in range(0,400):
            x = (i-200) * math.pi / 200.0
            expa = (x - angle)**2 / (2 * asigma**2)
            f[i] += 1/(2*asigma)**0.5 * math.exp(-1.0*expa)*icount 
    return f

def matching(data):
    array = []
    array.append(data)
    dist = []
    intdist = []
    for j in range(0,len(data)):
           dist1 = array[j][0]
           count = 0
           sub = []
           for dist2 in array[j]:
               if abs(dist1 - dist2) < tolerance:
                   count += 1 
               else:
                   sub.append(dist2)
           if len(sub) == 0:  
              break
           else:
              array.append(sub)
              dist.append(dist1)
              intdist.append(count)
 
    distnz1 = np.array(dist)
    intdistnz1 = np.array(intdist)
    distdata = np.c_[distnz1,intdistnz1]
    distdata_sorted = np.array(sorted(distdata, key=lambda distnz1:distnz1[0]))
    return distdata_sorted

def theta(cartposdata, rmin,rmax):
    theta = []
    phi = []
    for cartpos in cartposdata:
        tmpvec = np.array(cartpos)
        length = np.linalg.norm(tmpvec)
        if length >= rmin and length <= rmax:
           tmpvec = tmpvec / length
           if tmpvec[2] < 0:
              tmpvec = tmpvec * (-1.0)
           theta.append(math.acos(tmpvec[2]))
           if tmpvec[0]**2 + tmpvec[1]**2 > 0:
              phi.append(math.acos(tmpvec[0]/(tmpvec[0]**2+tmpvec[1]**2)**0.5))
    return theta,phi 



def caserun(filename,numlat,rmin,rmax):
    lattice, element, numatom, posatom, eleatom, atomid = parse_poscar(filename)
    exposatom, exeleatom, exatomid = expand_lattice(posatom,eleatom, atomid, numlat) 
    distdata,  cartposdata = count_bondlength(exposatom, exeleatom, exatomid, posatom, eleatom, atomid, lattice, element)
    thetaSiSi,phiSiSi = theta(cartposdata[0],rmin,rmax)
    thetadataSiSi = matching(thetaSiSi)
    phidataSiSi = matching(phiSiSi)
    print thetadataSiSi.shape
    print phidataSiSi.shape
    rdfSiSi = rdf(distdata[0])
    rdfthetaSiSi = rdfangle(thetadataSiSi) 
    rdfphiSiSi = rdfangle(phidataSiSi) 
    return(distdata,rdfSiSi,rdfthetaSiSi,rdfphiSiSi) 

def run():

    #call,cSiSi, cSiN, cNN = caserun(ccar,[2,2,3])
    #sall,sSiSi, sSiN, sNN = caserun(scar,[2,2,6])
    call,crdfSiSi,crdfthetaSiSi,crdfphiSiSi = caserun(ccar,[2,2,3],2.0,20)
    sall,srdfSiSi,srdfthetaSiSi,srdfphiSiSi = caserun(scar,[2,2,6],2.0,20)
    #plt.plot(cSiN[:,0],cSiN[:,1]/2,label="cSiN")
    #plt.plot(sSiN[:,0],sSiN[:,1],label="sSiN")
    #plt.plot(cSiSi[:,0],cSiSi[:,1]/2,label="cSiSi")
    #plt.plot(sSiSi[:,0],sSiSi[:,1],label="sSiSi")
    #plt.scatter(cNN[:,0],cNN[:,1]/2,label="cNN")
    #plt.scatter(sNN[:,0],sNN[:,1],label="sNN")
    #plt.bar(x,sySiN)
    rdfx = 0.1*np.arange(0,200)
    ardfx = math.pi*np.arange(-200,200) / 200.0 / math.pi * 180 
    plt.figure(figsize=(6,9))
    plt.subplot(3,1,1)
    plt.plot(rdfx,crdfSiSi/2,label="alpha_SiSi")
    plt.plot(rdfx,srdfSiSi,label="beta_SiSi")
    #plt.hist(call[0],bins=nbins,normed=True,histtype='step',rwidth=0.5,label="alpha_SiSi")
    #plt.hist(sall[0],bins=nbins,normed=True,histtype='step',rwidth=0.5,label="beta_SiSi")
    plt.legend(loc='upper right')
    plt.subplot(3,1,2)
    #plt.plot(rdfx,crdfSiN/2,label="alpha_SiN")
    #plt.plot(rdfx,srdfSiN,label="beta_SiN")
    #plt.hist(call[1],bins=nbins,normed=True,histtype='step',rwidth=0.5,label="alpha_SiN")
    #plt.hist(sall[1],bins=nbins,normed=True,histtype='step',rwidth=0.5,label="beta_SiN")
    plt.legend(loc='upper right')
    plt.subplot(3,1,3)
    #plt.plot(rdfx,crdfNN/2,label="alpha_NN")
    #plt.plot(rdfx,srdfNN,label="beta_NN")
    #plt.hist(call[3],bins=nbins,normed=True,histtype='step',rwidth=0.5,label="alpha_NN")
    #plt.hist(sall[3],bins=nbins,normed=True,histtype='step',rwidth=0.5,label="beta_NN")
    plt.legend(loc='upper right')
    #plt.figure(figsize=(6,6))
    #plt.subplot(2,1,1)
    #plt.bar(x,sySiN, width=0.17)
    #plt.scatter(cstereoSiSi[:,0],cstereoSiSi[:,1],label="alpha_SiSi")
    #plt.scatter(sstereoSiSi[:,0],sstereoSiSi[:,1],label="beta_SiSi")
    #plt.legend(loc='upper right')
    #plt.plot(ardfx,crdfthetaSiSi/2,label="alpha_SiSi_theta")
    #plt.plot(ardfx,srdfthetaSiSi,label="beta_SiSi_theta")
    #plt.legend(loc='upper left')
    #plt.subplot(2,1,2)
    #plt.plot(ardfx,crdfphiSiSi/2,label="alpha_SiSi_phi")
    #plt.plot(ardfx,srdfphiSiSi,label="beta_SiSi_phi")
    #plt.legend(loc='upper left')
    #plt.plot(srdf)

    #print cSiN[0:5,:]
    #print sSiN[0:5,:]


run()
plt.show()
