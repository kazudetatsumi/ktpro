#!/usr/bin/env python
from operator import itemgetter, attrgetter
import matplotlib.pyplot as plt
import numpy as np

#dirname="/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/"
dirname="/home/kazu/asi3n4/phono3py_112_fc2_334_sym_monk_shift/"
#dirname="/home/kazu/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift"
inifile_no=1
finfile_no=4


def sortirreps(fname):
   f=open(fname,'r')
   lines=[]
   cha=[]
   freq=[]
   indices=[]
   data=[]
   for l in f:
       lines.append(l.replace("360","  0"))
   f.close()

   for l in lines:
       if "q-position" in l:
           qpos=l[:-1]
       if "band_indices" in l:
           indices.append(l[:-1])
       if "frequency" in l:
           freq.append(l[:-1])
       if "characters" in l:
           cha.append(l[:-1])
   

   for i,f,c in zip(indices,freq,cha):
       data.append( (qpos,i,f,c) )

   data_sorted=sorted(data,key=itemgetter(3))
   data_final=assign_no(data_sorted)
   return(data_final)



def assign_no(data):
   data_assigned=[]
   origd=data[0][3]
   values=[]
   values=origd.replace("[","").replace("]","").replace(",","").split()
   i=0
   mag=[]
   ph=[]
   for v in values[1:]:
      if i % 2 == 0:
          mag.append(v)
      else:
          ph.append(v)
      i+=1
   origmag=mag
   origph=ph
   gid=0
   for d in data:
       values=[]
       mag=[]
       ph=[]
       i=0
       values=d[3].replace("[","").replace("]","").replace(",","").split()
       for v in values[1:]:
           if i % 2 == 0:
               mag.append(v)
           else:
               ph.append(v)
           i+=1

       #if d[3] != origd: 
       if mag != origmag:                        
           origmag = mag
           origph  = ph
           gid += 1
       else:
           difffound = 0
           for om,op,p in zip(origmag,origph,ph):
               if op !=p and om != "0":
                   difffound +=1
           if difffound > 0:
               #print origmag
               #print mag
               #print origph
               #print ph
               #print gid
               origmag = mag
               origph  = ph
               gid += 1

       values=d[0].split()
       qz=float(values[4])
       qy=float(values[3][:-1])
       values=d[2].split()
       omega=float(values[1])
       #print d[2], d[3]
       data_assigned.append((qy,qz,omega,gid))
       data_assigned.append((qy,qz,omega,gid))
   return(data_assigned)




def run():
   alldata=[]
   freqs=[]
   qs=[]
   gids=[]
   for i in range(inifile_no,finfile_no+1):
      filename=dirname+"/irreps-"+str(i)+".yaml"
      alldata.extend(sortirreps(filename))
   for d in alldata:
         print d
      #if d[2] == 2:
         qs.append(d[1])
         freqs.append(d[2])
         gids.append(d[3])
   freqall=np.array(freqs)
   qsall=np.array(qs)
   gidsall=np.array(gids)
   plt.figure(figsize=(6,12))
   plt.scatter(qs,freqs,c=gids,marker='x',linewidth=0.4,s=15, cmap='jet')
   plt.xlim(0.000,0.50)
   
   

run()
#plt.show()
plt.savefig("bands_irreps.eps")
