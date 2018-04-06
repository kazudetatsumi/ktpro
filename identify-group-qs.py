#!/usr/bin/env python
from operator import itemgetter, attrgetter

dirname="/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/"


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
   gid=0
   for d in data:
       if d[3] != origd: 
           origd = d[3]
           gid += 1
       values=d[0].split()
       qz=float(values[4])
       values=d[2].split()
       omega=float(values[1])
       data_assigned.append((qz,omega,gid))
   return(data_assigned)

def run():
   alldata=[]
   for i in [1,2,3,4]:
      filename=dirname+"/irreps-"+str(i)+".yaml"
      alldata.extend(sortirreps(filename))
   for d in alldata:
      print d

    

run()
