#!/usr/bin/env python
# This python script collects representative atom index from phonpy-log firstly,
# then creates "BORN" file from OUTCAR file containg born effective charge info.

import numpy as np
import sys
argvs = sys.argv

indx=[]
f = open('phonopy-log')
AtomFound=0
superFound=0

for l in f:
	if "Atomic positions" in l:
		AtomFound=1
	if "super cell" in l:
		superFound=1
	if AtomFound == 1 and superFound == 0:
		if "*" in l:
			values=l.split()
			#indx[na]=values[0]
			indx.append(int(values[0].lstrip('*')))
print "index number of representative atoms:"
for na in indx:
	print na
f.close()

lines=[]
g = open(argvs[1],'r')
for l in g:
	lines.append(l)
g.close()

dieFound=0

dca=np.zeros(9)
bca=np.zeros((len(indx),9))
linenumber=0
g = open(argvs[1],'r')
for l in g:
	linenumber=linenumber+1
	if "MACROSCOPIC STATIC DIELECTRIC TENSOR" in l and dieFound == 0:
		print linenumber
		print len(lines)
		print lines[linenumber+1]
		for i in range(1,4):
			values=lines[linenumber+i].split()
			print len(values)
			for j in range(0,3):
				dca[j+(i-1)*3]=float(values[j])
		dieFound=1
		for i in range(1,4):
			print i,"fuck"
	if "BORN EFFECTIVE CHARGES" in l:
		k=0
		for i in indx:
			print lines[linenumber+(i-1)*4+1]
			print lines[linenumber+(i-1)*4+2]
			print lines[linenumber+(i-1)*4+3]
			print lines[linenumber+(i-1)*4+4]
			values=lines[linenumber+(i-1)*4+2].split()
			bca[k,0]=float(values[1])
			bca[k,1]=float(values[2])
			bca[k,2]=float(values[3])
			values=lines[linenumber+(i-1)*4+3].split()
			bca[k,3]=float(values[1])
			bca[k,4]=float(values[2])
			bca[k,5]=float(values[3])
			values=lines[linenumber+(i-1)*4+4].split()
			bca[k,6]=float(values[1])
			bca[k,7]=float(values[2])
			bca[k,8]=float(values[3])
			k=k+1
			
print dca
print bca

h=open('born.txt','w')
h.write("14.400 \n")
h.close()
with open('born.txt','a') as fuck:
	np.savetxt(fuck,dca.reshape(1,9), fmt='%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f')
	for i in range(0,len(indx)):
		np.savetxt(fuck,bca[i,:].reshape(1,9), fmt='%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f')



