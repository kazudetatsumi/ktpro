#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import re

lines=[]
f = open('run2.dos1ev','r')
i=0
for l in f:
	lines.append(l.rstrip())
	i=i+1

endnum=i
start=3
num=endnum-start
ene=np.zeros((num))
dos=np.zeros((num))
for j in range(start,endnum):
	values=re.split(" +",lines[j])
	ene[j-start]=values[1]
	dos[j-start]=values[2]

intdos=0
for j in range(0,num):
	if ene[j] <= 0:
		print ene[j],dos[j]
		intdos=intdos+dos[j]*(ene[j+1]-ene[j])
#plt.figure(figsize=(6,9))
#plt.plot(ene,dos)

print intdos
plt.show()

