#!/usr/bin/env python
#coding: UTF-8
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('classic')
w=0.4
xdata1=np.arange(1,17)
labels=["personality", "originality","language ability","expertise","independence","toughness","team work ryoku", "logical thinking", "calcs", "PC skill", "business manner","common sense", "common culture", "communication","others"] 
xdata2=xdata1+w
ydata1=[3.5,5.5,0.4,1.0,20.4,5.5,13.3,4.5,4.8,0.1,0.2,3.8,11,3.5,10,1.2]
ydata2=[3.8,7.6,16.5,11.8,5.6,3.6,3.0,2.3,6.1,10.2,5.7,6.2,5.8,3.1,8.0,0.7]
plt.figure(figsize=(8,4))

plt.bar(xdata1,ydata1,width=w,label="company",color="red")
plt.bar(xdata2,ydata2,width=w,label="stu")
#plt.xlabel("Omega [THz]")
#plt.ylabel("DOS [1/THz]")
plt.xticks(xdata1+w, labels,rotation=90)
plt.ylim(0,22)
plt.xlim(1,17)
plt.yticks([0,10,20])
plt.legend(loc='upper right')
#plt.show()
plt.savefig("gakusei-kigyo.eps")
