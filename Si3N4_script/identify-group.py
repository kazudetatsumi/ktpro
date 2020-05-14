#!/usr/bin/env python
from operator import itemgetter, attrgetter
f=open("/home/kazu/bsi3n4_m/phono3py_113_fc2_338_sym_monk/irreps.yaml",'r')
lines=[]
cha=[]
freq=[]
indices=[]
data=[]
for l in f:
    lines.append(l.replace("360","  0"))
f.close()

for l in lines:
    if "band_indices" in l:
        indices.append(l[:-1])
    if "frequency" in l:
        freq.append(l[:-1])
    if "characters" in l:
        cha.append(l[:-1])

print len(indices)
print len(freq)
print len(cha)

for i,f,c in zip(indices,freq,cha):
    data.append( (i,f,c) )

data_sorted=sorted(data,key=itemgetter(2))


for d in data_sorted:
    print d[1]
    


