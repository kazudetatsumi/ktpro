#!/usr/bin/env python
import pickle
import numpy as np

#with open("dumpfile.pkl", 'rb') as f:
f = open("dump.pkl", 'rb')
nt = pickle.load(f)
lossfunc = pickle.load(f)
f.close()
print(nt.params['W1'])
print(lossfunc)
