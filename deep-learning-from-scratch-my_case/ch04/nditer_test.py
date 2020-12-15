#!/usr/bin/env python
import numpy as np
a = np.arange(10).reshape(2,5)
#for x in np.nditer(a):
#    print(x, end=" ")

it = np.nditer(a, flags=['multi_index'])
for x in it:
    print(x,"@", it.multi_index)
