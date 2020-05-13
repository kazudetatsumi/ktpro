#!/usr/bin/env python
import numpy as np
import time
n = 100000000
x = np.ones(n)
alpha = 0.01
start_time = time.perf_counter()
x = x * alpha
elapsed_time = time.perf_counter() - start_time 
print ("elapsed_time:{0}".format(elapsed_time))
