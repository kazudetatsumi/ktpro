#!/usr/bin/env python
from ctypes import *
import numpy as np 
import time

fmath = np.ctypeslib.load_library("fmath.so",".")
fmath.scal.argtypes = [
                       np.ctypeslib.ndpointer(dtype=np.float64),
                       POINTER(c_double), POINTER(c_int64)
                       ] 
fmath.scal.restype = c_void_p
n = 100000000
x = np.ones(n, dtype=np.float64) 
alpha = c_double(0.01)
len = c_int64(x.size)
start_time = time.perf_counter()
fmath.scal(x, alpha, len)
elapsed_time = time.perf_counter() - start_time 
print ("elapsed_time:{0}".format(elapsed_time))
