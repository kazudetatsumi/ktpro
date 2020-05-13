#!/usr/bin/env python
from ctypes import * 
import numpy as np
from mpi4py import MPI
fmath = np.ctypeslib.load_library("fmath.so",".") 
fmath.scal.argtypes = [
                      np.ctypeslib.ndpointer(dtype=np.float64),
                      POINTER(c_double), POINTER(c_int64)
                      ] 
fmath.scal.restype = c_void_p
comm = MPI.COMM_WORLD 
print(comm)
size = comm.Get_size()
print(size)
n = int(10000000/size)
x = np.ones(n, dtype=np.float64) 
alpha = c_double(0.01)
len = c_int64(x.size) 
fmath.scal(x, alpha, len)
