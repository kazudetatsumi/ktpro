#!/usr/bin/env python
from mpi4py import MPI 
import numpy as np 

comm = MPI.COMM_WORLD 

size = 10
itemsize = MPI.DOUBLE.Get_size() 
print(itemsize)
if comm.Get_rank() == 0: 
    nbytes = size * itemsize 
else: 
    nbytes = 0 
win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm) 

buf, itemsize = win.Shared_query(0) 
assert itemsize == MPI.DOUBLE.Get_size() 
buf = np.array(buf, dtype='B', copy=False) 
ary = np.ndarray(buffer=buf, dtype='d', shape=(size,)) 

if comm.rank == 1: 
    ary[:5] = np.arange(5) 
comm.Barrier() 
if comm.rank == 0: 
    print(ary[:10])
