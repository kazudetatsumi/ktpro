#!/usr/bin/env python
from mpi4py import MPI
comm = MPI.COMM_WORLD 
rank = comm.Get_rank()
size = comm.Get_size()
print("Hello World at rank {0}/{1}".format(rank, size))
