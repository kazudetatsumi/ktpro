#!/usr/bin/env python 
import  mpi4py
mpi4py.rc(initialize=False, finalize=False)
from mpi4py import MPI
MPI.Init()

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

print("Hello from rank %s of %s." % (rank, size))
#print(MPI.Is_initialized())
#print(MPI.Is_finalized())

#MPI.COMM_WORLD.Barrier()
MPI.Finalize()
print("Is Finalized?")

#print(MPI.Is_initialized())
#print(MPI.Is_finalized())

