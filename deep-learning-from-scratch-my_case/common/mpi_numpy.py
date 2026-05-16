#!/usr/bin/env python
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

sendbuf = np.zeros((10, 20), dtype='i') 
sendbuf[rank,0] = 1
#recvbuf = None
#recvbuf = np.empty([size, 10, 20], dtype='i')
recvbuf = np.array(comm.allgather(sendbuf))
if rank == 0:
    print(recvbuf)
    print(recvbuf.shape)
    print(np.sum(recvbuf, axis=0))
    print(np.sum(recvbuf, axis=0).shape)

    #for i in range(size):
    #    assert np.allclose(recvbuf[i,:], i)



