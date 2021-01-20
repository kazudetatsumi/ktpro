import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def numerical_gradient(f, x):
    h = 1e-4 # 0.0001
    sendbuf = np.zeros_like(x)
    mrank = size - rank - 1
    xshape = x.shape
    if x.ndim == 2:
        if xshape[0] >= mrank + 1:
            for ix in range(mrank + 0, xshape[0] - ((xshape[0]-mrank-1) % size) + size - 1, size):
                for iy in range(0, xshape[1]):
                    tmp_val = x[ix, iy]
                    x[ix, iy] = tmp_val + h
                    fxh1 = f(x)
                    x[ix, iy] = tmp_val - h
                    fxh2 = f(x)
                    sendbuf[ix, iy] = (fxh1 - fxh2) / (2.*h)
                    x[ix, iy] = tmp_val
    if x.ndim == 1:
        if xshape[0] >= mrank + 1:
            for ix in range(mrank + 0, xshape[0] - ((xshape[0]-mrank-1) % size) + size - 1, size):
                    tmp_val = x[ix]
                    x[ix] = tmp_val + h
                    fxh1 = f(x)
                    x[ix] = tmp_val - h
                    fxh2 = f(x)
                    sendbuf[ix] = (fxh1 - fxh2) / (2.*h)
                    x[ix] = tmp_val
    MPI.COMM_WORLD.barrier()
    recvbuf=comm.allgather(sendbuf)
    return np.sum(recvbuf, axis=0)

