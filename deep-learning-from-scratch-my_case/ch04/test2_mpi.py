#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.mnist import load_mnist
from two_layer_net_mpi import twoLNet
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def run():
    num_iters = 10000
    batch_size = 100
    learning_rate = 0.1
    lossfunc = []
    acc_test = []
    acc_train = []
    (x_train, t_train), (x_test, t_test) = load_mnist(
                                                      normalize=True,
                                                      one_hot_label=True
                                                      )
    train_size = x_train.shape[0]
    nt = twoLNet(dsize=784, hsize=50, osize=10)
    comm.Bcast(nt.params['W1'], root=0)
    comm.Bcast(nt.params['W2'], root=0)
    itersperepoch = train_size // batch_size
    if rank == 0:
        print("number of train data: {}".format(train_size))
        print("number of iteraions per epoch: {}".format(itersperepoch))

    for lidx in range(num_iters):
        if rank == 0:
            batch_mask = np.random.choice(train_size, batch_size)
        else:
            batch_mask = np.zeros_like(np.random.choice(train_size, batch_size))
        comm.Bcast(batch_mask, root=0)
        x_batch = x_train[batch_mask]
        t_batch = t_train[batch_mask]
        grad = nt.numerical_grad(x_batch, t_batch)
        for key in {'W1', 'W2', 'b1', 'b2'}:
            nt.params[key] -= learning_rate * grad[key]
        if rank == 0:
            print("iteration# {} finished, lossfunc: {}".format(lidx,nt.loss(x_batch, t_batch)))
        lossfunc.append(nt.loss(x_batch, t_batch))
        if lidx % itersperepoch == 0:
            acc_test.append(nt.acc(x_test, t_test))
            acc_train.append(nt.acc(x_train, t_train))
    if rank == 0:
        #plt.plot(lossfunc)
        #plt.savefig("lossfunc.png")
        pb = open('dump.pkl', 'wb')
        pickle.dump(nt, pb)
        pickle.dump(lossfunc, pb)
        pickle.dump(acc_test, pb)
        pickle.dump(acc_train, pb)
        pb.close()

    #plt.show()


run()
