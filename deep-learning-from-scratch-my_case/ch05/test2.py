#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.mnist import load_mnist
from two_layer_net import twoLNet
import pickle


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
    itersperepoch = train_size // batch_size
    for lidx in range(num_iters):
        batch_mask = np.random.choice(train_size, batch_size)
        x_batch = x_train[batch_mask]
        t_batch = t_train[batch_mask]
        #grad = nt.numerical_grad(x_batch, t_batch)
        grad = nt.analytical_grad(x_batch, t_batch)
        for key in grad.keys():
            nt.params[key] -= learning_rate * grad[key]
        print("iteration# {} finished, lossfunc: {}".format(lidx,nt.loss(x_batch, t_batch)))
        lossfunc.append(nt.loss(x_batch, t_batch))
        if lidx % itersperepoch == 0:
            acc_test.append(nt.acc(x_test, t_test))
            acc_train.append(nt.acc(x_train, t_train))

    #pb = open('dump.pkl', 'wb')
    #pickle.dump(nt, pb)
    #pickle.dump(lossfunc, pb)
    #pickle.dump(acc_test, pb)
    #pickle.dump(acc_train, pb)
    #pb.close()


run()
