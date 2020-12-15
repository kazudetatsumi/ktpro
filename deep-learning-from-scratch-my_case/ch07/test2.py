#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.mnist import load_mnist
from conv_net import convNet
from common.optimizer import Adam
import pickle


def run():
    num_iters = 10
    batch_size = 100
    lossfunc = []
    acc_test = []
    acc_train = []

    (x_train, t_train), (x_test, t_test) = load_mnist(
                                                      normalize=True,
                                                      flatten=False,
                                                      one_hot_label=False
                                                      )
    print(x_train.shape)
    train_size = x_train.shape[0]
    nt = convNet(input_dim=(1, 28, 28),
            conv_params={'filter_num':30, 'filter_size':5, 'pad':0, 'stride':1},
            hsize = 100, osize = 10, wstd=0.01)

    opt = Adam(lr=0.001)
    itersperepoch = max(train_size // batch_size, 1)
    for lidx in range(num_iters):
        batch_mask = np.random.choice(train_size, batch_size)
        x_batch = x_train[batch_mask]
        t_batch = t_train[batch_mask]
        #grad = nt.numerical_grad(x_batch, t_batch)
        grad = nt.analytical_grad(x_batch, t_batch)
        for key in grad.keys():
            params = nt.params
            opt.update(params, grad)
        print("iteration# {} finished, lossfunc: {}".format(lidx,nt.loss(x_batch, t_batch)))
        lossfunc.append(nt.loss(x_batch, t_batch))
        if lidx % itersperepoch == 0:
            acc_test.append(nt.acc(x_test[:1000], t_test[:1000]))
            acc_train.append(nt.acc(x_train[:1000], t_train[:1000]))
            #print(nt.acc(x_test[:1000], t_test[:1000]))
            #print(nt.acc(x_train[:1000], t_train[:1000]))

    #pb = open('dump2.pkl', 'wb')
    #pickle.dump(nt, pb)
    #pickle.dump(lossfunc, pb)
    #pickle.dump(acc_test, pb)
    #pickle.dump(acc_train, pb)
    #pb.close()


run()
