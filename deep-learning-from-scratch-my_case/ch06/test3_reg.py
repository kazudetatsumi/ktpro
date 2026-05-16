#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.cos_and_linear import load_cos
from three_layer_net_reg import threeLNet
from common.optimizer import Adagrad
import pickle
import matplotlib.pyplot as plt
import math


def run():
    num_iters = 10000
    batch_size = 200
    lossfunc = []
    loss_test = []
    loss_train = []

    (x_train, t_train), (x_test, t_test) = load_cos()
    #y_train = np.zeros_like(x_train)
    #y_test = np.zeros_like(x_test)

    train_size = x_train.shape[0]
    nt = threeLNet(dsize=1, hsize1=12, hsize2=12,  osize=1)
    opt = Adagrad()
    itersperepoch = train_size // batch_size
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
            loss_test.append(nt.loss(x_test, t_test))
            loss_train.append(lossfunc[-1])
        #    acc_test.append(nt.acc(x_test, t_test))
        #    acc_train.append(nt.acc(x_train, t_train))

    plt.scatter(x_test.reshape(-1)[:], nt.predict(x_test).reshape(-1)[:], marker='x')
    xref = np.arange(0, 4.0*math.pi, 0.01)
    plt.plot(xref, np.cos(xref)+0.05*xref)  
    #plt.scatter(x_batch.reshape(-1), t_batch.reshape(-1))
    #y_train = nt.predict(x_train)
    #y_test  = nt.predict(x_test)
    #print(x_train[0:10])
    #print(y_train[0:10])

    pb = open('dump3_reg.pkl', 'wb')
    pickle.dump(nt, pb)
    pickle.dump(lossfunc, pb)
    pickle.dump(loss_test, pb)
    pickle.dump(loss_train, pb)
    #pickle.dump(x_train, pb)
    #pickle.dump(y_train, pb)
    #pickle.dump(x_test, pb)
    #pickle.dump(y_test, pb)
    pb.close()
    plt.show()


run()
