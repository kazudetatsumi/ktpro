#!/usr/bin/env python
import sys
import os
import numpy as np
sys.path.append(os.pardir)
from dataset.sin import load_sin
from two_layer_net_reg import twoLNet
from common.optimizer import Adagrad
import pickle
import matplotlib.pyplot as plt
import math


def run():
    num_iters = 80000
    batch_size = 10
    lossfunc = []
    loss_test = []
    loss_train = []

    (x_train, t_train), (x_test, t_test) = load_sin()
    ## normalization: By normalization, the parameters are optimized over the whole data range.
    ave_x_train = np.average(x_train)
    ave_t_train = np.average(t_train)
    var_x_train = np.var(x_train)
    var_t_train = np.var(t_train)
    x_train = (x_train - ave_x_train) / np.sqrt(var_x_train)
    t_train = (t_train - ave_t_train) / np.sqrt(var_t_train)
    x_test = (x_test - ave_x_train) / np.sqrt(var_x_train)
    t_test = (t_test - ave_t_train) / np.sqrt(var_t_train)

    train_size = x_train.shape[0]
    nt = twoLNet(dsize=1, hsize=12,  osize=1)
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
        #if lidx % itersperepoch == 0:
        #    loss_test.append(nt.loss(x_test, t_test))
        #    loss_train.append(lossfunc[-1])

    
    predicted_t_test = (nt.predict(x_test) * np.sqrt(var_t_train) + ave_t_train) # denormaleze the predicted 
    x_test = x_test * np.sqrt(var_x_train) + ave_x_train                         # denormalize the x_test
    x_train = x_train * np.sqrt(var_x_train) + ave_x_train # denormaleze the predicted 
    t_train = t_train * np.sqrt(var_t_train) + ave_t_train                         # denormalize the x_test

    plt.scatter(x_test.reshape(-1)[:], predicted_t_test.reshape(-1)[:], marker='x')
    plt.scatter(x_train.reshape(-1)[:], t_train.reshape(-1)[:], marker='+', c='red', s=50)
    xref = np.arange(0, 3.0*math.pi, 0.01)
    plt.plot(xref, np.sin(xref))  
    #plt.scatter(x_batch.reshape(-1), t_batch.reshape(-1))
    #y_train = nt.predict(x_train)
    #y_test  = nt.predict(x_test)
    #print(x_train[0:10])
    #print(y_train[0:10])

    pb = open('dump_reg.pkl', 'wb')
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
