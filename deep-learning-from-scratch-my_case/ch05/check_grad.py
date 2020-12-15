#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
import numpy as np
from dataset.mnist import load_mnist
from two_layer_net import twoLNet


def run():
    (x_train, t_train), (x_test, t_test) = load_mnist(
                                                      normalize=True,
                                                      one_hot_label=True
                                                      )
    nt = twoLNet(dsize=784, hsize=50, osize=10)
    x_batch = x_train[:3]
    t_batch = t_train[:3]

    grads_numerical = nt.numerical_grad(x_batch, t_batch)
    grads_analytical = nt.analytical_grad(x_batch, t_batch)

    for key in grads_numerical.keys():
        diff = np.average(np.abs(grads_numerical[key] - grads_analytical[key]))
        print("mean absolute difference in two grad wrt {}: {}".format(key, diff))


run()


