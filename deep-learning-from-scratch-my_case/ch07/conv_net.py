#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
import numpy as np
from collections import OrderedDict
from common.layers import *
from common.gradient import numerical_gradient


class convNet:
    def __init__(self, input_dim=(1, 28, 28), 
            conv_params={'filter_num':30, 'filter_size':5, 'pad':0, 'stride':1},
            hsize = 100, osize = 10, wstd=0.01):

        filter_num = conv_params['filter_num']
        filter_size = conv_params['filter_size']
        filter_pad = conv_params['pad']
        filter_stride = conv_params['stride']
        input_size = input_dim[1]
        conv_output_size = (input_size - filter_size + 2*filter_pad)/filter_stride + 1
        pool_output_size = int(filter_num * ((conv_output_size + 2*0 - 2)/2 + 1) * ((conv_output_size + 2*0 - 2)/2 + 1))  # pad_size is fixed as 2.

        self.params = {}
        self.params['W1'] = wstd*np.random.randn(filter_num, input_dim[0], filter_size, filter_size)    # FN, C, FH, FW
        self.params['b1'] = np.zeros(filter_num)                                                                            # FN
        self.params['W2'] = wstd*np.random.randn(pool_output_size, hsize)                             # for affine from the output of pooling
        self.params['b2'] = np.zeros(hsize)
        self.params['W3'] = wstd*np.random.randn(hsize, osize)
        self.params['b3'] = np.zeros(osize)
        
        ## layers are ordered as Affine-Relu-Affine-SoftmaxWithLoss
        self.layers = OrderedDict()
        self.layers['Conv1'] = Convolution(self.params['W1'], self.params['b1'], conv_params['stride'], conv_params['pad'])
        self.layers['Relu1'] = Relu()
        self.layers['Pool1'] = Pooling(PH=2, PW=2, S=2, P=0)
        self.layers['Affine1'] = Affine(self.params['W2'], self.params['b2'])
        self.layers['Relu2'] = Relu()
        self.layers['Affine2'] = Affine(self.params['W3'], self.params['b3'])

        self.lastLayer = SoftmaxWithLoss()

    def predict(self, x):
        for layer in self.layers.values():
            print(layer)
            print("x", x.shape)
            x = layer.forward(x)
        return x

    def loss(self, x, t):
        y = self.predict(x)
        return self.lastLayer.forward(y, t)

    def acc(self, x, t):
        y = self.predict(x)
        if t.ndim != 1:
            t = np.argmax(t, axis=1)
        acc = np.sum(np.argmax(y, axis=1) == t) / float(t.shape[0])
        return acc

    #def numerical_grad(self, x, t):
    #    def f(dummy): return self.loss(x, t)
    #    grads = {}
    #    grads['W1'] = numerical_gradient(f, self.params['W1'])
    #    grads['b1'] = numerical_gradient(f, self.params['b1'])
    #    grads['W2'] = numerical_gradient(f, self.params['W2'])
    #    grads['b2'] = numerical_gradient(f, self.params['b2'])
    #    return grads

    def analytical_grad(self, x, t):
        # forward   Is this necessary? -> Yes. It conveys t.
        self.loss(x, t)

        # backward
        dout = 1
        dout = self.lastLayer.backward(dout)
        layers = list(self.layers.values())
        layers.reverse()
        for layer in layers:
            dout = layer.backward(dout)

        grads = {}
        grads['W1'] = self.layers['Conv1'].dW
        grads['b1'] = self.layers['Conv1'].db
        grads['W2'] = self.layers['Affine1'].dW
        grads['b2'] = self.layers['Affine1'].db
        grads['W3'] = self.layers['Affine2'].dW
        grads['b3'] = self.layers['Affine2'].db

        return grads

