#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
import numpy as np
from collections import OrderedDict
from common.layers import *


class convdpNet:
    def __init__(self, input_dim=(1, 1, 2299),
                 conv_params={'filter_num':30, 'filter_height':1, 
                              'filter_width':5, 'pad':0, 'stride':1},
                 hsize=100, osize=4, wstd=0.01):

        FN = conv_params['filter_num']
        FH = conv_params['filter_height']
        FW = conv_params['filter_width']
        P = conv_params['pad']
        S = conv_params['stride']
        C = input_dim[0]
        H = input_dim[1]
        W = input_dim[2]
        if (H + 2*0 - FH) % S != 0:
            print("data along axis 1 is partly not used in convolution layer!")
        if (W + 2*0 - FW) % S != 0:
            print("data along axis 2 is partly not used in convolution layer!")
        OH = (H - FH + 2*P)/S + 1
        OW = (W - FW + 2*P)/S + 1
        PH = 1
        PW = 2
        if (OH + 2*0 - PH) % PH != 0:
            print("data along axis 1 is partly not used in pooling layer!")
        if (OW + 2*0 - PW) % PW != 0:
            print("data along axis 2 is partly not used in pooling layer!")
        POH = int((OH + 2*0 - PH)/PH + 1)
        POW = int((OW + 2*0 - PW)/PW + 1)
        pool_output_size = int(FN * POH * POW)

        self.params = {}
        self.params['W1'] = wstd*np.random.randn(FN, C, FH, FW)
        self.params['b1'] = np.zeros(FN)
        self.params['W2'] = wstd*np.random.randn(pool_output_size, hsize)
        self.params['b2'] = np.zeros(hsize)
        self.params['W3'] = wstd*np.random.randn(hsize, osize)
        self.params['b3'] = np.zeros(osize)

        self.layers = OrderedDict()
        self.layers['Conv1'] = Convolution(self.params['W1'],
                                           self.params['b1'],
                                           conv_params['stride'],
                                           conv_params['pad'])
        self.layers['Relu1'] = Relu()
        self.layers['Pool1'] = Pooling(PH=1, PW=2, S=2, P=0)
        self.layers['Affine1'] = Affine(self.params['W2'], self.params['b2'])
        self.layers['Relu2'] = Relu()
        self.layers['Dropout1'] = Dropout()
        self.layers['Affine2'] = Affine(self.params['W3'], self.params['b3'])
        self.layers['Dropout2'] = Dropout()

        self.lastLayer = SoftmaxWithLoss()

    def predict(self, x):
        for layer in self.layers.values():
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

    def analytical_grad(self, x, t):
        # forward
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

