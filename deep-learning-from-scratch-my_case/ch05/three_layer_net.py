#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
import numpy as np
from collections import OrderedDict
from common.layers import *
from common.gradient import numerical_gradient


class threeLNet:
    def __init__(self, dsize, hsize1, hsize2, osize, wstd=0.1):

        self.params = {}
        self.params['W1'] = wstd * np.random.randn(dsize, hsize1)
        self.params['b1'] = np.zeros(hsize1)
        self.params['W2'] = wstd * np.random.randn(hsize1, hsize2)
        self.params['b2'] = np.zeros(hsize2)
        self.params['W3'] = wstd * np.random.randn(hsize2, osize)
        self.params['b3'] = np.zeros(osize)
        
        ## layers are ordered as Affine-Relu-Affine-SoftmaxWithLoss
        self.layers = OrderedDict()
        self.layers['Affine1'] = Affine(self.params['W1'], self.params['b1'])
        self.layers['Relu1'] = Relu()
        self.layers['Affine2'] = Affine(self.params['W2'], self.params['b2'])
        self.layers['Relu2'] = Relu()
        self.layers['Affine3'] = Affine(self.params['W3'], self.params['b3'])

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
        acc = np.sum(np.argmax(y, axis=1) == t) / t.shape[0] * 1.
        return acc

    def numerical_grad(self, x, t):
        def f(dummy): return self.loss(x, t)
        grads = {}
        grads['W1'] = numerical_gradient(f, self.params['W1'])
        grads['b1'] = numerical_gradient(f, self.params['b1'])
        grads['W2'] = numerical_gradient(f, self.params['W2'])
        grads['b2'] = numerical_gradient(f, self.params['b2'])
        grads['W3'] = numerical_gradient(f, self.params['W3'])
        grads['b3'] = numerical_gradient(f, self.params['b3'])
        return grads

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
        grads['W1'] = self.layers['Affine1'].dW
        grads['b1'] = self.layers['Affine1'].db
        grads['W2'] = self.layers['Affine2'].dW
        grads['b2'] = self.layers['Affine2'].db
        grads['W3'] = self.layers['Affine3'].dW
        grads['b3'] = self.layers['Affine3'].db

        return grads

