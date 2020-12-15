#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
import numpy as np
from collections import OrderedDict
from common.layers import *
from common.gradient import numerical_gradient


class twoLNet:
    def __init__(self, dsize, hsize, osize, wstd=0.01):

        self.params = {}
        #self.params['W1'] = wstd * np.random.randn(dsize, hsize)
        self.params['W1'] = np.random.randn(dsize, hsize) / np.sqrt(dsize)
        self.params['b1'] = np.zeros(hsize)
        #self.params['W2'] = wstd * np.random.randn(hsize, osize)
        self.params['W2'] = np.random.randn(hsize, osize) / np.sqrt(hsize)
        self.params['b2'] = np.zeros(osize)
        
        ## layers are ordered as Affine-Tanh-Affine-IdentityWithLoss
        self.layers = OrderedDict()
        self.layers['Affine1'] = Affine(self.params['W1'], self.params['b1'])
        self.layers['Tanh1'] = Tanh()
        self.layers['Affine2'] = Affine(self.params['W2'], self.params['b2'])

        self.lastLayer = IdentityWithLoss()

    def predict(self, x):
        for layer in self.layers.values():
            x = layer.forward(x)
        return x

    def loss(self, x, t):
        y = self.predict(x)
        return self.lastLayer.forward(y, t)

    def numerical_grad(self, x, t):
        def f(dummy): return self.loss(x, t)
        grads = {}
        grads['W1'] = numerical_gradient(f, self.params['W1'])
        grads['b1'] = numerical_gradient(f, self.params['b1'])
        grads['W2'] = numerical_gradient(f, self.params['W2'])
        grads['b2'] = numerical_gradient(f, self.params['b2'])
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

        return grads

