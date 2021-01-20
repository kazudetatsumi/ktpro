#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
import numpy as np
from common.functions import sigmoid, softmax, cross_entropy_error
from common.gradient_mpi import numerical_gradient



class twoLNet:
    def __init__(self, dsize, hsize, osize, wstd=0.1):
        self.params = {}
        self.params['W1'] = wstd * np.random.randn(dsize, hsize)
        self.params['b1'] = np.zeros(hsize)
        self.params['W2'] = wstd * np.random.randn(hsize, osize)
        self.params['b2'] = np.zeros(osize)

    def y(self, x):
        W1, W2, b1, b2 = self.params['W1'], self.params['W2'], \
                         self.params['b1'], self.params['b2']
        a1 = np.dot(x, W1) + b1
        z1 = sigmoid(a1)
        a2 = np.dot(z1, W2) + b2
        y = softmax(a2)
        return y

    def loss(self, x, t):
        y = self.y(x)
        return cross_entropy_error(y, t)

    def acc(self, x, t):
        y = self.y(x)
        acc = np.sum(np.argmax(y, axis=1) == np.argmax(t, axis=1)) /\
            t.shape[0] * 1.
        return acc

    def numerical_grad(self, x, t):
        def f(dummy): return self.loss(x, t)
        grads = {}
        grads['W1'] = numerical_gradient(f, self.params['W1'])
        grads['b1'] = numerical_gradient(f, self.params['b1'])
        grads['W2'] = numerical_gradient(f, self.params['W2'])
        grads['b2'] = numerical_gradient(f, self.params['b2'])
        return grads


