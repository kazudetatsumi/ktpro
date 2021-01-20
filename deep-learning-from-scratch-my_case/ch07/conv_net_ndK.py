#!/usr/bin/env python
import sys
import os
sys.path.append(os.pardir)
import numpy as np
from collections import OrderedDict
from common.layers import *


class convNetK:
    def __init__(self, input_dim=(1, 1, 2299),
            conv_params={'FN':30, 'FH':1, 'FW':5, 'P':0, 'S':1},
            pool_params={'PH':1, 'PW':5, 'P':0, 'S':1},
            hsize=100, osize=4, wstd=0.01):

        FN = conv_params['FN']
        FH = conv_params['FH']
        FW = conv_params['FW']
        P = conv_params['P']
        S = conv_params['S']
        C = input_dim[0]
        H = input_dim[1]
        W = input_dim[2]
        OW = (W - FW + 2*P)/S + 1
        PH = pool_params['PH']
        PW = pool_params['PW']
        POW = int((OW + 2*pool_params['P'] - PW)/PW + 1)
        pool_output_size = int(FN * POW)

        self.params = {}
        self.params['W1'] = (wstd*np.random.randn(FN, C, FH, FW)).astype(np.float32)    # FN, C, FH, FW
        self.params['b1'] = np.zeros(FN, dtype='float32')                                                               # FN
        self.params['W2'] = (wstd*np.random.randn(pool_output_size, hsize)).astype(np.float32)              
        self.params['b2'] = np.zeros(hsize, dtype='float32')
        self.params['W3'] = (wstd*np.random.randn(hsize, osize)).astype(np.float32)
        self.params['b3'] = np.zeros(osize, dtype='float32')
        
        self.layers = OrderedDict()
        self.layers['Conv1'] = ConvolutionK(self.params['W1'], self.params['b1'], conv_params['S'], conv_params['S'], 0, conv_params['P'])
        self.layers['Relu1'] = Relu()
        self.layers['Pool1'] = PoolingK(PH=PH, PW=PW, Sh=PH, Sw=PW, Ph=0, Pw=pool_params['P'])
        self.layers['Affine1'] = Affine(self.params['W2'], self.params['b2'])
        self.layers['Relu2'] = Relu()
        self.layers['Affine2'] = Affine(self.params['W3'], self.params['b3'])

        self.lastLayer = SoftmaxWithLoss()

    def predict(self, x):
        for layer in self.layers.values():
            #print("x", x.shape)
            #print(layer)
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

