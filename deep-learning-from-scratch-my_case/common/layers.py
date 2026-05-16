import numpy as np
from common.functions import sigmoid, cross_entropy_error, softmax,\
                             mean_squared_error
from common.util import im2col, col2im, im2colK, col2imK


class SoftmaxWithLoss:
    def __init__(self):
        self.loss = None
        self.y = None
        self.t = None

    def forward(self, x, t):
        self.t = t
        self.y = softmax(x)
        self.loss = cross_entropy_error(self.y, self.t)
        return self.loss

    def backward(self, dout=1):
        batch_size = self.t.shape[0]
        if self.t.size == self.y.size:
            dx = (self.y - self.t) / batch_size
        else:
            dx = self.y.copy()
            dx[np.arange(batch_size), self.t] -= 1
            dx = dx / batch_size

        return dx


class IdentityWithLoss:
    def __init__(self):
        self.loss = None
        self.y = None
        self.t = None

    def forward(self, x, t):
        self.t = t
        self.y = x
        self.loss = mean_squared_error(self.y, self.t)
        return self.loss

    def backward(self, dout=1):
        batch_size = self.t.shape[0]
        dx = (self.y - self.t) / batch_size
        return dx


class Tanh:
    def __init__(self):
        self.out = None

    def forward(self, x):
        out = np.tanh(x)
        self.out = out
        return out

    def backward(self, dout):
        dx = dout * (1. - self.out**2)
        return dx


class LeakyRelu:
    def __init__(self):
        self.mask = None
        self.alpha = None

    def forward(self, x):
        self.mask = x <= 0
        self.alpha = 0.01
        out = x.copy()
        out[self.mask] = self.alpha * x[self.mask]
        return out

    def backward(self, dout):
        dout[self.mask] *= self.alpha
        dx = dout
        return dx


class Sigmoid:
    def __init__(self):
        self.out = None

    def forward(self, x):
        out = sigmoid(x)
        self.out = out
        return out

    def backward(self, dout):
        dx = dout * self.out*(1. - self.out)
        return dx


class Affine:
    def __init__(self, W, b):
        self.W = W
        self.b = b
        self.x = None
        self.x_shape_org = None
        self.dW = None
        self.db = None

    def forward(self, x):
        self.x_shape_org = x.shape
        x = x.reshape(x.shape[0], -1)
        self.x = x

        out = np.dot(self.x, self.W) + self.b
        return out

    def backward(self, dout):
        dx = np.dot(dout, self.W.T)
        self.dW = np.dot(self.x.T, dout)
        self.db = np.sum(dout, axis=0)

        dx = dx.reshape(*self.x_shape_org)
        return dx


class Relu:
    def __init__(self):
        self.mask = None

    def forward(self, x):
        self.mask = x <= 0
        out = x.copy()
        out[self.mask] = 0
        return out

    def backward(self, dout):
        dout[self.mask] = 0
        dx = dout
        return dx


class Convolution:
    def __init__(self, W, b, S=1, P=0):
        self.W = W
        self.b = b
        self.S = S
        self.P = P

        self.x = None
        self.col = None
        self.col_W = None

        self.dW = None
        self.db = None

    def forward(self, x):
        FN, C, FH, FW = self.W.shape
        N, C, H, W = x.shape
        OH = int((H + 2*self.P - FH) / self.S + 1)
        OW = int((W + 2*self.P - FW) / self.S + 1)

        col = im2col(x, FH, FW, self.S, self.P)
        col_W = self.W.reshape(FN, -1).T

        out = np.dot(col, col_W) + self.b
        out = out.reshape(N, OH, OW, -1).transpose(0, 3, 1, 2)

        self.x = x
        self.col = col
        self.col_W = col_W

        return out

    def backward(self, dout):
        FN, C, FH, FW = self.W.shape
        dout = dout.transpose(0, 2, 3, 1).reshape(-1, FN)

        self.db = np.sum(dout, axis=0)
        self.dW = np.dot(self.col.T, dout)
        self.dW = self.dW.transpose(1, 0).reshape(FN, C, FH, FW)

        dcol = np.dot(dout, self.col_W.T)
        dx = col2im(dcol, self.x.shape, FH, FW, self.S, self.P)

        return dx


class ConvolutionK:
    def __init__(self, W, b, Sh=1, Sw=1, Ph=0, Pw=0):
        self.W = W
        self.b = b
        self.Sh = Sh
        self.Sw = Sw
        self.Ph = Ph
        self.Pw = Pw

        self.x = None
        self.col = None
        self.col_W = None

        self.dW = None
        self.db = None

    def forward(self, x):
        FN, C, FH, FW = self.W.shape
        N, C, H, W = x.shape
        OH = int((H + 2*self.Ph - FH) / self.Sh + 1)
        OW = int((W + 2*self.Pw - FW) / self.Sw + 1)

        col = im2colK(x, FH, FW, self.Sh, self.Sw,  self.Ph, self.Pw)
        col_W = self.W.reshape(FN, -1).T

        out = np.dot(col, col_W) + self.b
        out = out.reshape(N, OH, OW, -1).transpose(0, 3, 1, 2)

        self.x = x
        self.col = col
        self.col_W = col_W

        return out

    def backward(self, dout):
        FN, C, FH, FW = self.W.shape
        dout = dout.transpose(0, 2, 3, 1).reshape(-1, FN)

        self.db = np.sum(dout, axis=0)
        self.dW = np.dot(self.col.T, dout)
        self.dW = self.dW.transpose(1, 0).reshape(FN, C, FH, FW)

        dcol = np.dot(dout, self.col_W.T)
        dx = col2imK(dcol, self.x.shape, FH, FW, self.Sh, self.Sw, self.Ph, self.Pw)

        return dx


class Pooling:
    def __init__(self, PH, PW, S=2, P=0):
        self.PH = PH
        self.PW = PW
        self.S = S
        self.P = P
        self.x = None
        self.arg_max = None

    def forward(self, x):
        N, C, H, W = x.shape
        OH = int((H + 2*self.P - self.PH) / self.S + 1)
        OW = int((W + 2*self.P - self.PW) / self.S + 1)

        col = im2col(x, self.PH, self.PW, self.S, self.P)
        # reshaped to a 2D matrix with each row has a size  of PH*PW, from the
        # im2coled 2D matrix  with each row a C*PH*PW size.
        col = col.reshape(-1, self.PH*self.PW)
        arg_max = np.argmax(col, axis=1)
        out = np.max(col, axis=1)
        out = out.reshape(N, OH, OW, C).transpose(0, 3, 1, 2)

        self.x = x
        self.arg_max = arg_max

        return out

    def backward(self, dout):
        dout = dout.transpose(0, 2, 3, 1)     # N, OH, OW, C
        P_size = self.PH*self.PW
        dmax = np.zeros((dout.size, P_size))  # N*OH*OW*C, PH*PW
        dmax[np.arange(self.arg_max.size), self.arg_max.flatten()] \
            = dout.flatten()                  # N*OH*OW*C, PH*PW
        dmax = dmax.reshape((dout.shape + (P_size,)))  # N, OH, OW, C, P
        dcol = dmax.reshape(dmax.shape[0] * dmax.shape[1] * dmax.shape[2], -1)
                                              # N*OH*OW, C*PH*PW,
                                              # for input shape of col2im
        dx = col2im(dcol, self.x.shape, self.PH, self.PW, self.S, self.P)
                                              #  N, C, H, W

        return dx


class PoolingK:
    def __init__(self, PH, PW, Sh=2, Sw=2, Ph=0, Pw=0):
        self.PH = PH
        self.PW = PW
        self.Sh = Sh
        self.Sw = Sw
        self.Ph = Ph
        self.Pw = Pw
        self.x = None
        self.arg_max = None

    def forward(self, x):
        N, C, H, W = x.shape
        OH = int((H + 2*self.Ph - self.PH) / self.Sh + 1)
        OW = int((W + 2*self.Pw - self.PW) / self.Sw + 1)

        col = im2colK(x, self.PH, self.PW, self.Sh, self.Sw, self.Ph, self.Pw)
        col = col.reshape(-1, self.PH*self.PW)
        arg_max = np.argmax(col, axis=1)
        out = np.max(col, axis=1)
        out = out.reshape(N, OH, OW, C).transpose(0, 3, 1, 2)

        self.x = x
        self.arg_max = arg_max

        return out

    def backward(self, dout):
        dout = dout.transpose(0, 2, 3, 1)     # N, OH, OW, C
        P_size = self.PH*self.PW
        dmax = np.zeros((dout.size, P_size), dtype='float32')  # N*OH*OW*C, PH*PW
        dmax[np.arange(self.arg_max.size), self.arg_max.flatten()] \
            = dout.flatten()                  # N*OH*OW*C, PH*PW
        dmax = dmax.reshape((dout.shape + (P_size,)))  # N, OH, OW, C, P
        dcol = dmax.reshape(dmax.shape[0] * dmax.shape[1] * dmax.shape[2], -1)
        dx = col2imK(dcol, self.x.shape, self.PH, self.PW,
                    self.Sh, self.Sw, self.Ph, self.Pw)

        return dx


class BatchNormalization:
    def __init__(self, gamma, beta, momentum=0.9,
                 running_mean=None, running_var=None):
        self.gamma = gamma
        self.beta = beta
        self.momentum = momentum
        self.input_shape = None

        self.running_mean = running_mean
        self.running_var = running_var

        self.batch_size = None
        self.xp = None
        self.std = None
        self.dgamma = None
        self.dbeta = None

    def forward(self, x, train_flg=True):
        self.input_shape = x.shape
        if len(self.input_shape) != 2:
            N, C, H, W = x.shape
            # including H and W in the averaging
            x = x.transpose(0, 2, 3, 1).reshape(N*H*W, -1)
        out = self.__forward(x, train_flg)
        if len(self.input_shape) != 2:
            out = out.reshape(N, H, W, C).transpose(0, 3, 1, 2)
        return out

        #### original example
        #

        #if x.ndim != 2:
        #    N, C, H, W = x.shape
        #    x = x.reshape(N, -1)

        #out = self.__forward(x, train_flg)

        #return out.reshape(*self.input_shape)

    def __forward(self, x, train_flg):
        if self.running_mean is None:
            N, D = x.shape
            self.running_mean = np.zeros(D)
            self.running_var = np.zeros(D)

        if train_flg:
            mu = x.mean(axis=0)
            xp = x - mu
            var = np.mean(xp**2, axis=0)
            std = np.sqrt(var + 10e-7)
            xhat = xp/std

            self.batch_size = x.shape[0]
            self.xp = xp
            self.xhat = xhat
            self.std = std
            self.running_mean = self.momentum * self.running_mean \
                + (1-self.momentum) * mu
            self.running_var = self.momentum * self.running_var \
                + (1-self.momentum) * var

        else:
            xp = x - self.running_mean
            xhat = xp / ((np.sqrt(self.running_var + 10e-7)))

        out = self.gamma * xhat + self.beta
        return out

    def backward(self, dout):
        if len(self.input_shape) != 2:
            N, C, H, W = dout.shape
            dout = dout.transpose(0, 2, 3, 1).reshape(N*H*W, -1)

        dx = self.__backward(dout)
        if len(self.input_shape) != 2:
            dx = dx.reshape(N, H, W, C).transpose(0, 3, 1, 2)

        return dx

        ## original example 
        #if dout.ndim != 2:
        #    N, C, H, W = dout.shape
        #    dout = dout.reshape(N, -1)

        #dx = self.__backward(dout)

        #dx = dx.reshape(*self.input_shape)
        #return dx

    def __backward(self, dout):
        dbeta = dout.sum(axis=0)
        dgamma = np.sum(self.xhat*dout, axis=0)
        dxhat = self.gamma*dout
        dxp = dxhat / self.std
        dstd = -np.sum(dxhat/(self.std*self.std), axis=0)
        dvar = 0.5 * 1 / self.std * dstd
        dxp += (2.0 / self.batch_size) * self.xp * dvar
        dmu = -np.sum(dxp, axis=0)
        dx = dxp + dmu / self.batch_size

        self.dgamma = dgamma
        self.dbeta = dbeta

        return dx


class Dropout:
    def __init__(self, dropout_ratio=0.5):
        self.dropout_ratio = dropout_ratio
        self.mask = None

    def forward(self, x, train_flg=True):
        if train_flg:
            self.mask = np.random.rand(*x.shape) > self.dropout_ratio
            return x * self.mask
        else:
            return x * (1.0 - self.dropout_ratio)

    def backward(self, dout):
        return dout * self.mask


