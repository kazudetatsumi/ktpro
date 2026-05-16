import numpy as np


def numerical_gradient(f, x):
    h = 1e-4
    grad = np.zeros_like(x)

    for xidx in range(x.size):
        tmp_val = x[xidx]
        x[xidx] = tmp_val + h
        f_fwd = f(x)
        x[xidx] = tmp_val - h
        f_bwd = f(x)

        grad = (f_fwd - f_bwd) / (2.*h)
        x[xidx] = tmp_val

    return grad
