#!/usr/bin/env python
# A test program for Gaussian Process.
# This is an exercise 2.9 #1 of Rasmussen's text book.
# Kazuyoshi TATSUMI 2021/03/22
import numpy as np
import math
import matplotlib.pyplot as plt


class GaussianProcess:
    def __init__(self, num_inputs=80, min_x=-5, max_x=5, ms=".", ls="", c="k"):
        self.num_inputs = num_inputs
        self.min_x = min_x
        self.max_x = max_x
        self.ms = ms
        self.ls = ls
        self.c = c
        self.generate_input()
        self.create_fig()

    def create_fig(self, title="Gaussian Process exercise"):
        self.fig = plt.figure(figsize=(6, 4.5))
        self.fig.suptitle(title)

    def generate_input(self):
        self.x = np.linspace(self.min_x, self.max_x, self.num_inputs)

    def generate_train(self):
        self.x_train = np.array([-2.3, 1.0, 3.5, -1.0, -4.0])
        self.y_train = np.array([1.11, 3.00,  -2.00, 4.0, 1.0])
        #self.x_train = np.array([-2.3, 1.0, 3.5, -1.0])
        #self.y_train = np.array([1.11, 3.00,  -2.00, 4.0])
        #self.x_train = np.array([-2.3, 1.0, 3.5])
        #self.y_train = np.array([1.11, 3.00,  -2.00])
        #self.x_train = np.array([-2.3, 1.0])
        #self.y_train = np.array([1.11, 3.00])
        #self.x_train = np.array([-2.3])
        #self.y_train = np.array([1.11])

    def calc_alpha(self):
        self.alpha = np.linalg.lstsq(self.L.T,
                                     np.linalg.lstsq(self.L, self.y_train)[0]
                                     )[0]

    def calc_fx_and_V_fx(self):
        k = np.zeros((self.x_train.shape[0]))
        self.fx = np.zeros((self.x.shape[0]))
        self.V_fx = np.zeros((self.x.shape[0]))
        for q, xq in enumerate(self.x):
            for p, xp in enumerate(self.x_train):
                k[p] = math.exp(-0.5*abs(xq - xp)**2)
            self.fx[q] = np.dot(k.T, self.alpha)
            v = np.linalg.lstsq(self.L, k)[0]
            self.V_fx[q] = 1.0 - np.dot(v.T, v)

    def calc_K(self, x1, x2):
        K = np.zeros((x1.shape[0], x2.shape[0]))
        for p, xp in enumerate(x1):
            for q, xq in enumerate(x2):
                K[p, q] = math.exp(-0.5*abs(xp - xq)**2)
        return(K)

    def calc_L(self, x):
        K = self.calc_K(x, x)
        K += 0.0000000000001*np.identity(x.shape[0])
        self.L = np.linalg.cholesky(K)

    def calc_covariance(self):
        K = self.calc_K(self.x, self.x)
        K_train = self.calc_K(self.x_train, self.x_train)
        K_offd = self.calc_K(self.x_train, self.x)
        self.cov = K - np.dot(K_offd.T, np.dot(np.linalg.inv(K_train), K_offd))
        self.cov += 0.0000000000001*np.identity(self.cov.shape[0])
        self.Lcov = np.linalg.cholesky(self.cov)

    def calc_covariance_using_L(self):
        K_train = self.calc_K(self.x_train, self.x_train)
        L_train = np.linalg.cholesky(K_train +
                                     0.0000000000001*np.identity(K_train.shape[0]))
        K = self.calc_K(self.x, self.x)
        K_offd = self.calc_K(self.x_train, self.x)
        Km1K = np.zeros(K_offd.shape)
        for p in range(0, K_offd.shape[1]):
            Km1K[:, p] = np.linalg.lstsq(L_train.T,
                                         np.linalg.lstsq(L_train,
                                                         K_offd[:, p])[0]
                                         )[0]
        self.cov = K - np.dot(K_offd.T, Km1K)
        self.cov += 0.0000000000001*np.identity(self.cov.shape[0])
        self.Lcov = np.linalg.cholesky(self.cov)

    def draw_random_func(self):
        u = np.random.normal(size=self.num_inputs)
        f = np.dot(self.L, u)
        plt.plot(self.x, f)

    def draw_random_func_with_train(self):
        u = np.random.normal(size=self.Lcov.shape[0])
        f = np.dot(self.Lcov, u)
        plt.plot(self.x, f+self.fx)

    def output_with_train(self):
        self.sd = np.sqrt(self.V_fx)
        plt.fill_between(self.x, self.fx + self.sd*2.0,
                         self.fx - self.sd*2.0, fc='lightgray')
        plt.plot(self.x, self.fx, c="r")
        plt.plot(self.x_train, self.y_train, marker="+", linestyle="", c='k')


def samplerun():
    prj = GaussianProcess()
    prj.calc_L(prj.x)
    for i in range(3):
        prj.draw_random_func()
    plt.xlabel("input, x")
    plt.ylabel("output, f(x)")


def samplerun_test():
    prj = GaussianProcess()
    prj.generate_train()
    prj.calc_covariance_using_L()
    prj.calc_L(prj.x_train)
    prj.calc_alpha()
    prj.calc_fx_and_V_fx()
    prj.output_with_train()
    plt.xlabel("input, x")
    plt.ylabel("output, f(x)")
    for i in range(0):
        prj.draw_random_func_with_train()


samplerun_test()
plt.show()



