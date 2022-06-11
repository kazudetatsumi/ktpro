#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss


def run():
    #mean = np.zeros((4, 2))
    #cov = np.zeros((4, 25, 2))
    mean = np.array([50, 50])
    cov = np.array([[25,0], [0, 20]])

    x = np.linspace(0.0, 100.0, 300)
    y = np.linspace(0.0, 100.0, 300)
    X, Y = np.meshgrid(x, y)
    pos = np.dstack((X, Y))

    mean = np.array([50, 60])
    cov = np.array([[50,0], [0, 80]])
    pdf = ss.multivariate_normal(mean, cov).pdf(pos)
    mean = np.array([70, 60])
    cov = np.array([[1225,0], [0, 10]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)
    mean = np.array([0, 60])
    cov = np.array([[525,0], [0, 5]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)
    mean = np.array([0, 60])
    cov = np.array([[425,0], [0, 375]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)*3.
    mean = np.array([75, 60])
    cov = np.array([[50,0], [0, 55]])*10
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)*10.
    mean = np.array([135, 60])
    cov = np.array([[50,0], [0, 55]])*10
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)*10.
    mean = np.array([80, 60])
    cov = np.array([[5,0], [0, 1]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)*0.01
    mean = np.array([90, 60])
    cov = np.array([[5,0], [0, 1]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)*0.01
    mean = np.array([55, 60])
    cov = np.array([[25,3], [-10, 4800]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)
    mean = np.array([75, 60])
    cov = np.array([[25,3], [-10, 4800]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)
    mean = np.array([85, 60])
    cov = np.array([[25,3], [-10, 4800]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)
    mean = np.array([100, 60])
    cov = np.array([[25,3], [-10, 4800]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)
    mean = np.array([15, 60])
    cov = np.array([[25,3], [-10, 4800]])
    pdf = pdf + ss.multivariate_normal(mean, cov).pdf(pos)

    fig = plt.figure(figsize=(6,16))
    ax = fig.add_subplot(211)
    plt.pcolor(X, Y, pdf, cmap='gray', vmin = -1. *np.max(pdf)*0.5)
    ax = fig.add_subplot(212)
    counts =  np.random.poisson(pdf*50000)
    plt.pcolor(X, Y, counts, vmin = -1.*np.max(counts)*0.5, cmap='gray')
    plt.show()



run()
