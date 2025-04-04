#!/usr/bin/env python
import numpy as np
import math
import matplotlib.pyplot as plt
import sys
from scipy import stats
sys.path.append("/home/kazu/desktop/210108/AdaptiveKDE/adaptivekde")
import ssvkernel


def gaussian_func(x, y, sigma):
    exponent = - (x ** 2 + y ** 2) / (2 * pow(sigma, 2))
    denominator = 2 * math.pi * pow(sigma, 2)
    return pow(math.e, exponent) / denominator

def double_gaussian_func(x, y, sigma):
    exponent1 = - (x ** 2 + y ** 2) / (2 * pow(sigma, 2))
    denominator1 = 2 * math.pi * pow(sigma, 2)
    exponent2 = - ((x - 5*sigma) ** 2 + y ** 2) / (2 * pow(sigma, 2))
    denominator2 = 2 * math.pi * pow(sigma, 2)
    return pow(math.e, exponent1) / denominator1 + pow(math.e, exponent2) / denominator2

def single_gaussian_func(x, sigma):
    exponent = - (x ** 2) / (2 * pow(sigma, 2))
    denominator = (2 * math.pi * pow(sigma, 2))**0.5
    return pow(math.e, exponent) / denominator
np.random.seed(314)
M = 80
winparam = 5
WinFunc = "Boxcar"

dr = 0.01
x = np.arange(-3.0, 3.0, dr)
y = np.arange(-3.0, 3.0, dr)
X, Y = np.meshgrid(x, y)
Z = double_gaussian_func(X, Y, 0.1)
z = single_gaussian_func(x, 0.1)
#print('fuck', z.sum()*dr)
#plt.pcolor(X, Y, Z)
#plt.show()

#tmpdata = np.random.randn(1000,2)
#plt.scatter(tmpdata[:,0], tmpdata[:, 1])
#plt.show()
data = np.random.poisson(Z*dr**2*20000)

vec1 = []
vec2 = []
for ix in range(data.shape[0]):
    for iy in range(data.shape[1]):
        for irpt in range(data[ix, iy]):
            vec1.append(x[ix])
            vec2.append(y[iy])

vec1 = np.array(vec1)
vec2 = np.array(vec2)
#data1 = np.sum(data, axis=1)
#data2 = np.sum(data, axis=0)
#plt.plot(data1)
#plt.show()

#vec1 = np.array([idx for idx in range(0, x.shape[0]) for irpt in
#                range(0, int(data1[idx]))], dtype=float)
#vec2 = np.array([idx for idx in range(0, y.shape[0]) for irpt in
#                range(0, int(data2[idx]))], dtype=float)

#vec1_real = np.array([x[idx] for idx in range(0, x.shape[0]) for irpt in
#                     range(0, int(data1[idx]))], dtype=float)
#vec2_real = np.array([y[idx] for idx in range(0, y.shape[0]) for irpt in
#                     range(0, int(data2[idx]))], dtype=float)

shift1 = np.random.uniform(-0.5, 0.5, size=vec1.shape[0])
shift2 = np.random.uniform(-0.5, 0.5, size=vec2.shape[0])
vec1 += shift1*dr
vec2 += shift2*dr
print(vec1.shape)
#vec1_real += shift1*dr
#vec2_real += shift2*dr
#print(vec1_real)

#tmpdata  = np.random.randn(10000,2)
f1 = ssvkernel.ssvkernel(vec1, x, M=M, winparam=winparam, WinFunc=WinFunc)
f2 = ssvkernel.ssvkernel(vec2, y, M=M, winparam=winparam, WinFunc=WinFunc)
#f1 = ssvkernel.ssvkernel(tmpdata[:, 0], M=M, winparam=winparam, WinFunc=WinFunc)
#f2 = ssvkernel.ssvkernel(tmpdata[:, 1], M=M, winparam=winparam, WinFunc=WinFunc)
print('chk1', x[0:10])
print('chk1', f1[1][0:10])



#print(np.sum(f1[0])*dr)
plt.plot(f1[1], f1[0])
plt.plot(f2[1], f2[0])
plt.plot(x, z)
#plt.plot(x, np.sum(Z, axis=1)*dr)
#plt.xlim(-0.3,0.3)
plt.show()
#plt.plot(f2[1], f2[0])
#plt.plot(y, np.sum(Z, axis=0)*dr)
#plt.show()

u1 = []
u2 = []
for v1, v2 in zip(vec1, vec2):
    u1.append(f1[0][f1[1] <= v1].sum())
    u2.append(f2[0][f2[1] <= v2].sum())
u1 = np.array(u1)*dr
u2 = np.array(u2)*dr
u = np.vstack([u1, u2])
fig, ax = plt.subplots(3, 1)
ax[0].scatter(u1, u2, s=0.1)
ax[1].hist2d(u1, u2, bins=15)
UX, UY = np.mgrid[np.min(u1):np.max(u1):100j, np.min(u2):np.max(u2):100j]
copulad = stats.gaussian_kde(u)
positions = np.vstack([UX.ravel(), UY.ravel()])
#print(kde(np.vstack([UX.ravel(), UY.ravel()])).reshape(UX.shape).shape)
c = ax[2].pcolor(UX, UY, copulad(positions).reshape(UX.shape))
plt.colorbar(c)
plt.show()


u1x1 = []
u2x2 = []
for idx, (x1, x2) in enumerate(zip(X.ravel(), Y.ravel())):
    u1x1.append(f1[0][f1[1] <= x1].sum())
    u2x2.append(f2[0][f2[1] <= x2].sum())
u1x1 = np.array(u1x1)*dr
u2x2 = np.array(u2x2)*dr
positions2 = np.vstack([u1x1, u2x2])
f1f2 = np.outer(f1[0], f2[0])
f_cop = f1f2*copulad(positions2).reshape(X.shape)

fig, ax = plt.subplots(2, 1)
ax[0].pcolor(X, Y, f_cop)
ax[0].set_title('Estimated using couplad')
ax[1].pcolor(X, Y, Z)
ax[1].set_title('Ground truth')
plt.legend()
plt.show()



#kde1 = ssvkernel.ssvkernel(vec1, x, M=M, winparam=winparam, WinFunc=WinFunc)
#kde2 = ssvkernel.ssvkernel(vec2, y, M=M, winparam=winparam, WinFunc=WinFunc)

#    proj = qc.qens(datadir, save_file, qsel=True, winparam=winparam, M=M,
#                   WinFunc=WinFunc, figname=figname, showplot=False)
#    proj.select_spectra()
#    proj.add_shift()
#    proj.run_ssvkernel()


#    def select_spectra(self):
#        spectra = self.dataset['spectra']
#        if self.odata:  # case outgoing beam
#            if self.qsel:   # case spectra are already integrated over a
#                            # specific q range.
#                self.selected_spectra = self.dataset['spectra']
#                self.selected_energy = self.dataset['energy']
#            else:       # case dataset are distributed over 2-D PSD elements.
#                dp = self.dataset['detector_position']
#                mask = np.where((dp[:, 0] >= 10) & (dp[:, 1] <= 65))[0]
#                self.selected_spectra = np.sum(spectra[mask, 1, :], axis=0)
#                self.selected_energy = spectra[0, 0, :]
#        else:
#            spectra[0, 0, :] = spectra[0, 0, :] - 2.085
#            mergin = 0.001
#            xlim = np.array([-0.10025 - mergin, 0.14975 + mergin])
#            mask = np.where((spectra[0, 0, :] >= xlim[0]) &
#                            (spectra[0, 0, :] <= xlim[1]))[0]
#            self.selected_spectra = spectra[0, 1, mask]
#            self.selected_energy = spectra[0, 0, mask]
#        #print(self.selected_energy[np.argmax(self.selected_spectra)])
#        self.de = self.selected_energy[1] - self.selected_energy[0]
#        print("check!!!", self.selected_energy[0], self.selected_energy[-1],self.de)
#        self.get_xvec()

#    def get_xvec(self):
#        # To keep the accuracy, kde is executed on the channel numbers in the
#        # histogram data.
#        # The actual energies were retrieved by "_real" variables.
#        print("CHK", self.selected_spectra.shape[0])
#        self.xvec = np.array([idx for idx in
#                             range(0, self.selected_spectra.shape[0]) for
#                             num_repeat in
#                             range(0, int(self.selected_spectra[idx]))
#                              ], dtype=float)
#        self.xvec_real = np.array([self.selected_energy[idx] for idx in
#                                  range(0, self.selected_spectra.shape[0]) for
#                                  num_repeat in
#                                  range(0, int(self.selected_spectra[idx]))
#                                   ], dtype=float)
#
#    def add_shift(self):
#        #np.random.seed(314)
#        self.xvecorg = np.array(self.xvec)
#        self.shift = np.random.uniform(-0.5, 0.5, size=self.xvec.shape[0])
#        self.xvec += self.shift
#        self.xvec_real += self.shift*self.de
#        #print(self.xvec[0:30])

#    # Since get_qlist_nova_class.get_alldata() takes the bin bottom energy as the neutron energy,
#    # the random shift enegy in the following method is modified to be [0, 1]*de.
#    # However, the  q values are still inaccurate.
#    def add_shift_de(self):
#        #np.random.seed(314)
#        self.xvecorg = np.array(self.xvec)
#        self.shift = np.random.uniform(0., 1.0, size=self.xvec.shape[0])
#        self.xvec += self.shift
#        self.xvec_real += self.shift*self.de
#        #print(self.xvec[0:30])
#
#    def run_ssvkernel(self, num=800, isnovariablebw=False):
#        #self.optsm = True
#        if self.optsm:
#            tinmax = 10.**int(np.log10(self.selected_spectra.shape[0])+1.)
#            self.tin = np.linspace(0.0, tinmax, int(tinmax)*2+1)
#            print("Check parameters of horizontal axis")
#            print("de=", self.de, "selected_energy[0]=",
#                  self.selected_energy[0], "num channels=", self.tin.shape[0])
#            self.tin_real = self.tin*self.de + self.selected_energy[0]
#        else:
#            self.tin = np.arange(self.selected_energy.shape[0])
#            print("Check parameters of horizontal axis")
#            print("de=", self.de, "selected_energy[0]=",
#                  self.selected_energy[0], "num channels=", self.tin.shape[0])
#            self.tin_real = np.linspace(self.selected_energy[0],
#                                        self.selected_energy[-1],
#                                        num=self.selected_spectra.shape[0])
#                                        #num=800)
#                                        #num=8000)
#                                        #num=66670)
#                                        #num=200000)
#            print("num of tin_real element=",self.tin_real.shape)
#            print(self.tin_real[0:10])
#            print(self.selected_energy[0:20])
#            print(self.selected_spectra[0:20])
#        print('Total count number:', self.xvec_real.shape)
#        self.y = ssvkernel.ssvkernel(self.xvec_real, self.tin_real, M=self.M,
#                                     winparam=self.winparam,
#                                     WinFunc=self.WinFunc)
#        self.y_ = sskernel.sskernel(self.xvec_real, self.tin_real)
#
