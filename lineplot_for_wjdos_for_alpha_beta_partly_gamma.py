#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
#import scipy.optimize
plt.style.use('classic')
plt.rcParams['font.family'] = 'Times New Roman'


homedir = "/home/kazu/"
cdataA0 = np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g0-A-t300.dat',    comments='#', dtype='float')
sdataA0 = np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g0-A-t300.dat', comments='#', dtype='float')
gdataA0 = np.loadtxt(homedir + '/gamma-si3n4-unit/phono3py_111_fc2_222_sym_monk_k-shift/jdos-m222222-g0-t300.dat', comments='#', dtype='float')
cdataA1 = np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g784-05A-t300.dat', comments='#', dtype='float')
sdataA1 = np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g1568-05A-t300.dat', comments='#', dtype='float')
cdataA2 = np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g1568-A-t300.dat', comments='#', dtype='float')
sdataA2 = np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g3136-A-t300.dat', comments='#', dtype='float')
cdataK0 = np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g0-K-t300.dat', comments='#', dtype='float')
sdataK0 = np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g0-K-t300.dat', comments='#', dtype='float')
cdataK1 = np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g56-05K-t300.dat', comments='#', dtype='float')
sdataK1 = np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g56-05K-t300.dat', comments='#', dtype='float')
cdataK2 = np.loadtxt(homedir + '/asi3n4/phono3py_112_fc2_334_sym_monk_shift/jdos-m141416-g98-K-t300.dat', comments='#', dtype='float')
sdataK2 = np.loadtxt(homedir + '/bsi3n4_m/phono3py_113_fc2_338_sym_monk_shift/jdos-m141432-g98-K-t300.dat', comments='#', dtype='float')
cz = 4
sz = 2
gz = 2
c3na = 3*cz*7
s3na = 3*sz*7
g3na = 3*gz*7


plt.figure(figsize=(4.5, 6))
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.subplot(6, 1, 1)
plt.subplots_adjust(hspace=0, wspace=0)
plt.plot(cdataA0[:, 0], (cdataA0[:, 1]+cdataA0[:, 2])/c3na**2, color="k", label="aSi3N4", linewidth=1.3, linestyle=(0, (1.0, 1.5)))
plt.plot(sdataA0[:, 0], (sdataA0[:, 1]+sdataA0[:, 2])/s3na**2, color="k", label="bSi3N4", linewidth=0.8)
#plt.plot(gdataA0[:, 0], (gdataA0[:, 1]+gdataA0[:, 2])/gz**2, color="gray", label="gSi3N4", linewidth=0.1)
yz = 0*gdataA0[:, 0]
plt.fill_between(gdataA0[:, 0], yz, (gdataA0[:, 1]+gdataA0[:, 2])/g3na**2, facecolor='lightgray', linewidth=0.1)
plt.text(52, 19, 'q=(0, 0, 0)')
x = range(0, 70, 10)
plt.xticks(x, " ")
plt.ylim(0, 0.064)
plt.yticks([0, 0.03, 0.06])


plt.xlabel("Frequency (THz)")
plt.subplot(6, 1, 2)
plt.plot(cdataA1[:, 0], (cdataA1[:, 1]+cdataA1[:, 2])/c3na**2, color="k", label="aSi3N4", linewidth=1.3, linestyle=(0, (1.0, 1.5)))
plt.plot(sdataA1[:, 0], (sdataA1[:, 1]+sdataA1[:, 2])/s3na**2, color="k", label="bSi3N4", linewidth=0.8)
plt.text(52, 19, 'q=(0, 0, 1/4)')
plt.xticks(x, " ")
plt.ylim(0, 0.064)
plt.yticks([0, 0.03, 0.06])
#plt.xlabel("Frequency (THz)")
#plt.ylabel("JDOS [1/THz]")

plt.subplot(6, 1, 3)
plt.plot(cdataA2[:, 0], (cdataA2[:, 1]+cdataA2[:, 2])/c3na**2, color="k", label="aSi3N4", linewidth=1.3, linestyle=(0, (1.0, 1.5)))
plt.plot(sdataA2[:, 0], (sdataA2[:, 1]+sdataA2[:, 2])/s3na**2, color="k", label="bSi3N4", linewidth=0.8)
plt.text(52, 19, 'q=(0, 0, 1/2)')
plt.text(1, 21, r'$\alpha$ Si$_3$N$_4$')
plt.text(1, 14, r'$\beta$ Si$_3$N$_4$')
plt.text(1, 9, r'$\gamma$ Si$_3$N$_4$')
plt.xticks(x, " ")
plt.ylim(0, 0.064)
plt.yticks([0, 0.03, 0.06])
#plt.xlabel("Frequency (THz)")
#plt.ylabel("JDOS [1/THz]")

plt.subplot(6, 1, 4)
plt.plot(cdataK0[:, 0], (cdataK0[:, 1]+cdataK0[:, 2])/c3na**2, color="k", label="aSi3N4", linewidth=1.3, linestyle=(0, (1.0, 1.5)))
plt.plot(sdataK0[:, 0], (sdataK0[:, 1]+sdataK0[:, 2])/s3na**2, color="k", label="bSi3N4", linewidth=0.8)
plt.text(52, 19, 'q=(0, 0, 0)')
plt.xticks(x, " ")
plt.ylim(0, 0.064)
plt.yticks([0, 0.03, 0.06])
#plt.ylim(0, 28)
#plt.yticks([0, 10, 20])
#plt.xlabel("Frequency (THz)")
#plt.ylabel("JDOS [1/THz]")

plt.subplot(6, 1, 5)
plt.plot(cdataK1[:, 0], (cdataK1[:, 1]+cdataK1[:, 2])/c3na**2, color="k", label="aSi3N4", linewidth=1.3, linestyle=(0, (1.0, 1.5)))
plt.plot(sdataK1[:, 0], (sdataK1[:, 1]+sdataK1[:, 2])/s3na**2, color="k", label="bSi3N4", linewidth=0.8)
plt.text(52, 19, 'q=(0, 1/4, 0)')
plt.xticks(x, " ")
#plt.ylim(0, 28)
#plt.yticks([0, 10, 20])
plt.ylim(0, 0.064)
plt.yticks([0, 0.03, 0.06])
#plt.xlabel("Frequency (THz)")
#plt.ylabel("JDOS [1/THz]")

plt.subplot(6, 1, 6)
plt.plot(cdataK2[:, 0], (cdataK2[:, 1]+cdataK2[:, 2])/c3na**2, color="k", label="aSi3N4", linewidth=1.3, linestyle=(0, (1.0, 1.5)))
plt.plot(sdataK2[:, 0], (sdataK2[:, 1]+sdataK2[:, 2])/s3na**2, color="k", label="bSi3N4", linewidth=0.8)
plt.text(52, 19, 'q=(0, 1/2, 0)')
#plt.ylim(0, 28)
#plt.yticks([0, 10, 20])
plt.ylim(0, 0.064)
plt.yticks([0, 0.03, 0.06])
plt.xlabel("Frequency (THz)")
plt.ylabel(r'WJDOS (THz$^{-1}$f.u.$^{-2}$)')

plt.savefig("wjdos-aSi3N4_bSi3N4_joint.pdf")
plt.show()
