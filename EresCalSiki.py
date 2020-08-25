#!/usr/bin/env python

#
# EresCalSiki.py
# 
# Calculate the energy resolution on SIKI
#
# R.K., Jan. 25, 2011, updated by Y. I.
#       Feb. 27, 2014 added the term due to the uncertinitiy of L2 (DL2)
#       May. 29, 2014 display style, DE/Ei or DE, can be chosen
#       Mar. 26, 2015 the last plot can be saved as a csv file
#       Aug. 28, 2017 use average path length for a cylinder as DL2
#       Feb. 19, 2018 added the sample size in the header of graph. Added a header when saved in a csv file
#       Nov. 17, 2018 became compatible with python 3. Axes are autoscaled. Added legend(s). Fixed some bugs.
#       Jan. 20, 2020 changed the slit width to 0.44 mm, which is 10% larger than the designed value.
#
# Usage:
#   Input values of Ei (meV), Frequency (Hz) and Sample size (mm).
#   No input value ends the script.
#
### for python 3 compatibility
from __future__ import print_function
if hasattr(__builtins__, 'raw_input'):
    input = raw_input
###
from scipy import *
from scipy.interpolate import interp1d
#import matplotlib as mpl
#mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'
#mpl.rcParams['axes.xmargin'] = 0
#mpl.rcParams['axes.ymargin'] = 0
import matplotlib.pyplot as plt
import csv

# Ei vs Dtm (data from http://j-parc.jp/MatLife/ja/source/data/Pulse_0310_mod.xls)
Eis = [
    0.099999, 0.12589, 0.15849, 0.19953, 0.25119, 0.31623, 0.3981, 0.50118,
    0.63096,  0.79433, 0.99999, 1.2589,  1.5849,  1.9953,  2.5119, 3.1623,
    3.981,    5.0118,  6.3096,  7.9433,  9.9999,  12.589,  15.849, 19.953,
    25.119,   31.623,  39.81,   50.118,  63.096,  79.433,  99.999, 125.89,
    158.49,   199.53,  251.19,  316.23,  398.1,   501.18,  630.96, 794.33,
    999.99,   1258.9,  1584.9,  1995.3,  2511.9,  3162.3,  3981,   5011.8,
    6309.6,   7943.3,  9999.9
    ]
Dtms = [
    487.05,   487.71,  461.9,   448,     417.73,  383.09,  363.51, 337.71,
    307.18,   286.35,  262.2,   241.77,  227.26,  208.89,  188.66, 172.44,
    155.93,   137.95,  119.63,  105.82,  89.559,  76.105,  58.736, 42.069,
    31.527,   25.151,  20.432,  17.571,  15.06,   13.605,  12.149, 10.844,
    9.6387,   8.5343,  7.5804,  6.7772,  5.974,   5.472,   4.8193, 4.4679,
    4.0161,   3.5141,  3.1627,  2.8113,  2.5603,  2.2591,  2.0583, 1.8073,
    1.6566,   1.4558,  1.3052
    ]
Dtm0 = interp1d(Eis, Dtms) # linear interpolation

# Instrumental parameters
L1 = 18.0 # moderator-sample distance [m]
L2 = 2.5  # sample-detector distance [m]
L3 = 1.71 # chopper-sample distance [m]
wdet = 0.019 # width of a detector tube [m]

Lgs = 2200.    # guide-sample distance [mm]
wguide = 43.35 # width of the end of the guide [mm]
mvalue = 4.    # m-value of the supermirror
gslope = (56.74-43.35)/2./1755.*1000
               # slope of the last section of the guide [mrad]
wslit = 0.44    # width of each slit of the chopper [mm]
lslit = 20.0  # length of each slit of the chopper [mm]
# Note: for the old chopper, wslit = 2.0 mm and lslit = 100 mm, but wslit/lslit is the same.
# Note2: wslit is 10% larger than the designed value 0.4 mm.

# angular divergence defined by the chopper [mrad]
div_chopper = wslit/lslit*1000.

#
# calculate and plot the resolution
#
def CalcEresolution(Ei,freq,ssize):
    # pulse width at the moderator [usec]
    Dtm = Dtm0(Ei)

    # time when neutrons reach the chopper [usec]
    tch = 2286./sqrt(Ei)*(L1-L3)

    # angular divergence of the incident beam [mrad]
    #  (1) divergence by reflection by the supermirror guide:
    #    slope of the last section of the guide
    #     + critical angle of the mirror
    div_guide = 2.*(gslope + 1.73*mvalue*9.044/sqrt(Ei))
    #  (2) divergence defined by the sample [mrad]
    div_sample = 2.*(wguide/2.+ssize/2.)/Lgs*1000.
    div_beam = min(div_guide,div_sample)

    # broadening of the pulse width at the chopper
    # due to the light-house effect [cf. Windsor's textbook]
    u = div_beam/div_chopper
    if u <= 0.8:
        p = 1. + u/4.
    else:
        if u <= 2.:
            p = 2. + u - sqrt(4.*u-u*u)
        else:
            p = u

    # opening time of the chopper [usec]
    Dtch0 = wslit/(2.*pi*lslit*freq)*1e6   # nominal value
    Dtch = Dtch0*p   # broadening due to the lighthouse effect

    # uncertinity of L2 due to the sample size and the detector size
#    DL2 = sqrt(wdet**2+(ssize/1000.)**2)
    DL2 = sqrt((wdet/4.0*3.1415927)**2+(ssize/1000./4.0*3.1415927)**2)

    # array of energy transfer
    e = arange(0, Ei+Ei/100., Ei/100.)

    # energy resolution relative to Ei [cf. Windsor's textbook]
    De_Ei =sqrt((2.*Dtch/tch*(1.+L1/L2*(1.-e/Ei)**1.5))**2
                + (2.*Dtm/tch*(1.+L3/L2*(1.-e/Ei)**1.5))**2
                + (2.*DL2/L2*(1.-e/Ei))**2)

    return e, De_Ei

#
# plotting
#
def PlotRes(e, Eres, labelstring):
#    plt.clf()
    plt.title('Energy resolution of 4SEASONS: \nEi = ' + str(Ei)
          + ' meV, f = ' + str(freq) + ' Hz, sample diameter = ' + str(ssize) +' mm')
    plt.xlabel('Energy transfer (meV)')
    plt.grid(True)
    plt.ylabel(labelstring)
    plt.plot(e, Eres, label="%d meV, %d Hz, %d mm" %(Ei, freq, ssize))
    plt.axis('auto')
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.legend()
    plt.draw()  ##[inamura]
#    plt.draw()  ##[inamura] # Only one draw() seems to be enough.

#
# main loop
#

while(True):
    disp_flag = input("Show (0) DE/Ei [%] or (1) DE [meV] ? = ")
    if disp_flag not in ["0", "1"]:
        print("### Input 0 or 1. ###")
        continue
    disp_flag = int(disp_flag)
    break

#plt.ion() ##[inamura] # ion() here causes very high CPU usage with matplotlib 2.1.1 & qt4.

while(True):
    # read parameters
    inputs = input("Ei (meV), Frequency (Hz), Sample size (mm) = ")
    if inputs == "":
        plt.ioff()  ##[inamura]
        break
    inputs = inputs.split(",")
    if len(inputs) < 3:
        print("### Too small number of inputs. ###")
        continue

    Ei = float(inputs[0])
    freq = float(inputs[1])
    ssize = float(inputs[2])

    if (Ei < 0.099999) or (Ei > 9999.9):
        print("### Ei is out of range. ###")
        continue
    if freq <= 0:
        print("### Frequency must be positive. ###")
        continue
    if ssize < 0:
        print("### Sample size must be positive or zero. ###")
        continue

    e, De_Ei = CalcEresolution(Ei,freq,ssize)
    if disp_flag == 1:
        labelstring = "Energy resolution (meV)"
        header = ["E(meV)", "DE(meV)"]
        Eres = De_Ei*Ei
    else:
        labelstring = "Energy resolution / Ei (%)"
        header = ["E(meV)", "DE/E(%)"]
        Eres = De_Ei*100.

    plt.ion()
    PlotRes(e, Eres, labelstring)
    plt.ioff()

#plt.show()  ##[inamura] # This show() is not required anymore to keep graph, because another process to save the result was added below.
    
# save the result in a csv file.
save = input("Save the last plot as text? [y/n] = ")
if save == "y":
    rows = array([e, Eres])
    with open('Eresolution.csv', 'w') as f:
        writer = csv.writer(f,lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows.T)
    print("Saved as ./Eresolution.csv")
        

