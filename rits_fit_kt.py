import sys
import os
import numpy as np

np.random.seed(222)

def get_sim_spectrum(inpfile='rits_initial.inp'):
    sys.path.append("/home/kazu/desktop/240424/rits/src/cpp/rietveldcc/Linux")
    from AdvRits import AdvRietveld
    rits = AdvRietveld()
    if not os.path.exists("spgra"):
        print("'spgra' dose not exist.")
        return

    rits.SetInitialFileName(inpfile)
    rits.ReadInitialFile()
    rits.SetDataFileName('fit_results.dat.tmp')
    rits.ReadDataFile()
    rits.SetFittingResultDataFileName('fit_func.data')
    rits.SetFittingResultParamFileName('fit_para.data')

    rits.Fit()

    return np.array(rits.PutObservedLambda()),\
        np.array(rits.PutFittedIntensity())


def get_sim_edgespectrum(inpfile='edge_3.inp'):
    sys.path.append("/home/kazu/desktop/240424/rits/src/cpp/edgecc/Linux")
    from AdvEdge import AdvEdge
    rits = AdvEdge()
    rits.SetInitialFileName(inpfile)
    rits.ReadInitialFile()
    rits.SetDataFileName('fit_results.dat.tmp')
    rits.ReadDataFile()
    rits.SetFittingResultDataFileName('fit_func.data')
    rits.SetFittingResultParamFileName('fit_para.data')

    rits.Fit()

    return np.array(rits.PutObservedLambda()),\
        np.array(rits.PutFittedIntensity())

#x, y = get_sim_spectrum()
#x = np.array(x)
#y = np.array(y)

#timescale = 3600000
#y *= timescale
#sensitivity = 0.01
#import matplotlib.pyplot as plt
#_y = np.random.poisson(y*sensitivity*x)/(sensitivity*x)
#plt.scatter(x, _y, s=15, fc='none', ec='b')
#plt.plot(x, y, c='r')
#plt.plot(x, y - _y) 
#plt.xlim([1, 5.1])
#plt.show()
