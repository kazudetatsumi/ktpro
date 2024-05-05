import sys
import os
import numpy as np
sys.path.append("/home/kazu/desktop/240424/rits/src/cpp/rietveldcc/Linux")
from AdvRits import AdvRietveld

os.environ['LD_LIBRARY_PATH'] = '/home/kazu/desktop/240424/rits/src/cpp/rietveldcc/Linux'
sys.path.append("/home/kazu/desktop/240424/rits/src/cpp/rietveldcc/Linux")
np.random.seed(222)

def get_sim_spectrum():
    sys.stdout = open(os.devnull, 'w')
    rits = AdvRietveld()
    if not os.path.exists("spgra"):
        print("'spgra' dose not exist.")
        return

    rits.SetInitialFileName('rits_initial.inp')
    rits.ReadInitialFile()
    rits.SetDataFileName('fit_results.dat')
    rits.ReadDataFile()
    rits.SetFittingResultDataFileName('fit_func.data')
    rits.SetFittingResultParamFileName('fit_para.data')

    rits.Fit()

    return np.array(rits.PutObservedLambda()), np.array(rits.PutFittedIntensity())

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
