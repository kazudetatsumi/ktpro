#!?usr/bin/env python
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os

def get_tdata():
    tdatafile = 'tdata.pkl'
    if not os.path.exists(tdatafile):
        for i in range(0, 5):
            with open('./seed' + str(i) + '/valtesttot.pkl', 'rb') as f:
                data = pickle.load(f)
                var = pickle.load(f)
                if i == 0:
                    tdata = data.squeeze()[np.newaxis]
                    tvar = var.squeeze()[np.newaxis]
                else:
                    tdata = np.vstack((tdata, data.squeeze()[np.newaxis]))
                    tvar = np.vstack((tvar, var.squeeze()[np.newaxis]))
        with open('tdata.pkl', 'wb') as f:
            pickle.dump(tdata.astype('float32'), f, 4)
            pickle.dump(tvar.astype('float32'), f, 4)
    else:
        print('reading ' + tdatafile)
        with open(tdatafile, 'rb') as f:
            tdata = pickle.load(f)
            tvar = pickle.load(f)
    return tdata, tvar


def check_std():
    tdata, tvar = get_tdata()
    with open('/home/kazu/desktop/240424/uNID_data_KO/255/openbeamstrd.pkl',
              'rb') as f:
        openbeam_expt = pickle.load(f)
        openbeam_expt_noisy = pickle.load(f)
        sample_expt = pickle.load(f)
        sample_expt_noisy = pickle.load(f)
        etc = pickle.load(f)
    sample_expt_noisy = openbeam_expt_noisy.transpose((2, 1, 0))
    #std = np.std(tdata, axis=0)
    M = tdata.shape[0]
    std = (1/M*(tvar.sum(axis=0) + (tdata**2).sum(axis=0)) - np.mean(tdata, axis=0)**2)**0.5
    maxstd = np.max(std[:-2], axis=0)
    poslists = [[40, 130], [30, 50], [30, 30], [40, 15], [66, 79]]
    numsets = len(poslists)
    #print(np.argmax(maxstd[40, 130]))
    #diff = (tdata - np.average(tdata, axis=0))**2
    #print(diff.shape)
    #print(np.argmax(np.sum(diff[:, :, 40, 130, 82], axis=0)))
    #print(np.argmin(np.sum(diff[:, :, 40, 130, 82], axis=0)))
    #plt.plot(tdata[0, 276, 40, 130])
    #plt.plot(tdata[0, 1145, 40, 130])
    #plt.show()
    fig, ax = plt.subplots(2, numsets, figsize=(20, 5)) 
    for iset in range(numsets):
        vidx = poslists[iset][0]
        hidx = poslists[iset][1]
        ax[0, iset].imshow(tdata[0, -2].sum(axis=-1))
        ax[0, iset].set_xlabel('x / ch')
        ax[0, iset].set_ylabel('y / ch')
        ax[0, iset].annotate('sample point', xy=(hidx, vidx),
                             xytext=(hidx+5, vidx+5), textcoords='data',
                             color='white',
                             arrowprops=dict(arrowstyle="->", color='white'))
        ax[1, iset].plot(maxstd[vidx, hidx], label='max std of test data sets')
        #ax[1, iset].scatter(np.arange(std.shape[-1]), std[-1, vidx, hidx])
        ax[1, iset].plot(std[-2, vidx, hidx], label='std of sample_expt_noisy')
        ax[1, 0].set_ylabel('std of inference')
        ax[1, 2].set_xlabel('tof / ch')
        #ax[2, iset].plot(np.arange(std.shape[-1]), sample_expt_noisy[vidx, hidx], label='op sample spectrum')
        #for i in range(10):
        #    ax[2, iset].plot(tdata[0, i, vidx, hidx], label='sample spectrum in test data set ' + str(i))
        #    #ax[2, iset].plot(tdata[0, i, vidx, hidx])
        #ax[2, iset].plot(np.max(np.max(tdata[:, :, vidx, hidx], axis=0), axis=0), label='max intensities of test data sets')
        #ax[2, iset].plot(np.min(np.min(tdata[:, :, vidx, hidx], axis=0), axis=0), label='min intensities of test data sets')
        #ax[2, 0].set_ylabel('neutron count')
        #ax[2, 2].set_xlabel('tof / ch')
    #plt.legend()
    #plt.savefig('std_op.png')
    plt.show()


def check_exptspectra():
    tdata = get_tdata()
    print(tdata.shape)
    for i in range(20):
        icol = np.random.randint(tdata.shape[2])
        irow = np.random.randint(tdata.shape[3])
        plt.plot(tdata[0, -2, icol, irow]/tdata[0, -1, icol, irow]*315715./553690)
    plt.savefig('exptspectra_samples.png')


check_std()
#check_exptspectra()
