#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from keras.models import load_model
from keras.layers import Input, Dense, Conv2D, MaxPooling2D, UpSampling2D,\
        BatchNormalization, Activation
from keras.models import Model
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.optimizers import Adam
from skimage import metrics
import os
import pickle
import numpy as np
import pickle
import sys
sys.path.append("/home/kazu")
import nr2d
sys.path.append("/home/kazu/decnn_sources/DnCNN-tf2")
from DnCNN import DnCNN

def evaluate_PSNRandSSIM(_noise, _orig, _denoise):
    _noise = _noise
    _orig = _orig
    _denoise = _denoise
    val_min = _orig.min()
    val_range = _orig.max() - val_min
    __orig = (_orig - val_min)/val_range
    __noise = (_noise - val_min)/val_range
    __denoise = (_denoise - val_min)/val_range
    print("PSNR between __orig and __noise:",
          metrics.peak_signal_noise_ratio(__orig, __noise,
                                          data_range=val_range))
    print("PSNR between __orig and __denoise:",
          metrics.peak_signal_noise_ratio(__orig, __denoise,
                                          data_range=val_range))
    print("SSIM between __orig and __noise:",
          metrics.structural_similarity(__orig, __noise,
                                        data_range=val_range))
    print("SSIM between __orig and __denoise:",
          metrics.structural_similarity(__orig, __denoise,
                                        data_range=val_range))

os.environ["KERAS_BACKEND"] = "tensorflow"
kerasBKED = os.environ["KERAS_BACKEND"]
print(kerasBKED)

batch_size = 32
num_classes = 10
epochs = 100
saveDir = "/home/kazu/decnn/"
dataDir = "/home/kazu/2D-NR_2frame/"
gathered_data = dataDir + "gather.pkl"
if not os.path.isdir(saveDir):
    os.makedirs(saveDir)

np.random.seed(314)
if not os.path.exists(gathered_data):
    files = sorted(os.listdir(dataDir))
    num_data = int(len(files)/2)
    for ifdx, file in enumerate(files):
        print(ifdx, file)
        if ifdx == 0:
            datasize = np.load(dataDir + file).shape
            data = np.zeros((num_data*2, datasize[0], datasize[1]))
        data[ifdx] = np.load(dataDir + file)
    data = np.expand_dims(data, 3)
    with open(gathered_data, 'wb') as f:
        datasets = {}
        datasets['noisy'] = data[::2]
        datasets['target'] = data[1::2]
        pickle.dump(datasets, f, 4)
else:
    with open(gathered_data, 'rb') as f:
        datasets = pickle.load(f)
    num_data = datasets['target'].shape[0]

datasets['target'] = np.log(datasets['target'] + 1.)
datasets['noisy'] = np.log(datasets['noisy'] + 1.)
_test = datasets['target'][9]
_ntest = datasets['noisy'][9]
#for i in range(10):
#    plt.plot(_test[0+int(i/10*_test.shape[0])])
#plt.imshow(datasets['target'][9])
#plt.show()
#plt.imshow(datasets['noisy'][9])
#plt.show()
#for i in range(10):
#    plt.plot(_ntest[0+int(i/10*_ntest.shape[0])])
#plt.show()
#print("max_intensity (log)", max_intensity)
#print("Number of data:", num_data)
train_size = 15000
test_size = num_data - train_size
train_idx = np.random.choice(num_data, train_size,  replace=False)
test_idx = np.setdiff1d(np.arange(num_data), train_idx)
x_train = datasets['target'][train_idx]
x_test = datasets['target'][test_idx]
x_train_noisy = datasets['noisy'][train_idx]
x_test_noisy = datasets['noisy'][test_idx]


x_val = x_test[:7000]
x_test = x_test[7000:]
x_val_noisy = x_test_noisy[:7000]
x_test_noisy = x_test_noisy[7000:]


x_train = x_train.astype('float32')
x_test = x_test.astype('float32')
x_val = x_val.astype('float32')
x_train_noisy = x_train_noisy.astype('float32')
x_test_noisy = x_test_noisy.astype('float32')
x_val_noisy = x_val_noisy.astype('float32')


print("validation data: {0} \ntest data: {1}".format(x_val.shape, x_test.shape)
      )

input_img = Input(shape=x_train[0].shape)
model = DnCNN(depth=7).model()
#model.load_weights(saveDir + "DnCNN_nr_denoise_weights.17-0.00-0.00.hdf5")

model.compile(optimizer='adam', loss='mean_squared_error')
#model.summary()


es_cb = EarlyStopping(monitor='val_loss', patience=2, verbose=1, mode='auto')
chkpt = saveDir + 'DnCNN_nr_denoise_weights_7d_bnrm.' +\
        '{epoch:02d}-{loss:.2f}-{val_loss:.2f}.hdf5'
cp_cb = ModelCheckpoint(filepath=chkpt, monitor='val_loss', verbose=1,
                        save_best_only=True, mode='auto')

history = model.fit(x_train_noisy, x_train,
                    batch_size=batch_size,
                    epochs=epochs,
                    verbose=1,
                    validation_data=(x_val_noisy, x_val),
                    callbacks=[es_cb, cp_cb],
                    shuffle=True)

score = model.evaluate(x_test_noisy, x_test, verbose=1)
print(score)
nrtest = model.predict(x_test_noisy)
nrval = model.predict(x_val_noisy)
print("nr_test: {0}\nnr_val: {1}".format(np.average(nrtest), np.average(nrval))
      )


def showOrigDec(orig, noise, denoise, num=2):
    #import matplotlib.pyplot as plt 
    n = num 
    plt.figure(figsize=(8, 8)) 

    #for i in range(n):
    #    # display original
    #    ax = plt.subplot(3, n, i+1)
    #    print(orig.shape)
    #    _orig = (10**orig[i] - 1.).reshape(x_train[0].shape)
    #    #_orig = orig[i].reshape(x_train[0].shape)
    #    for j in range(10):
    #        plt.plot(_orig[int(j/10*_orig.shape[0])])
    #    #plt.imshow(datasets['target'][9])
    #    #plt.imshow((10**orig[i]-1.).reshape(x_train[0].shape))
    #    #plt.gray()
    #    #ax.get_xaxis().set_visible(False)
    #    #ax.get_yaxis().set_visible(False)
    #    #for j in range(10):
    #    #plt.plot(10**(orig[i]-1.).reshape(x_train[0].shape)[int(j/x_train[0].shape[0])])
#
#        # display noisy image
#        ax = plt.subplot(3, n, i + 1 + n)
#        _noise = (10**noise[i] - 1.).reshape(x_train[0].shape)
#        #_noise = noise[i].reshape(x_train[0].shape)
#        for j in range(10):
#            plt.plot(_noise[int(j/10*_noise.shape[0])])
#        #plt.imshow(10**(noise[i]-1.).reshape(x_train[0].shape))
#        #plt.gray()
#        #ax.get_xaxis().set_visible(False)
#        #ax.get_yaxis().set_visible(False)
#
#        # display denoised image
#        ax = plt.subplot(3, n, i + 1 + n + n)
#        _denoise = (10**denoise[i] - 1.).reshape(x_train[0].shape)
#        #_denoise = denoise[i].reshape(x_train[0].shape)
#        for j in range(10):
#            plt.plot(_denoise[int(j/10*_denoise.shape[0])])
#        #plt.imshow(10**(denoise[i]-1.).reshape(x_train[0].shape)) 
#        #plt.gray()
#        #ax.get_xaxis().set_visible(False)
##        #ax.get_yaxis().set_visible(False)
#    plt.show()

    for i in range(n):
        # display original
        ###ax = plt.subplot(n, 3, i*3 + 1)
        _noise = (10**noise[i] - 1.).squeeze()
        #_orig = orig[i].reshape(x_train[0].shape)
        #for j in range(10):
        #    plt.plot(_orig[int(j/10*_orig.shape[0])])
        ###plt.imshow(_noise.T, origin='lower', norm=LogNorm())
        #plt.imshow((10**orig[i]-1.).reshape(x_train[0].shape))
        #plt.gray()
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)
        #for j in range(10):
        #plt.plot(10**(orig[i]-1.).reshape(x_train[0].shape)[int(j/x_train[0].shape[0])])

        # display noisy image
        ###ax = plt.subplot(n, 3, i*3 + 2)
        _orig = (10**orig[i] - 1.).squeeze()
        #_noise = noise[i].reshape(x_train[0].shape)
        #for j in range(10):
        #    plt.plot(_noise[int(j/10*_noise.shape[0])])
        ###plt.imshow(_orig.T, origin='lower', norm=LogNorm())
        #plt.imshow(10**(noise[i]-1.).reshape(x_train[0].shape))
        #plt.gray()
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)

        # display denoised image
        ###ax = plt.subplot(n, 3, i*3  + 3)
        _denoise = (10**denoise[i] - 1.).squeeze()
        #_denoise = denoise[i].reshape(x_train[0].shape)
        #for j in range(10):
        #    plt.plot(_denoise[int(j/10*_denoise.shape[0])])
        ###plt.imshow(_denoise.T, origin='lower', norm=LogNorm())
        #plt.imshow(10**(denoise[i]-1.).reshape(x_train[0].shape)) 
        #plt.gray()
        #ax.get_xaxis().set_visible(False)
        #ax.get_yaxis().set_visible(False)
    ###plt.show()
        nr2d.plot(_noise, _orig, _denoise, savefig='test_dncnn_'+str(i)+'.png')
        evaluate_PSNRandSSIM(np.log(_noise[50]+1.), np.log(_orig[50]+1.), np.log(_denoise[50]+1.))
    #_noise = _noise[50]
    #_orig = _orig[50]
    #_denoise = _denoise[50]
    #val_min = _orig.min()
    #val_range = _orig.max() - val_min
    #__orig = (_orig - val_min)/val_range
    #__noise = (_noise - val_min)/val_range
    #__denoise = (_denoise - val_min)/val_range
    #print("PSNR between __orig and __noise:",
    #      metrics.peak_signal_noise_ratio(__orig, __noise,
    #                                      data_range=val_range))
    #print("PSNR between __orig and __denoise:",
    #      metrics.peak_signal_noise_ratio(__orig, __denoise,
    #                                      data_range=val_range))
    #print("SSIM between __orig and __noise:",
    #      metrics.structural_similarity(__orig, __noise,
    #                                    data_range=val_range))
    #print("SSIM between __orig and __denoise:",
    #      metrics.structural_similarity(__orig, __denoise,
    #                                    data_range=val_range))


showOrigDec(x_test, x_test_noisy, nrtest, num=5)


