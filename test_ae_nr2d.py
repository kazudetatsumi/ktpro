#!/usr/bin/env python
# This script is a test denoising auto encoder code. 
# Sorry for dirty coding!
# Probably, I took a similar model from
# https://elix-tech.github.io/ja/2016/07/17/autoencoder.html
# and modified the model so as to use the data containing raw counts without
# any normalization.
# For denoising neutron reflectometry 2D data.
# Kazuyoshi TATSUMI 2024/02/17
#import matplotlib.pyplot as plt
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
# importing  a data plotting script written by Hiroyuki AOKI.
import nr2d


def evaluate_PSNRandSSIM(_noise, _orig, _denoise):
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
np.random.seed(314)

batch_size = 32
num_classes = 10
epochs = 100
saveDir = "/home/kazu/ae/"
dataDir = "/home/kazu/2D-NR_2frame/"
gathered_data = dataDir + "gather.pkl"
if not os.path.isdir(saveDir):
    os.makedirs(saveDir)

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

datasets['target'] = np.log1p(datasets['target'])
datasets['noisy'] = np.log1p(datasets['noisy'])
_test = datasets['target'][9]
_ntest = datasets['noisy'][9]
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
x = Conv2D(28, (3, 3), padding='same')(input_img)
#x = BatchNormalization()(x)
x = Activation('relu')(x)
x = MaxPooling2D((2, 2), padding='same')(x)
x = Conv2D(28, (3, 3), padding='same')(x)
#x = BatchNormalization()(x)
x = Activation('relu')(x)
encoded = MaxPooling2D((2, 2), padding='same')(x)

x = Conv2D(28, (3, 3), padding='same')(encoded)
#x = BatchNormalization()(x)
x = Activation('relu')(x)
x = UpSampling2D((2, 2))(x)
x = Conv2D(28, (3, 3), padding='same')(x)
#x = BatchNormalization()(x)
x = Activation('relu')(x)
x = UpSampling2D((2, 2))(x)
x = Conv2D(1, (3, 3), padding='same')(x)
#x = BatchNormalization()(x)
# 256 gray scale is handled with sigmoid, but not for neutron count data.
#decoded = Activation('sigmoid')(x)
decoded = Activation('linear')(x)

model = Model(input_img, decoded)
model.load_weights(saveDir +
                   "AutoEncoder_nr_denoise_weights.18-0.00-0.00.hdf5")

model.compile(optimizer='adam', loss='mean_squared_error')
model.summary()


es_cb = EarlyStopping(monitor='val_loss', patience=2, verbose=1, mode='auto')
chkpt = saveDir + 'AutoEncoder_nr_denoise_weights.' +\
        '{epoch:02d}-{loss:.2f}-{val_loss:.2f}.hdf5'
cp_cb = ModelCheckpoint(filepath=chkpt, monitor='val_loss', verbose=1,
                        save_best_only=True, mode='auto')

# If u need training, uncomment these lines:
#history = model.fit(x_train_noisy, x_train,
#                    batch_size=batch_size,
#                    epochs=epochs,
#                    verbose=1,
#                    validation_data=(x_val_noisy, x_val),
#                    callbacks=[es_cb, cp_cb],
#                    shuffle=True)

score = model.evaluate(x_test_noisy, x_test, verbose=1)
print(score)
nrtest = model.predict(x_test_noisy)
nrval = model.predict(x_val_noisy)
print("nr_test: {0}\nnr_val: {1}".format(np.average(nrtest), np.average(nrval))
      )


def showOrigDec(orig, noise, denoise, num=2):
    for i in range(num):
        _noise = np.expm1(noise[i]).squeeze()
        _orig = np.expm1(orig[i]).squeeze()
        _denoise = np.expm1(denoise[i]).squeeze()
        nr2d.plot(_noise, _orig, _denoise, savefig='test_ae_'+str(i)+'.png')
        evaluate_PSNRandSSIM(np.log1p(_noise), np.log1p(_orig),
                             np.log1p(_denoise))


showOrigDec(x_test, x_test_noisy, nrtest, num=5)
