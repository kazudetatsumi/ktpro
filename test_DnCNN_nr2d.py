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
    n = num
    plt.figure(figsize=(8, 8))
    for i in range(n):
        _noise = np.expm1(noise[i]).squeeze()
        _orig = np.expm1(orig[i]).squeeze()
        _denoise = np.expm1(denoise[i]).squeeze()
        nr2d.plot(_noise, _orig, _denoise, savefig='test_dncnn_'+str(i)+'.png')
        evaluate_PSNRandSSIM(np.log1p(_noise), np.log1p(_orig),
                             np.log1p(_denoise))


showOrigDec(x_test, x_test_noisy, nrtest, num=5)
