#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from keras.models import load_model
from keras.layers import Input, Dense, Conv2D, MaxPooling2D, UpSampling2D,\
        BatchNormalization, Activation
from keras.models import Model
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.optimizers import Adam
import tensorflow as tf
from skimage import metrics
import os
import pickle
import numpy as np
import pickle
import sys, glob, tqdm
from datetime import datetime, timedelta
import nr2d
# from DnCNN2 import DnCNN

os.environ["KERAS_BACKEND"] = "tensorflow"
kerasBKED = os.environ["KERAS_BACKEND"]
print(kerasBKED)

"""
https://danijar.com/structuring-models/
"""
class DnCNN(tf.keras.Model):
  
    def __init__(self, depth=5, grayscale=True):
        super(DnCNN,self).__init__()
        # Network params
        self.channels = 1 if grayscale else 3
        self.depth = depth

    def call(self,input_tensor,training=True):
        # First Convolution Layer with Conv and ReLU: kernel (q,t)
        x = tf.keras.layers.Conv2D(64,(3,3),padding="same",kernel_initializer='Orthogonal')(input_tensor)
        x = tf.keras.activations.relu(x)

        # Add Conv+Batch_Norm+ReLU for layers 2 to (depth-1)
        for i in range(self.depth - 1):
            x = tf.keras.layers.Conv2D(64,(5,9),padding="same",kernel_initializer='Orthogonal')(x)
            x = tf.keras.layers.BatchNormalization(epsilon=0.0001)(x,training=training)
            x = tf.keras.activations.relu(x)

        # The final conv layer will use only 1 filter to recontruct the original image
        x = tf.keras.layers.Conv2D(1,(3,3),padding="same",kernel_initializer='Orthogonal')(x)

        # Subtract the predicted noise from the noisy input image
        x = tf.keras.layers.Subtract()([input_tensor,x]) #input - noise

        return x
    
    def model(self):
        # Funtion to build the model
        x = tf.keras.Input(shape=(None,None,self.channels))
        return tf.keras.Model(inputs=[x],outputs= self.call(x) )


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

def eval_quality(truth, data, pred):
    max_val   = np.max([np.max(truth), np.max(data), np.max(pred)])
    min_val   = np.min([np.min(truth), np.min(data), np.min(pred)])
    psnr_data = metrics.peak_signal_noise_ratio(truth, data, data_range=max_val-min_val)
    psnr_pred = metrics.peak_signal_noise_ratio(truth, pred, data_range=max_val-min_val)
    ssim_data = metrics.structural_similarity(truth, data, data_range=max_val-min_val)
    ssim_pred = metrics.structural_similarity(truth, pred, data_range=max_val-min_val)
    return psnr_data, psnr_pred, ssim_data, ssim_pred

def showOrigDec(orig, noise, denoise, num=2, path='.', prefix='test'):
    n = num 
    plt.figure(figsize=(8, 8)) 
    for i in range(n):
        _noise = (10**noise[i] - 1.).squeeze()
        _orig = (10**orig[i] - 1.).squeeze()
        _denoise = (10**denoise[i] - 1.).squeeze()
        nr2d.plot(_noise, _orig, _denoise, savefig=f'{path}/{prefix}_'+str(i).zfill(3)+'.png')
        # evaluate_PSNRandSSIM(np.log(_noise[50]+1.), np.log(_orig[50]+1.), np.log(_denoise[50]+1.))
        eval_quality(np.log(_noise[50]+1.), np.log(_orig[50]+1.), np.log(_denoise[50]+1.))

def showOrigDec2(orig, noise, denoise, num=2, path='.', prefix='test'):
    plt.figure(figsize=(8, 8)) 
    print(f'      {"PSNR".ljust(16)}|  {"SSIM".ljust(8)}')
    print( '------'+'-'*16+'+'+'-'*18)
    print(f'      {"Data".ljust(8)}{"Pred".ljust(8)}|  {"Data".ljust(9)}{"Pred".ljust(9)}')
    print( '------'+'-'*16+'+'+'-'*18)
    result = []
    for i in range(len(orig)):
        qual = eval_quality(orig[i].squeeze(), noise[i].squeeze(), denoise[i].squeeze())
        result.append(list(qual))
        if i<num:
            print(f'{str(i).rjust(3)})  {qual[0]:.2f}   {qual[1]:.2f}   |  {qual[2]:.4f}   {qual[3]:.4f}')
            _noise   = np.expm1(noise[i]).squeeze()
            _orig    = np.expm1(orig[i]).squeeze()
            _denoise = np.expm1(denoise[i]).squeeze()
            nr2d.plot(_noise, _orig, _denoise, savefig=f'{path}/{prefix}_'+str(i).zfill(3)+'.png')
        elif i==num:
            print('         ........')
    res = np.mean(result, axis=0)
    print(f'Avg.  {res[0]:.2f}   {res[1]:.2f}   |  {res[2]:.4f}   {res[3]:.4f}')

def choice_idx(num, train_ratio, test=0.1):
    train_idx = np.random.choice(num, int(num*train_ratio), replace=False)
    if test>=1:
        test_num = max(int(test), 10)
    else:
        test_num = max(int(num*(1-train_ratio)*test), 10)
    val_num  = num - len(train_idx) - test_num
    digit = int(np.log10(num*train_ratio))+1
    print(f'Training:   {str(int(num*train_ratio)).rjust(digit)}')
    print(f'Validation: {str(val_num).rjust(digit)}')
    print(f'Test:       {str(test_num).rjust(digit)}')
    tmp      = np.setdiff1d(np.arange(num), train_idx)
    val_idx  = np.random.choice(tmp, val_num, replace=False)
    test_idx = np.setdiff1d(tmp, val_idx)
    return train_idx, val_idx, test_idx


def train(dataDir, batch_size=1024, epochs=128, saveDir='./'):
    num_classes = 10
    # saveDir = f'./result_{datetime.strftime(datetime.now()+timedelta(hours=9),"%y%m%d-%H%M%S")}/'
    # dataDir = '/home/aoki/storage1/nfs_share/2D-NR_2frame/'
    if not os.path.isdir(saveDir):
        os.makedirs(saveDir)

    # np.random.seed(314)
    files = sorted(os.listdir(dataDir))
    num_data = int(len(files)/2)
    for i, f in enumerate(tqdm.tqdm(files, desc='Loading data')):
        if i == 0:
            datasize = np.load(dataDir + f).shape
            data = np.zeros((num_data*2, datasize[0], datasize[1]))
        data[i] = np.load(dataDir + f)
    data = np.expand_dims(data, 3)
    datasets = {}
    datasets['noisy']  = data[::2]
    datasets['target'] = data[1::2]

    datasets['target'] = np.log1p(datasets['target'])
    datasets['noisy']  = np.log1p(datasets['noisy'])

    train_idx, val_idx, test_idx = choice_idx(num_data, 0.8, test=0.02)
    x_train_target = datasets['target'][train_idx]
    x_val_target   = datasets['target'][val_idx]
    x_test_target  = datasets['target'][test_idx]
    x_train_noisy  = datasets['noisy'][train_idx]
    x_val_noisy    = datasets['noisy'][val_idx]
    x_test_noisy   = datasets['noisy'][test_idx]

    x_train_target = x_train_target.astype('float32')
    x_test_target  = x_test_target.astype('float32')
    x_val_target   = x_val_target.astype('float32')
    x_train_noisy  = x_train_noisy.astype('float32')
    x_test_noisy   = x_test_noisy.astype('float32')
    x_val_noisy    = x_val_noisy.astype('float32')


    print("validation data: {0} \ntest data: {1}".format(x_val_target.shape, x_test_target.shape))

    gpus = tf.config.list_logical_devices('GPU')
    # gpus = [gpus[0]]  # GPUを指定
    strategy = tf.distribute.MirroredStrategy(gpus)
    with strategy.scope():
        input_img = Input(shape=x_train_target[0].shape)
        model = DnCNN(depth=7).model()
        #　model.load_weights(saveDir + "DnCNN_nr_denoise_weights.17-0.00-0.00.hdf5")
        optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
        model.compile(optimizer=optimizer, loss='mean_squared_error')
    #model.summary()

    es_cb = EarlyStopping(monitor='val_loss', patience=8, min_delta=0.0001, verbose=1, mode='auto')
    chkpt = saveDir + 'weights.{epoch:02d}-{loss:.2f}-{val_loss:.2f}.hdf5'
    cp_cb = ModelCheckpoint(filepath=chkpt, monitor='val_loss', verbose=1,
                            save_best_only=True, mode='auto')

    result = model.fit(x_train_noisy, x_train_target,
                        batch_size=batch_size,
                        epochs=epochs,
                        verbose=1,
                        validation_data=(x_val_noisy, x_val_target),
                        # callbacks=[es_cb, cp_cb],
                        callbacks=[es_cb],
                        shuffle=True)

    plt.plot(range(1,len(result.history['loss'])+1),     result.history['loss'], label='Training')
    plt.plot(range(1,len(result.history['val_loss'])+1), result.history['val_loss'], label='Validation')
    plt.xlabel('Epoch', fontsize='large')
    plt.ylabel('Loss', fontsize='large')
    plt.yscale('log')
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f'{saveDir}/learning_curve.pdf')

    score = model.evaluate(x_test_noisy, x_test_target, verbose=1)
    print(score)
    nrtest = model.predict(x_test_noisy)
    nrval  = model.predict(x_val_noisy)
    print("nr_test: {0}\nnr_val: {1}".format(np.average(nrtest), np.average(nrval)))

    showOrigDec2(x_test_target, x_test_noisy, nrtest, num=min(20,len(x_test_target)), path=saveDir)

def predict():
    pass

if __name__ == '__main__':
    train(
        dataDir    = '/home/aoki/storage1/nfs_share/2D-NR_2frame/',
        batch_size = 1024,
        epochs     = 256,
        saveDir    = f'./result_{datetime.strftime(datetime.now()+timedelta(hours=9),"%y%m%d-%H%M%S")}/',
    )
