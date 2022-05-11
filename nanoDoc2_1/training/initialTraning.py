import pyarrow.parquet as pq
import matplotlib.pyplot as plt
import numpy as np
# from tensorflow.keras.utils.training_utils import multi_gpu_model
import tensorflow
import pandas as pd
from tensorflow.keras.callbacks import ModelCheckpoint
from numba import jit

DATA_LENGTH = 1024

def loadpq(path, samplesize):

    df = pq.read_table(path).to_pandas()
    return df

def prepData(df1):

    train_x = []
    test_x = []
    train_y = []
    test_y = []
    totalcnt = 0
    labelidx = df1["fmer"].unique().tolist()
    print(labelidx)
    print(len(labelidx))

    signallen = 1024

    totalcnt = 0
    for idx, row in df1.iterrows():

        flg = labelidx.index(row[1])
        signal = np.array(row[3])
        # count = len(signal)/signallen
        # print(count)

        if totalcnt%12 > 10:
            test_x.extend(signal)
            test_y.append(flg)
        #elif totalcnt%12 < 5 :
        else:
            train_x.extend(signal)
            train_y.append(flg)

        # if cnt == 30000:
        #     break
        totalcnt += 1

    print("totalcnt", totalcnt)
    train_x = np.array(train_x)
    test_x = np.array(test_x)
    train_y = np.array(train_y)
    test_y = np.array(test_y)
    num_classes = np.unique(train_y).size
    print(type(train_x))

    print("train_x.shape", train_x.shape)
    print("test_x.shape", test_x.shape)

    print(num_classes, 'classes')

    print('y_train shape:', train_y.shape)
    print('y_test shape:', test_y.shape)
    train_x = train_x /255
    test_x = test_x /255

    train_x = np.reshape(train_x, (-1, DATA_LENGTH, 1))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH, 1))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    train_y = tensorflow.keras.utils.to_categorical(train_y, num_classes)
    test_y = tensorflow.keras.utils.to_categorical(test_y, num_classes)

    print('train_x:', train_x.shape)
    print('train_y:', train_y.shape)
    print('test_x shape:', test_x.shape)
    print('test_y shape:', test_y.shape)

    return train_x, test_x, train_y, test_y, num_classes

import tensorflow as tf
import os


def main(s_data, s_out):

    with tf.device('/GPU:1'):

        os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        #tf.keras.mixed_precision.experimental.set_policy('mixed_float16')
        _main(s_data, s_out)

import nanoDoc2_1.network.CnnWavenet as CnnWavenet
def _main(s_data, s_out):

    batch_size = 1024
    gpu_count = 1
    samplesize = 3000

    shape1 = (None, DATA_LENGTH, 1)
    df = loadpq(s_data, samplesize)
    train_x, test_x, train_y, test_y, num_classes = prepData(df)

    print(len(train_x), len(test_x), len(train_y), len(test_y), num_classes)

    # with tf.device("/cpu:0"):

    model = CnnWavenet.build_network(shape=shape1, num_classes=num_classes)
    model.summary()

    outweight = s_out+"/weightwn_keras.hdf"
    modelCheckpoint = ModelCheckpoint(filepath=outweight,
                                      monitor='val_accuracy',
                                      verbose=1,
                                      save_best_only=True,
                                      save_weights_only=True,
                                      mode='max',
                                      period=1)

    model.compile(loss='categorical_crossentropy',
                  optimizer=tensorflow.keras.optimizers.Adam(lr=0.0002, beta_1=0.9, beta_2=0.999, epsilon=None,
                                                             decay=0.0,
                                                             amsgrad=False),
                  # optimizer=opt,
                  metrics=['accuracy'])

    history = model.fit(train_x, train_y, epochs=100, batch_size=batch_size, verbose=1,
              shuffle=True, validation_data=(test_x, test_y),callbacks=[modelCheckpoint])

    historypath ="/data/nanopore/nanoDoc2_1/testrun/lealent_log.txt"
    hist_df = pd.DataFrame(history.history)
    hist_df.to_csv(historypath)


if __name__ == '__main__':

    #s_data = "/data/nanopore/nanoDoc2/5000each.pq"
    s_data = "/data/nanopore/nanoDoc2_1/1200signal.pq"
    s_out = "/data/nanopore/nanoDoc2_1/testrun/weight/"
    main(s_data, s_out)