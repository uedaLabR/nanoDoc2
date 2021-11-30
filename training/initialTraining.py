import pyarrow.parquet as pq
import matplotlib.pyplot as plt
import numpy as np
# from tensorflow.keras.utils.training_utils import multi_gpu_model
import tensorflow
import src.nanoDoc2.cnnnetwork
import os
import pandas as pd
from tensorflow.keras.callbacks import ModelCheckpoint
import src.nanoDoc2.cnnnetwork as cnn_network
from src.nanoDoc2 import cnnwavenet, cnnwavenet_keras

DATA_LENGTH =  420

def loadpq(path, samplesize):

    df = pq.read_table(path).to_pandas()
    return df


def prepData(df1):

    train_x = []
    test_x = []
    train_y = []
    test_y = []
    totalcnt = 0
    labelidx = df1["sixmer"].unique().tolist()

    cnt = 0
    for idx, row in df1.iterrows():

        #flg = labelidx.index(row[0])
        #signal = np.array(list(row[1]))
        flg = row[1]
        signal = np.array(list(row[2]))
        signal = np.array(signal)
        signal = signal.astype('float16') / 255.


        testidx = (cnt % 10 == 0)
        if testidx:
            test_x.extend(signal)
            test_y.append(flg)
        else:
            train_x.extend(signal)
            train_y.append(flg)

        totalcnt = totalcnt + 1
        if totalcnt % 1200 == 0:
            print(totalcnt, idx,row[0],row[1],len(row[2]))
            #print(totalcnt, idx, row[0], len(row[1]))
        cnt = cnt+1
        # if cnt == 30000:
        #     break

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

def _main(s_data, s_out):

    batch_size = 512
    gpu_count = 1
    samplesize = 2400

    shape1 = (None, DATA_LENGTH, 1)
    df = loadpq(s_data, samplesize)
    train_x, test_x, train_y, test_y, num_classes = prepData(df)

    print(len(train_x), len(test_x), len(train_y), len(test_y), num_classes)

    # with tf.device("/cpu:0"):

    #model = cnn_network.build_network(shape=shape1, num_classes=num_classes)
    model = cnnwavenet_keras.build_network(shape=shape1, num_classes=num_classes)
    model.summary()
    # model = multi_gpu_model(model, gpus=gpu_count)  # add
    outweight = s_out+"/weightwn_keras.hdf"
    modelCheckpoint = ModelCheckpoint(filepath=outweight,
                                      monitor='val_accuracy',
                                      verbose=1,
                                      save_best_only=True,
                                      save_weights_only=True,
                                      mode='max',
                                      period=1)

    model.compile(loss='categorical_crossentropy',
                  optimizer=tensorflow.keras.optimizers.Adam(lr=0.0003, beta_1=0.9, beta_2=0.999, epsilon=None,
                                                             decay=0.0,
                                                             amsgrad=False),
                  # optimizer=opt,
                  metrics=['accuracy'])

    history = model.fit(train_x, train_y, epochs=200, batch_size=batch_size, verbose=1,
              shuffle=True, validation_data=(test_x, test_y),callbacks=[modelCheckpoint])

    historypath ="/data/nanopore/IVT/lealent_log.txt"
    hist_df = pd.DataFrame(history.history)
    hist_df.to_csv(historypath)


if __name__ == '__main__':

    s_data = "/data/nanopore/IVT/2400each.pq"
    s_out = "/data/nanopore/IVT/weight/"
    main(s_data, s_out)