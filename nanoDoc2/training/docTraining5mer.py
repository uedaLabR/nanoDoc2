import pyarrow.parquet as pq

import numpy as np
import tensorflow
from keras.layers import GlobalAveragePooling1D,Dense
from keras.models import Model
from keras.engine.network import Network
from keras.optimizers import SGD
from keras import backend as K
import itertools
import keras

from nanoDoc2.network import cnnwavenettrace_keras

DATA_LENGTH = 40


def getNearExculudeSelf(nuc):

    first = nuc[0]
    last = nuc[4]
    nucs = ('A', 'T', 'C', 'G')
    gen = []
    for n1, n2, n3 in itertools.product(nucs, nucs, nucs):

        nc = first + n1 + n2 + n3 + last
        if nc != nuc:
            gen.append(nc)

    return gen

import os
def prepDataNear(s_data,samplesize, nuc):


    train_x = []
    test_x = []
    train_y = []
    test_y = []
    p_flg = 0
    flg = 0
    totalcnt = 0

    nucs = getNearExculudeSelf(nuc)
    s_data ="/data/nanopore/nanoDoc2/100traceeach"
    paths = list(map(lambda x: s_data + "/" + x + ".pq", nucs))

    samplecnt = 0
    for path in paths:

        if not os.path.exists(path):
            continue
        print("reading", path)

        #         print(df)
        df = pq.read_table(path).to_pandas()
        #         print(df)
        cnt = 0
        for idx, row in df.iterrows():

            flg = row[1]
            trace = np.array(row[2])


            testidx = (idx % 5 == 0)
            if testidx:
                test_x.append(trace)
                test_y.append(samplecnt)
            else:
                train_x.append(trace)
                train_y.append(samplecnt)

            cnt = cnt + 1
            totalcnt = totalcnt + 1

            if cnt % 12000 == 0:
                print(samplecnt, totalcnt, path, totalcnt, idx, row)
            if cnt == samplesize:
                break

        samplecnt = samplecnt + 1


    train_x = np.array(train_x)
    test_x = np.array(test_x)
    train_y = np.array(train_y)
    test_y = np.array(test_y)
    num_classes = np.unique(train_y).size

    train_x = train_x /255
    test_x = test_x /255
    DATA_LENGTH = 40
    train_x = np.reshape(train_x, (-1, DATA_LENGTH, 4))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH, 4))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    num_classes = 63
    train_y = tensorflow.keras.utils.to_categorical(train_y, num_classes)
    test_y = tensorflow.keras.utils.to_categorical(test_y, num_classes)

    return train_x, test_x, train_y, test_y, num_classes



def prepData(s_data,samplesize, nuc):

    train_x = []
    test_x = []
    train_y = []
    test_y = []
    totalcnt = 0
    totalcnt = 0
    path = s_data + "/" + nuc + ".pq"
    print("reading", path)
    df = pq.read_table(path).to_pandas()
    totalcnt = 0
    for idx, row in df.iterrows():

        trace = np.array(row[2])
        flg = row[1]
        if totalcnt%5 == 0:
            test_x.extend(trace)
            test_y.append(flg)
        else:
            train_x.extend(trace)
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

    DATA_LENGTH_UNIT = 8*5
    train_x = np.reshape(train_x, (-1, DATA_LENGTH_UNIT, 4))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH_UNIT, 4))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    # train_y = tensorflow.keras.utils.to_categorical(train_y, num_classes)
    # test_y = tensorflow.keras.utils.to_categorical(test_y, num_classes)

    print('train_x:', train_x.shape)
    print('train_y:', train_y.shape)
    print('test_x shape:', test_x.shape)
    print('test_y shape:', test_y.shape)

    return train_x, test_x, train_y, test_y, num_classes




batchsize = 128
feature_out = 1280  # secondary network out for MobileNet
alpha = 0.5  # for MobileNetV2
lambda_ = 0.1  # for compact loss

import tensorflow as tf
def train(s_data, s_out, nuc,bestwight ,samplesize , epoch_num):

    with tf.device('/GPU:1'):

        os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        #tf.keras.mixed_precision.experimental.set_policy('mixed_float16')
        _train(s_data, s_out, nuc,bestwight ,samplesize , epoch_num)


import nanoDoc2.network.cnnwavenet_decfilter as cnnwavenet_decfilter
def _train(s_data, s_out, nuc,bestwight ,samplesize , epoch_num):

    if not os.path.exists(s_out + '/' + nuc):
        os.makedirs(s_out + '/' + nuc)

    logf = open(s_out + '/' + nuc+'/log.txt', mode='w')
    pqpath = s_data + "/" + nuc + ".pq"
    if not os.path.exists(pqpath):
      print("no pq for",pqpath)
      return

    num_classes_org = 1024
    num_classes = 63
    shape1 = (None, DATA_LENGTH, 4)
    optimizer = SGD(lr=5e-5, decay=0.00005)

    model = cnnwavenet_decfilter.build_network(shape=shape1, num_classes=num_classes_org)
    model.load_weights(bestwight)
    #model.layers.pop()  # remove last layer
    #model.layers.pop()  # remove last layer
    #model.layers.pop()  # remove last layer


    for layer in model.layers:
        if layer.name == "conv1d_68":
            break
        else:
            layer.trainable = False

    flat = GlobalAveragePooling1D()(model.layers[-4].output)
    model_t = Model(inputs=model.input, outputs=flat)
    model_r = Network(inputs=model_t.input,
                      outputs=flat,
                      name="shared_layer")

    prediction = Dense(num_classes, activation='softmax')(model_t.output)
    model_r = Model(inputs=model_r.input, outputs=prediction)


    model_r.compile(optimizer=optimizer, loss="categorical_crossentropy")
    model_t.compile(optimizer=optimizer, loss=original_loss)

    model_r.summary()
    model_t.summary()

    print("reading data")
    x_ref, test_x_r, y_ref, test_y_r, num_classes_r = prepDataNear(s_data,samplesize, nuc)
    x_target, test_x, train_y, test_y, num_classes = prepData(s_data,samplesize, nuc)

    ref_samples = np.arange(x_ref.shape[0])

    loss, loss_c = [], []
    epochs = []
    print("training...")

    for epochnumber in range(epoch_num):

        x_r, y_r, lc, ld = [], [], [], []

        np.random.shuffle(x_target)
        np.random.shuffle(ref_samples)

        for i in range(len(x_target)):
            x_r.append(x_ref[ref_samples[i]])
            y_r.append(y_ref[ref_samples[i]])
        x_r = np.array(x_r)
        y_r = np.array(y_r)

        for i in range(int(len(x_target) / batchsize)):

            batch_target = x_target[i * batchsize:i * batchsize + batchsize]
            batch_ref = x_r[i * batchsize:i * batchsize + batchsize]
            batch_y = y_r[i * batchsize:i * batchsize + batchsize]
            # target data
            lc.append(model_t.train_on_batch(batch_target, np.zeros((batchsize, feature_out))))

            # reference data
            ld.append(model_r.train_on_batch(batch_ref, batch_y))

        loss.append(np.mean(ld))
        loss_c.append(np.mean(lc))
        epochs.append(epochnumber)

        print("epoch : {} ,Descriptive loss : {}, Compact loss : {}".format(epochnumber + 1, loss[-1], loss_c[-1]))
        logf.write("epoch : {} ,Descriptive loss : {}, Compact loss : {}".format(epochnumber + 1, loss[-1], loss_c[-1])+"\n")
        if epochnumber % 1 == 0:
            model_t.save_weights(s_out +'/'+ nuc + '/model_t_ep_{}.h5'.format(epochnumber))
            model_r.save_weights(s_out +'/'+ nuc + '/model_r_ep_{}.h5'.format(epochnumber))

    logf.close()

batchsize = 64
# feature_out = 512 #secondary network out for VGG16
# feature_out = 1280 #secondary network out for MobileNet
feature_out = 1536  # secondary network out for Inception Resnetv2
alpha = 0.5  # for MobileNetV2
lambda_ = 0.1  # for compact loss
classes = 63


def original_loss(y_true, y_pred):

    lc = 1 / (classes * batchsize) * batchsize ** 2 * K.sum((y_pred - K.mean(y_pred, axis=0)) ** 2, axis=[1]) / (
                (batchsize - 1) ** 2)
    return lc


import sys

if __name__ == '__main__':

    #train(sys.argv[2], sys.argv[3], sys.argv[1], 10)
    nucs = ('A', 'T', 'C', 'G')
    cnt = 0
    for n1, n2, n3, n4, n5 in itertools.product(nucs, nucs, nucs, nucs, nucs):

        s_data = "/data/nanopore/nanoDoc2/5000traceeach"
        s_out = "/data/nanopore/IVT_Trace/weight/docweight"
        nuc = n1+n2+n3+n4+n5
        bestwight = "/data/nanopore/IVT_Trace/weight_dec/weightwn_keras_dec.hdf"
        epoch_num = 10
        samplesize = 5000
        print(nuc,cnt)
        train(s_data, s_out, nuc, bestwight, samplesize, epoch_num)
        cnt+=1
        #break