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
from nanoDoc2_1.network import CnnWavenetDecDimention

DATA_LENGTH = 1024


def getNearExculudeSelf(nuc):
    first = nuc[0]
    last = nuc[5]
    nucs = ('A', 'T', 'C', 'G')
    gen = []
    for n1, n2, n3, n4 in itertools.product(nucs, nucs, nucs, nucs):

        nc = first + n1 + n2 + n3 + n4 + last
        if nc != nuc:
            gen.append(nc)

    return gen

import os
def prepDataNear(s_data2,samplesize, nuc):


    train_x = []
    test_x = []
    train_y = []
    test_y = []
    p_flg = 0
    flg = 0
    totalcnt = 0

    nucs = getNearExculudeSelf(nuc)
    print(nucs)
    paths = list(map(lambda x: s_data2 + "/" + x + ".pq", nucs))

    samplecnt = 0
    for path in paths:

        if not os.path.exists(path):
            continue
        print("reading", path)

        #         print(df)
        df = pq.read_table(path).to_pandas()
        df = df[0:100]
        #         print(df)
        cnt = 0
        for idx, row in df.iterrows():

            flg = row[1]
            signal = np.array(list(row[3]))
            signal = np.array(signal)
            signal = signal.astype('float16') / 255.

            testidx = (idx % 12 >= 10)
            train_x.append(signal)
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


    train_x = np.reshape(train_x, (-1, DATA_LENGTH, 1))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH, 1))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    num_classes = 255
    train_y = tensorflow.keras.utils.to_categorical(train_y, num_classes)
    test_y = tensorflow.keras.utils.to_categorical(test_y, num_classes)

    return train_x, test_x, train_y, test_y, num_classes


def prepData(s_data,samplesize, nuc):

    train_x = []
    test_x = []
    train_y = []
    test_y = []
    p_flg = 0
    flg = 0
    totalcnt = 0
    path = s_data + "/" + nuc + ".pq"
    print("reading", path)
    df = pq.read_table(path).to_pandas()
    cnt = 0
    for idx, row in df.iterrows():

        flg = 1
        signal = np.array(list(row[3]))
        signal = signal.astype('float16') / 255.
        # print(signal)
        # print(type(signal))
        #testidx = (idx % 12 >= 10)
        train_x.append(signal)
        train_y.append(flg)

        cnt = cnt + 1
        totalcnt = totalcnt + 1

        if cnt % 12000 == 0:
            print(totalcnt, path, totalcnt, idx, row)
        if cnt == samplesize:
            break


    train_x = np.array(train_x)
    test_x = np.array(test_x)
    train_y = np.array(train_y)
    test_y = np.array(test_y)
    num_classes = np.unique(train_y).size


    train_x = np.reshape(train_x, (-1, DATA_LENGTH, 1))
    test_x = np.reshape(test_x, (-1, DATA_LENGTH, 1))
    train_y = np.reshape(train_y, (-1, 1,))
    test_y = np.reshape(test_y, (-1, 1,))

    test_y = test_y - 1
    train_y = train_y - 1
    # train_y = tensorflow.keras.utils.to_categorical(train_y, num_classes)
    # test_y = tensorflow.keras.utils.to_categorical(test_y, num_classes)

    return train_x, test_x, train_y, test_y, num_classes


batchsize = 128
feature_out = 1280  # secondary network out for MobileNet
alpha = 0.5  # for MobileNetV2
lambda_ = 0.1  # for compact loss

import tensorflow as tf
def train(s_data1,s_data2, s_out, nuc,weightdir ,samplesize , epoch_num,device):

    with tf.device(device):

        os.environ["TF_FORCE_GPU_ALLOW_GROWTH"] = "true"
        os.environ["TF_GPU_ALLOCATOR"] = "cuda_malloc_async"
        _train(s_data1,s_data2, s_out, nuc,weightdir ,samplesize , epoch_num)


def _train(s_data1,s_data2, s_out, nuc,weightdir,samplesize , epoch_num):

    if not os.path.exists(s_out + '/' + nuc):
        os.makedirs(s_out + '/' + nuc)

    logf = open(s_out + '/' + nuc+'/log.txt', mode='w')
    pqpath = s_data1 + "/" + nuc + ".pq"
    if not os.path.exists(pqpath):
      print("no pq for",pqpath)
      return

    num_classes_org = 4078
    num_classes = 255
    shape1 = (None, DATA_LENGTH, 1)
    optimizer = SGD(lr=5e-5, decay=0.00005)

    bestwight = weightdir + "/weightwn_dec.hdf"

    model = CnnWavenetDecDimention.build_network(shape=shape1, num_classes=num_classes_org)
    model.load_weights(bestwight)


    for layer in model.layers:
        if layer.name == "conv1d_99":
            break
        else:
            layer.trainable = False

    flat = model.layers[-11].output
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
    x_ref, test_x_r_o, y_ref, test_y_r_o, num_classes_r = prepDataNear(s_data2,samplesize, nuc)
    x_target, test_x, train_yo, test_yo, num_classes = prepData(s_data1,samplesize, nuc)

    ref_samples = np.arange(x_ref.shape[0])
    print(ref_samples)

    loss, loss_c = [], []
    epochs = []
    print("training...")

    for epochnumber in range(epoch_num):

        x_r, y_r, lc, ld = [], [], [], []

        np.random.shuffle(x_target)
        np.random.shuffle(ref_samples)

        for i in range(len(x_target)):
            try:
                x_r.append(x_ref[ref_samples[i]])
                y_r.append(y_ref[ref_samples[i]])
            except IndexError:
                print('Index Error')

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
    #s_out = "/data/nanopore/IVT/weight/"
    s_out = "/data/nanopore/IVT/weight_dec/"
    bestwight = s_out + "weightwn_keras_dec2.hdf"
    epoch_num = 3
    #train(sys.argv[2], sys.argv[3], sys.argv[1], 10)


    nucs = ('A', 'T', 'C', 'G')
    cnt = 0
    for n1, n2, n3, n4, n5 ,n6 in itertools.product(nucs, nucs, nucs, nucs, nucs, nucs):

        s_data1 = "/data/nanopore/nanoDoc2_1/10000signal"
        s_data2 = "/data/nanopore/nanoDoc2_1/50signal"
        s_out = "/data/nanopore/nanoDoc2_1/weight/docweight"
        nuc = n1+n2+n3+n4+n5+n6
        epoch_num = 5
        samplesize = 12750
        dir = s_out + "/" + nuc
        if os.path.exists(dir):
            continue

        train(s_data1,s_data2, s_out, nuc, bestwight, samplesize, epoch_num)
        cnt+=1