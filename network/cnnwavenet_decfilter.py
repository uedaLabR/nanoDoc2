from funcy import rcompose
from keras.models import Model
from keras.layers import *
from keras.regularizers import l2
from keras.engine.base_layer import Layer
from keras import backend as K

from src.nanoDoc2 import cnnwavenet_keras


class Mish(Layer):
    def __init__(self, **kwargs):
        super(Mish, self).__init__(**kwargs)
        self.supports_masking = True

    def call(self, inputs):
        return inputs * K.tanh(K.softplus(inputs))

    def get_config(self):
        base_config = super(Mish, self).get_config()
        return dict(list(base_config.items()))

    def compute_output_shape(self, input_shape):
        return input_shape



def conv1D(filters, kernel_size):
    return Conv1D(filters, kernel_size, padding='same', kernel_initializer='he_normal',
                  kernel_regularizer=l2(0.0001))


def conv1D_halve(filters, kernel_size):
    return Conv1D(filters, kernel_size, padding='same', strides=2, kernel_initializer='he_normal',
                  kernel_regularizer=l2(0.0001))


def convBlock(f1, k1, f2, k2, f3, k3, do_r):
    return rcompose(conv1D(f1, k1),
                    conv1D(f2, k2),
                    conv1D(f3, k3),
                    MaxPooling1D(pool_size=2),
                    BatchNormalization(),
                    Dropout(do_r))

def convBlockEqual(f1, k1, f2, k2, f3, k3, do_r):

    return rcompose(conv1D(f1, k1),
                    conv1D(f2, k2),
                    conv1D(f3, k3),
                    BatchNormalization(),
                    Dropout(do_r))


def WaveNetResidualConv1D(num_filters, kernel_size, stacked_layer):

    def build_residual_block(l_input):

        actv = Mish()
        resid_input = l_input
        for dilation_rate in [2**i for i in range(stacked_layer)]:
            l_sigmoid_conv1d = Conv1D(
              num_filters, kernel_size, dilation_rate=dilation_rate,
              padding='same', activation='sigmoid')(l_input)
            l_tanh_conv1d = Conv1D(
             num_filters, kernel_size, dilation_rate=dilation_rate,
             padding='same', activation=actv)(l_input)
            l_input = Multiply()([l_sigmoid_conv1d, l_tanh_conv1d])
            l_input = Conv1D(num_filters, 1, padding='same')(l_input)
            resid_input = Add()([resid_input ,l_input])
        return resid_input
    return build_residual_block


def build_network(shape, num_classes,do_r = 0.2):

    model = cnnwavenet_keras.build_network(shape=shape, num_classes=num_classes)


    do_r = 0.2
    b4 = Dropout(do_r)(model.layers[-4].output)

    num_filters_ = [80, 60, 40, 20]
    for n in range(4):
        if n > 0:
            b4 = Dropout(do_r)(x)
        x = convBlockEqual(num_filters_[n], 1, num_filters_[n], 3, num_filters_[n], 1, do_r)(b4)  #
        x = Concatenate()([b4, x])
        x = MaxPooling1D(pool_size=2)(x)

    x = Conv1D(20, 1, padding='same')(x)
    x = GlobalAveragePooling1D()(x)
    x = Dense(256, activation='relu')(x)
    l_output = Dense(num_classes, activation='softmax')(x)
    model_n = Model(inputs=model.inputs, outputs=[l_output])
    return model_n