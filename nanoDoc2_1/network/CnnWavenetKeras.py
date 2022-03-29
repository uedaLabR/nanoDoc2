from funcy import rcompose
from keras.models import Model
from keras.layers import *
from keras.regularizers import l2
from keras.engine.base_layer import Layer
from keras import backend as K
from funcy import concat, identity, juxt, partial, rcompose, repeat, take



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

def se_block(x_in, layer_n):
    x = GlobalAveragePooling1D()(x_in)
    x = Dense(layer_n//8, activation="relu")(x)
    x = Dense(layer_n, activation="sigmoid")(x)
    x_out=Multiply()([x_in, x])
    return x_out

def build_network(shape, num_classes,do_r = 0.2):


    def conv1D_halve(filters, kernel_size):
        return Conv1D(filters, kernel_size, padding='same', strides=2, kernel_initializer='he_normal',
                      kernel_regularizer=l2(0.0001))


    num_filters_ = 12
    kernel_size_ = 3

    stacked_layers_ = [12, 8, 4, 1]

    l_input = Input(batch_shape=shape)
    x = GaussianNoise(stddev=0.002)(l_input)
    x = conv1D_halve(24, 3)(x) #
    x = se_block(x, 24)
    x = BatchNormalization()(x)
    b4 = Dropout(do_r)(x)
    x = convBlockEqual(12, 1, 24, 3, 12, 1, do_r)(b4) #
    x = Concatenate()([b4,x])
    x = conv1D_halve(12, 1)(x)  #
    x = se_block(x,12)

    x = WaveNetResidualConv1D(num_filters_, kernel_size_, stacked_layers_[0])(x)
    x = Conv1D(num_filters_*2, 1, padding='same')(x)
    x = WaveNetResidualConv1D(num_filters_*2, kernel_size_, stacked_layers_[1])(x)
    x = Conv1D(num_filters_*4, 1, padding='same')(x)
    x = WaveNetResidualConv1D(num_filters_*4, kernel_size_, stacked_layers_[2])(x)
    x = Conv1D(num_filters_*8, 1, padding='same')(x)
    x = WaveNetResidualConv1D(num_filters_*8, kernel_size_, stacked_layers_[3])(x)

    #l_output = Dense(num_classes, activation='softmax')(x)
    #add
    x = Conv1D(num_filters_ * 16, 1, padding='same')(x)
    x = se_block(x, num_filters_ * 16)
    x = GlobalAveragePooling1D()(x)
    l_output = Dense(num_classes, activation='softmax')(x)


    model = Model(inputs=[l_input], outputs=[l_output])

    return model