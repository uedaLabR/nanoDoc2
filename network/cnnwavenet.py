import tensorflow as tf
from tensorflow.keras.layers import Activation
from tensorflow.keras.utils import get_custom_objects
from funcy import concat, identity, juxt, partial, rcompose, repeat, take

from tensorflow.keras.models import Sequential, Model, save_model
from tensorflow.keras.layers import *
from tensorflow.keras.regularizers import l2


class Mish(Activation):
    '''
    Mish Activation Function.
    .. math::
        mish(x) = x * tanh(softplus(x)) = x * tanh(ln(1 + e^{x}))
    Shape:
        - Input: Arbitrary. Use the keyword argument `input_shape`
        (tuple of integers, does not include the samples axis)
        when using this layer as the first layer in a nnmodels.
        - Output: Same shape as the input.
    Examples:
        >>> X = Activation('Mish', name="conv1_act")(X_input)
    '''

    def __init__(self, activation, **kwargs):
        super(Mish, self).__init__(activation, **kwargs)
        self.__name__ = 'Mish'


def mish(inputs):
    return inputs * tf.math.tanh(tf.math.softplus(inputs))

get_custom_objects().update({'mish': Mish(mish)})


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
        resid_input = l_input
        for dilation_rate in [2**i for i in range(stacked_layer)]:
            l_sigmoid_conv1d = Conv1D(
              num_filters, kernel_size, dilation_rate=dilation_rate,
              padding='same', activation='sigmoid')(l_input)
            l_tanh_conv1d = Conv1D(
             num_filters, kernel_size, dilation_rate=dilation_rate,
             padding='same', activation='mish')(l_input)
            l_input = Multiply()([l_sigmoid_conv1d, l_tanh_conv1d])
            l_input = Conv1D(num_filters, 1, padding='same')(l_input)
            resid_input = Add()([resid_input ,l_input])
        return resid_input
    return build_residual_block


def build_network(shape, num_classes,do_r = 0.2):

    def ljuxt(*fs):
        return rcompose(juxt(*fs), list)

    def conv1D(filters, kernel_size):
        return Conv1D(filters, kernel_size, padding='same', kernel_initializer='he_normal',
                      kernel_regularizer=l2(0.0001))

    def conv1D_halve(filters, kernel_size):
        return Conv1D(filters, kernel_size, padding='same', strides=2, kernel_initializer='he_normal',
                      kernel_regularizer=l2(0.0001))

    def dense(units, activation):
        return Dense(units, activation=activation, kernel_regularizer=l2(0.0001))

        # Define SqueezeNet.

    def fire_module(filters_squeeze, filters_expand):
        return rcompose(BatchNormalization(),
                        Activation('mish'),
                        conv1D(filters_squeeze, 1),
                        BatchNormalization(),
                        Activation('mish'),
                        ljuxt(conv1D(filters_expand // 2, 1),
                              conv1D(filters_expand // 2, 3)),
                        Concatenate())

    def fire_module_with_shortcut(filters_squeeze, filters_expand):
        return rcompose(ljuxt(fire_module(filters_squeeze, filters_expand),
                              identity),
                        Add())

        # InceptionModule.

    def inception():
        u1 = rcompose(AveragePooling1D(pool_size=3, strides=1, padding='same'),
                      conv1D(48, 1))
        u2 = conv1D(48, 1)
        u3 = rcompose(conv1D(16, 1),
                      conv1D(48, 3))
        u4 = rcompose(conv1D(16, 1),
                      conv1D(48, 3),
                      conv1D(48, 3))

        return rcompose(ljuxt(u1, u2, u3, u4),
                        Concatenate(axis=2))

    num_filters_ = 12
    kernel_size_ = 3

    stacked_layers_ = [12, 8, 4, 1]
    #stacked_layers_ = [20,16,12, 8, 4, 1]
    l_input = Input(batch_shape=shape)

    # nnBlock = rcompose(GaussianNoise(stddev=0.02),
    #                    conv1D_halve(16, 3),
    #                    BatchNormalization(),
    #                    Dropout(do_r),
    #                    convBlock(48, 3, 48, 3, 48, 3, do_r),
    #                    convBlock(16, 1, 48, 3, 48, 3, do_r),
    #                    convBlock(16, 1, 48, 3, 48, 3, do_r),
    #                    fire_module(8, 64),
    #                    fire_module_with_shortcut(8, 64),
    #                    fire_module(16, 128),
    #                    Activation('mish'))
    #
    # x = nnBlock(l_input)
    #input = Input(batch_shape=shape)
    x = GaussianNoise(stddev=0.002)(l_input)
    x = conv1D_halve(32, 3)(x) #
    x = BatchNormalization()(x)
    b4 = Dropout(do_r)(x)
    x = convBlockEqual(16, 1, 32, 3, 16, 1, do_r)(b4) #
    x = Concatenate()([b4,x])
    # b4 = MaxPooling1D(pool_size=2)(x)
    # x = convBlockEqual(24, 1, 48, 3, 24, 1, do_r)(b4) #
    # x = Concatenate()([b4,x])
    # b4 = MaxPooling1D(pool_size=2)(x)
    # x = convBlockEqual(40, 1, 80, 3, 40, 1, do_r)(b4) #
    # x = Concatenate()([b4,x])

    x = Conv1D(num_filters_, 1, padding='same')(x)
    x = WaveNetResidualConv1D(num_filters_, kernel_size_, stacked_layers_[0])(x)
    x = Conv1D(num_filters_*2, 1, padding='same')(x)
    x = WaveNetResidualConv1D(num_filters_*2, kernel_size_, stacked_layers_[1])(x)
    x = Conv1D(num_filters_*4, 1, padding='same')(x)
    x = WaveNetResidualConv1D(num_filters_*4, kernel_size_, stacked_layers_[2])(x)
    x = Conv1D(num_filters_*8, 1, padding='same')(x)
    x = WaveNetResidualConv1D(num_filters_*8, kernel_size_, stacked_layers_[3])(x)
    #
    # x = Conv1D(num_filters_*16, 1, padding='same')(x)
    # x = WaveNetResidualConv1D(num_filters_*16, kernel_size_, stacked_layers_[4])(x)
    # x = Conv1D(num_filters_*32, 1, padding='same')(x)
    # x = WaveNetResidualConv1D(num_filters_*32, kernel_size_, stacked_layers_[5])(x)

    #l_output = Dense(num_classes, activation='softmax')(x)
    #add
    x = Conv1D(num_filters_ * 16, 1, padding='same')(x)
    x = GlobalAveragePooling1D()(x)
    l_output = Dense(num_classes, activation='softmax')(x)

    # x = Flatten()(x)
    # x = Conv1D(num_classes, 1)(x)
    # x = GlobalAveragePooling1D()(x)
    # l_output = Softmax()(x)

    model = Model(inputs=[l_input], outputs=[l_output])
    #opt = Adam(lr=LR)
    #opt = tfa.optimizers.SWA(opt)
    #nnmodels.compile(loss=losses.CategoricalCrossentropy(), optimizer=opt, metrics=['accuracy'])
    return model