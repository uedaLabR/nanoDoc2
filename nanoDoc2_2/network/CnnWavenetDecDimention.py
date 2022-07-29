from funcy import rcompose
from keras.models import Model
from keras.layers import *
from keras.regularizers import l2
from keras.engine.base_layer import Layer
from keras import backend as K

from nanoDoc2_2.network import CnnWavenetKeras


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


def cbr(x, out_layer, kernel, stride, dilation):
    x = Conv1D(out_layer, kernel_size=kernel, dilation_rate=dilation, strides=stride, padding="same")(x)
    x = BatchNormalization()(x)
    x = Activation("relu")(x)
    return x

def se_block(x_in, layer_n):
    x = GlobalAveragePooling1D()(x_in)
    x = Dense(layer_n//8, activation="relu")(x)
    x = Dense(layer_n, activation="sigmoid")(x)
    x_out=Multiply()([x_in, x])
    return x_out

def resblock(x_in, layer_n, kernel, dilation, use_se=True):
    x = cbr(x_in, layer_n, kernel, 1, dilation)
    x = cbr(x, layer_n, kernel, 1, dilation)
    if use_se:
        x = se_block(x, layer_n)
    x = Add()([x_in, x])
    return x


def build_network(shape, num_classes,inweight=None,do_r = 0.2):

    model = CnnWavenetKeras.build_network(shape=shape, num_classes=num_classes)
    if inweight != None:
        model.load_weights(inweight)

    do_r = 0.3
    x = model.layers[-3].output
    #dec filter
    num_filters_ = [128,84,56,36,24]
    for n in range(5):
        x = Conv1D(num_filters_[n], 1, padding='same')(x)
        b4 = Dropout(do_r)(x)
        x = convBlockEqual(num_filters_[n], 1, num_filters_[n], 3, num_filters_[n], 1, do_r)(b4)  #
        x = Add()([b4, x])
        x = se_block(x,num_filters_[n])

    x = BatchNormalization()(x)
    x = Activation("relu")(x)  # 4
    x = Flatten()(x)
    x = Dense(1024)(x)
    x = BatchNormalization()(x)
    x = Activation("relu")(x)
    x = Dropout(do_r)(x)
    x = Dense(512)(x)
    x = BatchNormalization()(x)
    x = Activation("relu")(x)
    x = Dropout(do_r)(x)
    x = Dense(256)(x)
    x = Dense(192)(x)
    x = BatchNormalization()(x)
    x = Activation("relu")(x)
    x = Dense(96)(x)
    x = Dense(48)(x)
    x = BatchNormalization()(x)
    x = Activation("relu")(x)
    x = Dense(24)(x)
    x = BatchNormalization()(x)
    x = Activation("relu")(x)
    x = Dense(20)(x)
    x = Activation("tanh", name="latent-features")(x)
    x = Dense(24)(x)
    x = Dense(48)(x)
    x = Dense(96)(x)
    x = BatchNormalization()(x)
    x = Activation("relu")(x)
    x = Dense(192)(x)
    x = Dense(256)(x)
    x = BatchNormalization()(x)
    x = Activation("relu")(x)
    l_output = Dense(num_classes, activation='softmax')(x)
    model_n = Model(inputs=model.inputs, outputs=[l_output])

    for layer in model_n.layers:
        if layer.name == "add_24":
            break
        else:
            layer.trainable = False

    return model_n
