import os

import pandas as pd

import tensorflow as tf
import keras
from keras.layers import Bidirectional, Dense, Embedding, Flatten, LSTM

from nbslpy.aacodes import _aa1lst as aa1lst

from deephlaffylib.Regularizers import (
    UnitNormMaxRegularizer,
    UnitNormRegularizer,
)

def AminoEmbLayer(embdim, name, reg = 'maxl2norm', **kwargs):
    '''Wrapper function to produce amino acid embedding layer

    Args:
        embdim : Embedding dimension
        name   : Embedding layer name
        reg    : Embedding regularizer. Can be a conventional Keras regularizer
                 or one of 'l2norm' or 'maxl2norm', the latter two of which
                 are defined by our library
    '''
    if reg == 'l2norm':
        embreg = UnitNormRegularizer()
    elif reg == 'maxl2norm':
        embreg = UnitNormMaxRegularizer()
    else:
        embreg = None
    return Embedding(input_dim=len(aa1lst), output_dim=embdim,
                     name=name, embeddings_regularizer=embreg,
                     **kwargs)



#------------------------------------------------------------------------------#
#                                  Operations                                  #
#------------------------------------------------------------------------------#
# The layers defined in this section do standard operations such as summming a #
# matrix, dot product and concatenation. Although these can be defined as      #
# functions, there were certain issues (at least on Tensorflow 2.18) when      #
# trying to contruct a secondary model based on another model's components     #
# wherein outputs from certain model layers could not be directly used         #
# for operations such as dot product.                                          #
#------------------------------------------------------------------------------#
@keras.saving.register_keras_serializable("deephlaffylib")
class DotProduct(keras.layers.Layer):
    '''Used to compute pHLA-I interaction strength'''
    def __init__(self, name):
        super().__init__(name = name)

    def call(self, hemb, pemb):
        pp = hemb @ tf.transpose(pemb, (0, 2, 1))
        return pp

    def get_config(self):
        return {'name' : self.name}

@keras.saving.register_keras_serializable("deephlaffylib")
class MatrixSum(keras.layers.Layer):
    '''Used to sum interaction strength'''
    def __init__(self, name):
        super().__init__(name = name)

    def call(self, inputs):
        return tf.reduce_sum(tf.reduce_sum(inputs, axis=-1), axis=-1, keepdims=True)

    def get_config(self):
        return {'name' : self.name}

@keras.saving.register_keras_serializable("deephlaffylib")
class Int32Cast(keras.layers.Layer):
    '''Cast input into 32-bit Integer'''
    def __init__(self, name):
        super().__init__(name = name)

    def call(self, inputs):
        return tf.cast(inputs, tf.int32)

    def get_config(self):
        return {'name' : self.name}

@keras.saving.register_keras_serializable("deephlaffylib")
class ConcatenateD1(keras.layers.Layer):
    '''Concatenate at the first dimension'''
    def __init__(self, name):
        super().__init__(name = name)

    def call(self, inputs):
        return tf.concat(inputs, axis = 1)

    def get_config(self):
        return {'name' : self.name}


@keras.saving.register_keras_serializable("deephlaffylib")
class SubsetpHLA(keras.layers.Layer):
    '''Concatenate at the first dimension'''
    def __init__(self, hlen, name):
        super().__init__(name = name)
        self.hlen = hlen

    def call(self, inputs):
        return inputs[:, :self.hlen:, self.hlen:]

    def get_config(self):
        return {'name' : self.name, 'hlen' : self.hlen}


#------------------------------------------------------------------------------#

@keras.saving.register_keras_serializable("deephlaffylib")
class ClipValues(tf.keras.layers.Layer):
    '''Layer to clip output values'''
    def __init__(self, name, minval, maxval):
        super().__init__(name=name)
        self.minval = minval
        self.maxval = maxval

    def call(self, inputs):
        return tf.clip_by_value(
            inputs, clip_value_min=self.minval, clip_value_max=self.maxval
        )

    def get_config(self):
        return {
            'name'   : self.name,
            'minval' : self.minval,
            'maxval' : self.maxval
        }


@keras.saving.register_keras_serializable("deephlaffylib")
class Lookup(tf.keras.layers.Layer):
    '''Create a dictionary layer

    Args:
        keys  : List of keys
        values: List of values to correspond with keys
        defval: Default value to present to user if key is not present
    '''
    def __init__(self, keys, values, defval, **kwargs):
        super().__init__(**kwargs)
        self.keys   = keys
        self.values = values
        self.defval = defval
        self.vtable = None

    def call(self, inputs):
        return self.vtable.lookup(inputs)

    def build(self, input_shape):
        self.vtable = tf.lookup.StaticHashTable(
            tf.lookup.KeyValueTensorInitializer(keys=tf.constant(self.keys), values=tf.constant(self.values)),
            self.defval,
        )

    def get_config(self):
        config = {
            "keys"   : self.keys,
            "values" : self.values,
            "defval" : self.defval,
        }
        return config
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                      ANN-related functions and classes                       #
#------------------------------------------------------------------------------#
#  Since our models use ANNs in various places and since hyperparameter        #
#  optimization is an important task during our model tuning, we define some   #
#  functions to facilitate the construction of ANNs.                           #
#------------------------------------------------------------------------------#

def parse_dense_conf(x):
    '''Parse dense layer configuration string into a usable dictionary

    Format of the configuration string:
        <actv><nodes><kern-and-bias-reg:optional>:<activity-reg:optional>

        Example: 'sigmoid16:l2' constructs a Dense layer with 16 nodes
        that uses a sigmoid activation function with an l2-norm activity
        regularizer but no regularization for the kernal and bias.

    Args:
        x: Configuration string

    Returns:
        Dense layer with the required configuration

    '''

    # Parse the activity regularization
    sp = x.split(":")
    if len(sp) == 1:
        areg = None
    elif len(sp) == 2:
        areg = sp[1]
    else:
        raise ValueError("Invalid configuration")

    s = sp[0]
    i = 0
    slen = len(s)

    # Parse until the first number to get the activation function
    while i < slen:
        if not s[i].isalpha():
            break
        i = i + 1
    actv = s[:i]

    # Parse till the end of the number to get the number of nodes
    j = i
    while j < slen:
        if s[j].isalpha():
            break
        j = j + 1
    units = int(s[i:j])

    # If anything else is remaining, it is the regularizer
    reg = s[j:] if j < slen else None

    # Parsing custom regularizers
    if areg == "maxl2norm":
        areg = UnitNormMaxRegularizer()
    if areg == "l2norm":
        areg = UnitNormRegularizer()

    return {
        'units'                : units,
        'activation'           : actv,
        'kernel_regularizer'   : reg,
        'bias_regularizer'     : reg,
        'activity_regularizer' : areg,
    }

def dense_from_conf(x):
    '''Construct Dense layer from a configuration string

    Args:
        x: Configuration string

    Returns:
        Dense layer with the required configuration

    '''
    conf = parse_dense_conf(x)
    return Dense(**conf)

def ANN(conf, name, flatten=False, clip=None):
    '''Construct an ANN from a configuration string

    ANN Conf format:
        <layer1conf>_<layer2conf>..."

    For the layer configuration, refer to the dense_from_conf function
    '''
    mod = tf.keras.Sequential(name=name)
    if flatten:
        mod.add(Flatten())

    conf_sp = conf.split("_")
    for layer_conf in conf_sp:
        mod.add(dense_from_conf(layer_conf))

    if clip is not None:
        lower, upper = clip
        mod.add(ClipValues(f"{name}clip", lower, upper))
    return mod


#------------------------------------------------------------------------------#
#                               Stand-in classes                               #
#------------------------------------------------------------------------------#
# In cases where a model is built without a specific module, the classes       #
# defined in this section help to produce stand-in values that can be used     #
# without having to define a different Model class that accounts for the       #
# absence of the module. The outputs of the classes defined here are designed  #
# to not interfere with computations being done elsewhere in the model.        #
#------------------------------------------------------------------------------#

@keras.saving.register_keras_serializable("deephlaffylib")
class EmitOnesMat(tf.keras.layers.Layer):
    '''Emits a square matrix filled with ones.

    Args:
        seqlen: Used to set the size of the square matrix sides.
        name  : Name of the layer

    Returns:
        Returns an object that emits matrices of dimension [1, seqlen, seqlen]

    Notes:
        This will be used in DeepHLAffy for models that consider all possible
        contacts ("allpairs" contacts type).
    '''
    def __init__(self, seqlen, name):
        super().__init__(name = name)
        self.seqlen = seqlen
        self.mat = tf.ones([1, self.seqlen, self.seqlen], dtype = tf.float32)

    def call(self, inputs):
        hlainp, pepinp = inputs
        return self.mat

    def get_config(self):
        return {'seqlen' : self.seqlen, 'name' : self.name}


@keras.saving.register_keras_serializable("deephlaffylib")
class EmitZeros(tf.keras.layers.Layer):
    '''Emits a matrix filled with zeros, with the shape set by inputs.'''
    def __init__(self, name):
        super().__init__(name = name)

    def call(self, inputs):
        inp = inputs
        zerovec = tf.expand_dims(tf.zeros_like(inp, dtype = tf.float32), -1)
        return zerovec[:,0,:]

    def get_config(self):
        return {'name': self.name}



#------------------------------------------------------------------------------#
#                               Contacts module                                #
#------------------------------------------------------------------------------#

# This snippet is used to conveniently load the Cα contact distance matrix.
libdir     = os.path.dirname(__file__)
projroot   = os.path.join(libdir, '..')
prereq_dir = os.path.join(projroot, 'prereq')
cmat_file  = os.path.join(prereq_dir, "contact-calpha.tsv.xz")
cmat_df    = pd.read_csv(cmat_file, sep="\t")
cmat       = cmat_df.values

@keras.saving.register_keras_serializable("deephlaffylib")
def RadiiMatrix(rad1, rad2):
    '''
    Given two vectors of Cα radii, it produces a matrix composed of all pairs
    of sums constructed from these vectors.
    '''
    len1 = rad1.shape[1]
    len2 = rad2.shape[1]
    tmp1 = tf.tile(rad1[:, :, tf.newaxis], [1, 1, len2])
    tmp2 = tf.transpose(tf.tile(rad2[:, :, tf.newaxis], [1, 1, len1]), (0, 2, 1))
    radmat = tmp1 + tmp2
    return radmat


@keras.saving.register_keras_serializable("deephlaffylib")
class BinaryContacts(tf.keras.layers.Layer):
    '''
    This layer is used to produce a mask to nullify interactions occuring at
    infeasible locations. The feasible contacts are determined during layer
    instantiation when a Cα distance matrix is given as an input,
    and no further calcultions are done during the layer calling.
    '''
    def __init__(self, cmat, name):
        super().__init__(name=name)
        self.cmat = cmat
        self.contact_mask = tf.expand_dims(tf.cast(self.cmat > 0, tf.float32), 0)

    def call(self, inputs):
        _ = inputs
        return self.contact_mask

    def get_config(self):
        conf = {'name' : self.name, 'cmat' : self.cmat}


@keras.saving.register_keras_serializable("deephlaffylib")
class ExpDistanceDecay(tf.keras.layers.Layer):
    '''
    This layer is used to produce a mask to nullify interactions occuring at
    infeasible locations. The feasible contacts are determined using a
    Cα distance matrix (embedded during layer initialization) and the specific
    Cα radii for the resuidues in the HLA pocket and peptide (given as input
    during layer call). Contacts that are not feasible are nullified.
    '''
    def __init__(self, cmat, name):
        super().__init__(name=name)
        self.cmat = cmat
        self.scale = None
        self.bias = None

    def build(self, input_shape):
        self.bias = self.add_weight(
            shape=(1,),
            initializer=tf.keras.initializers.RandomUniform(minval=4, maxval=8),
            constraint=tf.keras.constraints.NonNeg(),
            trainable=True,
            name="Bias",
        )

        self.scale = self.add_weight(
            shape=(1,),
            initializer=tf.keras.initializers.RandomUniform(minval=1, maxval=4),
            constraint=tf.keras.constraints.MinMaxNorm(1.0, 10.0),
            trainable=True,
            name="Scale",
        )

    def call(self, inputs):
        hlarad, peprad = inputs
        seqrad = tf.concat([hlarad, peprad], axis = 1)
        radmat = RadiiMatrix(seqrad, seqrad)
        radmat_biased = radmat + self.bias
        difmat = self.cmat - radmat_biased

        # This ensures that there aren't values that are too large or too small
        # when doing the exponential.
        difmat_scaled = tf.clip_by_value(self.scale * difmat, -20.0, 20.0)
        omat = tf.ones(difmat_scaled.shape[1:])
        divmat = tf.math.divide_no_nan(omat, tf.exp(difmat_scaled))
        resmat0 = tf.clip_by_value(divmat, 0.0, 1.0)
        resmat = resmat0 * tf.cast(self.cmat > 0, tf.float32)
        return resmat

    def get_config(self):
        conf = {'name': self.name, 'cmat' : self.cmat}
        return conf

def ContactsLayer(contacts_type, name = 'contacts', seqlen = None):
    """
    Wrapper function to instantiate the appropriate contact layer
    """
    if (contacts_type is None) or (contacts_type == "allpairs"):
        if seqlen is None:
            raise ValuerError("Require seqlen for allpairs contacts")
        ones_mat = tf.ones([seqlen, seqlen])
        contactsl = BinaryContacts(ones_mat, name="contacts")
    elif contacts_type == "simple":
        contactsl = BinaryContacts(cmat, name="contacts")
    elif contacts_type == "expdecay":
        contactsl = ExpDistanceDecay(cmat, name="contacts")
    else:
        raise ValueError("Unknown contacts type")

    return contactsl

#------------------------------------------------------------------------------#
#                         Positional importance module                         #
#------------------------------------------------------------------------------#
@keras.saving.register_keras_serializable("deephlaffylib")
class PosImpBlock(tf.keras.layers.Layer):
    '''Positional importance module

    Args:
        embdim  : Embedding size for amino acids in this module
        annconf : Configuration used to instantiate the ANN that will produce
                  the positional importance vectors. The final layer should
                  have size equal to hlalen + peplen.
    '''
    def __init__(self, embdim, annconf, name):
        super().__init__(name=name)

        self.embdim  = int(embdim)
        self.annconf = annconf

        self.hla_embmaker = AminoEmbLayer(self.embdim, name = f"{name}_hlaposemb")
        self.pep_embmaker = AminoEmbLayer(self.embdim, name = f"{name}_pepposemb")
        self.ann = ANN(self.annconf, name = f"{name}_ann", flatten = True)

        self.dotp = DotProduct(f"{name}_dotp")

    def call(self, inputs):

        hlainp, pepinp = inputs

        # Generate HLA-I and peptide embeddings
        hlaemb, pepemb = self.hla_embmaker(hlainp), self.pep_embmaker(pepinp)
        posemb = tf.concat([hlaemb, pepemb], axis = 1)

        # Calculate residue-wise importance
        posvec = tf.expand_dims(self.ann(posemb), axis = -1)

        # Compute importance of each residue-pair interaction
        posmat = self.dotp(posvec, posvec)
        return posmat

    def config(self):
        config = {
            'embdim'  : self.embdim,
            'annconf' : self.annconf,
            'name'    : self.name,
        }
        return config


#------------------------------------------------------------------------------#
#                                Miscellaneous                                 #
#------------------------------------------------------------------------------#
def default_output_layers():
    '''Default output layers for all DeepHLAffy models'''
    regl = tf.keras.Sequential(
        [Dense(1, activation="relu", name="regdense"), ClipValues("regclip", 0.0, 1.0)],
        name="reg",
    )
    bindl = tf.keras.Sequential(
        [
            Dense(1, activation="sigmoid", name="binddense"),
            ClipValues(f"bindclip", 0.0, 1.0),  # Redundant
        ],
        name="bind",
    )
    return regl, bindl

