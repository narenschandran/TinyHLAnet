from hashlib import shake_256, md5

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score

import tensorflow as tf
import keras

def model_load(modfn, model_f):
    mod = tf.keras.models.load_model(model_f, custom_objects={modfn.__name__: modfn})
    return mod

def md5sum_compute(f):
    fcon = open(f, 'rb')
    chksum = md5(fcon.read()).hexdigest()
    fcon.close()
    return chksum

def make_uuid(s):
    """
    This function is used to generate a unique identification code
    for each model instance based on the hyperparameters used for the model.
    """
    if not isinstance(s, str):
        raise ValueError("Wrong type")
    return shake_256(s.encode('utf-8')).hexdigest(4)

def writeDict(dct, f):
    """
    Write out dictionaries into plain text format
    """
    with open(f, 'w') as fcon:
        for k, v in dct.items():
            fcon.write(f"{k}\t{v}\n")
    return f

def readDict(f):
    """
    Read dictionaries written out by writeDict
    """
    dct = {}
    with open(f, 'r') as fcon:
        for line in fcon:
            k, v = line.strip().split('\t')
            dct[k] = v
    return dct

def layernames(tfobj):
    """
    Get names of layers in a given object
    """
    return [obj.name for obj in tfobj.layers]


def _layer_get(tfobj, name):
    """
    Subset a single layer by name
    """
    ind = layernames(tfobj).index(name)
    return tfobj.layers[ind]


def layer_get(tfobj, names):
    """
    Subset one or more layers by name. If a single layer is requested (if
    'names' is a single input), then provides the layer directly. Otherwise,
    provides the layers in dictionary format.
    """
    strtype = isinstance(names, str)
    if strtype:
        return _layer_get(tfobj, names)
    else:
        return [_layer_get(tfobj, name) for name in names]


#------------------------------------------------------------------------------#
#                            Performance assessment                            #
#------------------------------------------------------------------------------#
def subset_available(ac, pr, ex_val):
    """
    Given actual values (ac) and predicted values (pr), it subsets only
    cases where the actual value is known (i.e. if the actual value is not
    equal to the exclusion value (ex)).
    """
    def to1d(x):
        if x.ndim > 2:
            raise ValueError("Too many dimensions")
        elif x.ndim == 2:
            return np.squeeze(x)
        elif x.ndim == 1:
            return x
        else:
            raise Value("Catastrophic error")

    pr1 = to1d(pr)
    ac1 = to1d(ac)
    l = ac1 != ex_val
    if not np.any(l):
        return [None, None]
    else:
        return [ac1[l], pr1[l]]


def MSE(ac, pr):
    """
    Convenient function that computes mean-squared error only for datapoints
    where the true value is known.
    """
    ac1, pr1 = subset_available(ac, pr, -1.0)
    if ac1 is None:
        return np.nan
    return np.round(np.sqrt(np.mean(np.square(pr1 - ac1))), 3)


def BCE(ac, pr):
    """
    Convenient function that computes binary crossentropy only for datapoints
    where the true value is known.
    """

    ac1, pr1 = subset_available(ac, pr, -1.0)
    if ac1 is None:
        return np.nan
    pr1 = np.clip(pr1, 0.0001, 0.9999)  # Prevent issues with log and zero values
    tmp0 = -1.0 * ((ac1 * tf.math.log(pr1)) + ((1 - ac1) * tf.math.log(1 - pr1)))
    return np.round(np.mean(tmp0), 3)


def AUC(ac, pr):
    """
    Convenient function that computes AUC only for datapoints where the true
    value is known.
    """
    ac1, pr1 = subset_available(ac, pr, -1.0)
    if ac1 is None:
        return np.nan
    val = roc_auc_score(ac1, pr1)
    return np.round(val, 3)


def PRAUC(ac, pr):
    """
    Convenient function that computes PRAUC only for datapoints where the
    true value is known.
    """
    ac1, pr1 = subset_available(ac, pr, -1.0)
    if ac1 is None:
        return np.nan
    val = average_precision_score(ac1, pr1)
    return np.round(val, 3)


def SRCC(ac, pr):
    """
    Convenient function that computes Spearmen Correlation Coefficient only for
    datapoints the true value is known.
    """
    ac1, pr1 = subset_available(ac, pr, -1.0)
    if ac1 is None:
        return np.nan
    return np.round(spearmanr(ac1, pr1).correlation, 3)

