import pandas as pd
import pickle
import re


def SingleVectorize(fn, data, **kwargs):
    if isinstance(data, list) or isinstance(data, pd.core.series.Series):
        return [fn(x, **kwargs) for x in data]
    else:
        return fn(data, **kwargs)


def strlencheck(data, length):
    def _strlencheck(x, length):
        if not isinstance(x, str):
            raise ValueError("Input must be a string")
        if len(x) != length:
            raise ValueError("Input [%s] is not of size [%d]" % (x, length))

    return SingleVectorize(_strlencheck, data, length=length)


def grepl(pattern, strings, flags=0):
    def _grepl(s, pat, fl):
        if not isinstance(s, str):
            raise ValueError("Input must be a string")
        return re.search(pattern=pat, string=s, flags=fl) is not None

    return SingleVectorize(_grepl, strings, pat=pattern, fl=flags)


def which(data):
    if data is None:
        return None

    if isinstance(data, bool):
        data = [data]

    if isinstance(data, list):
        if not all([isinstance(x, bool) for x in data]):
            raise ValueError("The input must only be boolean values")
    return [i for i in range(len(data)) if data[i]]


def savePickle(obj, file):
    fcon = open(file, "wb")
    pickle.dump(obj, fcon)
    fcon.close()
    return file


def readPickle(file):
    fcon = open(file, "rb")
    obj = pickle.load(fcon)
    fcon.close()
    return obj
