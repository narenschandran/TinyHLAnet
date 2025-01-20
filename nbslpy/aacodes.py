from nbslpy._utils import SingleVectorize, strlencheck
from nbslpy.blosum import BLOSUM50, BLOSUM60, BLOSUM80
import nbslpy.aaradii
from itertools import chain
import random
import numpy as np

# The order of amino acids in _aa1lst will be used to give it a unique index.
# It is organized so that a symbol for unknown/missing amino acids is defined,
# followed by the standard 20 amino acids. This allows for non-standard amino
# acids to be added later without distrubting the original order. The
# unknown/missing symbol is given first so that it is convenient for masking
# with Tensorflow's Embedding layer.
_aa1lst = [
    # 1 unknown
    "X",

    # 20 standard Amino Acids
    "A", "C", "D", "E", "F",
    "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R",
    "S", "T", "V", "W", "Y",
]

_aa1to3dct = {
    "X": "UNK",

    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}


def aa1to3(aa1):
    def _aa1to3(aa1):
        strlencheck(aa1, 1)
        if aa1 not in _aa1to3dct.keys():
            aa1 = "X"
        return _aa1to3dct[aa1]

    return SingleVectorize(_aa1to3, aa1)


missing_aa1 = [aa1 for aa1 in _aa1lst if aa1 not in _aa1to3dct.keys()]

if len(missing_aa1) > 0:
    msg = (
        "The following amino acids do not have a mapping between 1 and 3 letter codes:\n"
        "\n".join(missing_aa1)
    )
    raise ValueError(msg)

_aa3lst = [aa1to3(aa1) for aa1 in _aa1lst]
_aa3to1dct = {_aa1to3dct[aa1]: aa1 for aa1 in _aa1lst}


def aa3to1(aa3):
    def _aa3to1(aa3):
        strlencheck(aa3, 3)
        if aa3 not in _aa3to1dct.keys():
            aa3 = "UNK"
        return _aa3to1dct[aa3]

    return SingleVectorize(_aa3to1, aa3)


def toaa1(aa):
    """
    User facing function for getting single-letter amino acid code
    """
    def _toaa1(aa):
        if len(aa) == 1:
            aa1 = aa if aa in _aa1to3dct.keys() else "X"
        elif len(aa) == 3:
            aa1 = aa3to1(aa)
        else:
            raise ValueError(
                "Input [%s] is neither a one-letter or three-letter code" % (aa)
            )
        return aa1

    return SingleVectorize(_toaa1, aa)


def toaa3(aa):
    """
    User facing function for getting three-letter amino acid code
    """
    return aa1to3(toaa1(aa))


def aaind(aa):
    def _aaind(aa):
        return _aa1lst.index(toaa1(aa))

    return np.asarray(SingleVectorize(_aaind, aa))


def indaa(ind):
    def _indaa(ind):
        return _aa1lst[ind]

    return SingleVectorize(_indaa, ind)


def aaohe(aa, dim=21, flatten=False):
    def _aaohe(aa):
        ind = aaind(aa)
        if not (ind < dim):
            raise ValueError("dim must be at least %d for given input" % (ind + 1))
        ohe = [0] * dim
        ohe[ind] = 1
        return ohe

    ohes = SingleVectorize(_aaohe, aa)
    if flatten:
        ohes = list(chain.from_iterable(ohes))
    return np.asarray(ohes)

def aablosum(aa, blosumdict, flatten=False):
    def _aablosum(aa1, blosumdict):
        aa1 = toaa1(aa1)
        if aa1 not in blosumdict:
            raise ValueError("The amino acid %s is not known to BLOSUM" % (aa1))
        aa1dct = blosumdict[aa1]
        keys = [aa for aa in _aa1lst if aa in aa1dct.keys()]
        return [blosumdict[aa1][key] for key in keys]

    bls = SingleVectorize(_aablosum, aa, blosumdict=blosumdict)
    if flatten:
        bls = list(chain.from_iterable(bls))
    return np.asarray(bls)

def aablosum50(aa, flatten=False):
    return aablosum(aa, BLOSUM50, flatten=flatten)

def aablosum60(aa, flatten=False):
    return aablosum(aa, BLOSUM60, flatten=flatten)

def aablosum80(aa, flatten=False):
    return aablosum(aa, BLOSUM80, flatten=flatten)

def aaradii(aa, radiidict):
    allvals = list(radiidict.values())
    avgval = sum(allvals) / len(allvals)

    def _aaradii(aa1, radiidict):
        aa1 = toaa1(aa1)
        if aa1 not in radiidict:
            print(
                "the amino acid %s is not known to aaradii, returning an average of all known values"
                % (aa1)
            )
            return avgval
        return radiidict[aa1]

    return SingleVectorize(_aaradii, aa, radiidict=radiidict)

def aaradii_pmhc_ca(aa):
    return aaradii(aa, nbslpy.aaradii.pmhc_ca)
