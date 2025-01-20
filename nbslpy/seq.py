from nbslpy._utils import SingleVectorize
from nbslpy.blosum import BLOSUM50, BLOSUM60, BLOSUM80
from nbslpy.aacodes import _aa1to3dct, _aa1lst
from nbslpy.aacodes import aaind, aaohe, aablosum, aaradii_pmhc_ca
import random
import numpy as np


def seq2ind(seq):
    def _seq2ind(seq):
        if not isinstance(seq, str):
            raise ValueError("Input must be a sequence string")
        seqlst = list(seq)
        return aaind(seqlst)

    return np.asarray(SingleVectorize(_seq2ind, seq))


def seq2ohe(seq, dim=21, flatten=False):
    def _seq2ohe(seq, dim, flatten):
        if not isinstance(seq, str):
            raise ValueError("Input must be a sequence string")
        seqlst = list(seq)
        return aaohe(seqlst, dim=dim, flatten=flatten)

    return np.asarray(SingleVectorize(_seq2ohe, seq, dim=dim, flatten=flatten))


def seq2blosum(seq, blosumdict, flatten=False):
    def _seq2blosum(seq, blosumdict, flatten):
        if not isinstance(seq, str):
            raise ValueError("Input must be a sequence string")
        seqlst = list(seq)
        return aablosum(seqlst, blosumdict, flatten)

    return np.asarray(
        SingleVectorize(_seq2blosum, seq, blosumdict=blosumdict, flatten=flatten)
    )


def seq2blosum50(seq, flatten=False):
    return seq2blosum(seq, blosumdict=BLOSUM50, flatten=flatten)


def seq2blosum60(seq, flatten=False):
    return seq2blosum(seq, blosumdict=BLOSUM60, flatten=flatten)


def seq2blosum80(seq, flatten=False):
    return seq2blosum(seq, blosumdict=BLOSUM80, flatten=flatten)

def seq2oheblosum50(seqs, flatten=False):
    if isinstance(seqs, str):
        seqs = [seqs]
    ohe = seq2ohe(seqs, flatten=False)
    bls50 = seq2blosum50(seqs, flatten=False)
    mat = np.concatenate([ohe, bls50], -1)
    if flatten:
        flatmat = np.asarray(
            [
                list(itertools.chain.from_iterable(mat[i, :, :]))
                for i in range(mat.shape[0])
            ]
        )
        # If this isn't done, this screws up the model later one because an
        # unnecessary dimension is added for single sequences, and this affects
        # how our processing happens later on.
        if flatmat.shape[0] == 1:
            flatmat = np.squeeze(flatmat)
        return flatmat
    else:
        if mat.shape[0] == 1:
            mat = np.squeeze(mat)
        return mat


def seq_sanitize(seq):
    def _seq_sanitize(seq):
        return "".join([s if s in _aa1to3dct else "X" for s in list(seq)])

    return SingleVectorize(_seq_sanitize, seq)


def seq_posmutate(seq, pos, seed=1):
    random.seed(1)
    # We don't want X in the final list
    aa1lst = _aa1lst.copy()
    aa1lst.remove("X")

    def _seq_posmutate(seq, pos, lst):
        if not isinstance(seq, str):
            raise TypeError("Input must be a string")
        seqlen = len(seq)
        poschange = pos - 1
        if seqlen == 0:
            raise ValueError("Empty sequences are not accepted")
        elif (pos < 1) or (pos > seqlen):
            raise ValueError(
                "Input sequence %s is %d characters long. Accepted pos input is range [1, %d]"
                % (seq, seqlen, seqlen)
            )
        newseq = list(seq)
        newseq[poschange] = random.choice(lst)
        return "".join(newseq)

    return SingleVectorize(_seq_posmutate, seq, pos=pos, lst=aa1lst)


def seq_posunknown(seq, pos, seed=1):
    def _seq_posmutate(seq, pos):
        if not isinstance(seq, str):
            raise TypeError("Input must be a string")
        seqlen = len(seq)
        poschange = pos - 1
        if seqlen == 0:
            raise ValueError("Empty sequences are not accepted")
        elif (pos < 0) or (pos > seqlen):
            raise ValueError(
                "Input sequence %s is %d characters long. Accepted pos input is range [1, %d]"
                % (seq, seqlen, seqlen)
            )
        newseq = list(seq)
        newseq[poschange] = "X"
        return "".join(newseq)

    return SingleVectorize(_seq_posmutate, seq, pos=pos)


def seq2aaradii_pmhc_ca(seq):
    def _seq2aaradii_pmhc_ca(seq):
        if not isinstance(seq, str):
            raise ValueError("Input must be a sequence string")
        seqlst = list(seq)
        return aaradii_pmhc_ca(seqlst)

    return np.asarray(SingleVectorize(_seq2aaradii_pmhc_ca, seq))
