import os
from copy import deepcopy
import re
import numpy as np
import pandas as pd

import tensorflow as tf
import keras

from nbslpy.aacodes import _aa1lst
standard_aminos = deepcopy(_aa1lst)
if "X" in standard_aminos:
    standard_aminos.remove("X")


from nbslpy.seq import seq2ind, seq2aaradii_pmhc_ca
from nbslpy._utils import readPickle
from nbslpy.pam import PAM70 as PAM

from deephlaffylib.DeepHLAffy import DeepHLAffy as Mod
from deephlaffylib.Utils import *

libdir     = os.path.dirname(__file__)
projroot   = os.path.join(libdir, '..')
prereq_dir = os.path.join(projroot, 'prereq')
mdl_dir    = os.path.join(projroot, 'models')


mdl_f      = os.path.join(mdl_dir, 'deephlaffy', 'both',
                          'deephlaffy-simple-posmodel-effects-1fc66164',
                          'seed-02', 'model.keras')


#------------------------------------------------------------------------------#
#                              Precomputed inputs                              #
#------------------------------------------------------------------------------#
# The following are pre-computed files which are loaded if available.          #
#------------------------------------------------------------------------------#
hlainp_f   = os.path.join(projroot, 'run-data', 'hlainp.pkl')
triplets_f = os.path.join(projroot, 'run-data', 'triplets.pkl')

if os.path.exists(hlainp_f):
    hlainp_dct = readPickle(hlainp_f)

if os.path.exists(triplets_f):
    triplets_dct = readPickle(triplets_f)

pkts_file = os.path.join(prereq_dir, "hlapockets.tsv.xz")
hla_pkts  = pd.read_csv(pkts_file, sep = '\t')
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#                       Input processing and validation                        #
#------------------------------------------------------------------------------#


#------------------------------------------------#
#                Input processing                #
#------------------------------------------------#

def seq_process(s):
    '''Convert sequence to DeepHLAffy input'''
    indrepr = np.expand_dims(seq2ind(s), -1)
    radrepr = np.expand_dims(seq2aaradii_pmhc_ca(s), -1)
    reprvec = np.concatenate([indrepr, radrepr], -1)
    return reprvec

def hla_inputs(x, dct = None):
    # Coerce to list if numpy or string
    if isinstance(x, np.ndarray):
        x = list(x)
    if isinstance(x, str):
        x = [x]

    # If the data is not of type list at this point, assume that
    # it is a DataFrame.
    if isinstance(x, list):
        hlas = x
    else:
        hlas = x.hla

    if dct is None:
        dct = hlainp_dct
    hlainp = np.asarray([dct[hla] for hla in hlas])
    return hlainp

def pep_inputs(x):
    # Coerce to list if numpy or string
    if isinstance(x, np.ndarray):
        x = list(x)
    if isinstance(x, str):
        x = [x]

    # If the data is not of type list at this point, assume that
    # it is a DataFrame.
    if isinstance(x, list):
        peptides = x
    else:
        peptides = x.peptide

    t1 = np.asarray([triplets_dct[pep[:3]]  for pep in peptides])
    t2 = np.asarray([triplets_dct[pep[3:6]] for pep in peptides])
    t3 = np.asarray([triplets_dct[pep[6:9]] for pep in peptides])
    pepinp  = np.concatenate([t1, t2, t3], axis = 1)
    return pepinp


#------------------------------------------------#
#               Allele validation                #
#------------------------------------------------#
def proc_allele(al):
    ''' Process allele name into the format used by our program

    Args:
        al - Single HLA-I allele ID

    Returns:
        Either the HLA-I allele ID as accpeted by DeepHLAffy, or "invalid"
        if the input ID is not recognized by the function.

    Note:
        We expect the allele names to be in the format HLA-[ABC]*<num1>:<num2>.
        We deal with cases where extra resolution is provided and in case the
        'HLA-' prefix is not provided.
    '''
    al_sp1 = al.split(':')
    if len(al_sp1) < 2:
        return "Invalid"
    al2 = ':'.join(al_sp1[:2])
    if "*" not in al2:
        return "Invalid"
    if al2[:4] != "HLA-":
        al2 = "HLA-" + al2
    return al2

def _validate_alleles(als):
    ''' Helper function to parse and report bad allele inputs

    Args:
        als - List containing allele IDs

    Returns:
        out - Input list whose elements are the full allele ID corresponding
              to the input wherever possible. If unable to parse the input
              allele into the full allele, it will be annotated as "Invalid"
        dct - The dictionary containing input to output maping
        inv - List of alleles IDs marked as invalid
    '''
    dct = {}
    out = []
    inv = []
    for al in als:
        if al not in dct:
            dct[al] = proc_allele(al)
            if dct[al] == "Invalid":
                inv.append(al)
        out.append(dct[al])
    return out, dct, inv

def validate_alleles(als, dct = None):
    ''' User facing function for validating input alleles

    Args:
        als - List of alleles to be validated
        dct - Optional dictionary that can be used to check if
              the input alleles have preprocessed data in a dictionary.
              Typically used to check if precomputed inputs are available
              for the validated alleles.
    '''
    proc_als, _, inv = _validate_alleles(als)
    if len(inv) > 0:
        inv_al = '\n'.join(inv)
        raise ValueError(f"The following alleles are invalid:\n{inv_al}")

    if dct is not None:
        hset = set(proc_als)
        missing_h = [h for h in hset if h not in dct.keys()]
        if len(missing_h) > 0:
            missing_str = '\n'.join(missing_h)
            raise ValueError(f"The following alleles are not known:\n{missing_str}")

    return proc_als

#------------------------------------------------#
#               Allele validation                #
#------------------------------------------------#
def validate_peptides(peptides):
    ''' Wrapper function to check if peptides are valid

    This function ensures that the input peptides are 9-mers
    '''
    proc_pep = [pep.strip() for pep in peptides]
    invalid_peps = [pep for pep in proc_pep if len(pep) != 9]
    if len(invalid_peps) > 0:
        inv_pep = '\n'.join(invalid_peps)
        raise ValueError(f"The following peptides are not 9-mers:\n{inv_pep}")
    return proc_pep


#------------------------------------------------#
#             Default model helpers              #
#------------------------------------------------#
def load_and_check_input(f):
    ''' Wrapper function to check whether an input DataFrame can be used as a
        valid input to the default DeepHLAffy model '''
    x = pd.read_csv(f, sep = '\t')
    return check_input(x)

def check_input(x):
    ''' Wrapper function to check whether an input DataFrame can be used as a
        valid input to the default DeepHLAffy model '''

    if ('hla' in x.columns) and ('allele' not in x.columns):
        cls = []
        for col in list(x.columns):
            if col == 'hla':
                cls.append('allele')
            else:
                cls.append(col)
        x.columns = cls

    # Ensure that the necessary columns are present to compute inputs
    if 'allele' not in x.columns:
        raise ValueError("Input needs to have an 'allele' column")

    if 'peptide' not in x.columns:
        raise ValueError("Input needs to have a 'peptide' column")

    # Ensure that the input HLA alleles are valid and that all alleles have
    # pockets associated with them.
    x.loc[:, "hla_map"] = validate_alleles(x.allele.values.tolist(), hlainp_dct)
    x.loc[:,"proc_pep"] = validate_peptides(x.loc[:,"peptide"].values.tolist())

    cols = ["hla_map", "proc_pep"]
    other_cols = ['regressand', 'binder', 'set', 'label', 'pam', 'orig_peptide', 'mut_peptide', 'mutcode']
    for col in other_cols:
        if col in x.columns:
            cols.append(col)
    datf = x.loc[:,cols]
    cnames = list(datf.columns)
    cnames[0] = "hla"
    cnames[1] = "peptide"
    datf.columns = cnames
    return datf

def ModelLoad():
    mod = model_load(Mod, mdl_f)
    return mod

def configure_model(conf = 'default', precomputed_pepeffects = False):
    '''Reconstruct required modules from the full model'''
    mod  = ModelLoad()

    # Operation layers
    i32_1       = layer_get(mod, 'i32_1')
    i32_2       = layer_get(mod, 'i32_2')
    dp          = layer_get(mod, 'dotp')
    msum        = layer_get(mod, 'matsum')
    concat      = layer_get(mod, 'concat')
    subset_phla = layer_get(mod, 'subset_phla')

    hlen = 44
    plen = 9
    hlainpfull = keras.layers.Input((hlen, 2), name = 'hlainp')
    pepinpfull = keras.layers.Input((plen, 2), name = 'pepinp')
    hlainp, pepinp = i32_1(hlainpfull[:,:,0]), i32_2(pepinpfull[:,:,0])
    hlarad, peprad = hlainpfull[:,:,1], pepinpfull[:,:,1]

    hlaembl    = layer_get(mod, "hla")
    pepembl    = layer_get(mod, "pep")
    contactsl  = layer_get(mod, "contacts")
    posimpl    = layer_get(mod, "posimp")
    pepeffectl = layer_get(mod, "pepeffect")
    regl       = layer_get(mod, "reg")
    bindl      = layer_get(mod, "bind")

    # [1]: Interaction strength between different residue pairs
    #      contingent on no other factor.
    hlaemb = hlaembl(hlainp)
    pepemb = pepembl(pepinp)
    ppmat  = dp(hlaemb, pepemb)

    # [2] : Contact quantization
    contacts0 = contactsl([hlarad, peprad])
    contacts  = subset_phla(contacts0)

    # [3] : Positional importance
    posimp0   = posimpl([hlainp, pepinp])
    posimp    = subset_phla(posimp0)

    # Intermediates that feed the output layer
    intrmat   = ppmat * contacts * posimp
    bstrength = msum(intrmat)

    if precomputed_pepeffects:
        pepeffect = Input((1,), name = 'pepeffect')
    else:
        pepeffect = pepeffectl(pepinp)

    pepeffect = pepeffectl(pepinp)
    bind_vec  = concat([bstrength, pepeffect])

    # Outputs
    regout    = regl(bstrength)
    bindout   = bindl(bind_vec)

    if conf == 'default':
        outputs = [regout, bindout]
    elif conf == 'extended':
        outputs = [regout, bindout, pepeffect, bstrength]
    elif conf == 'full':
        outputs = [regout, bindout, pepeffect, bstrength, intrmat]
    else:
        raise ValueError("Unknown conf")

    inps = [hlainpfull, pepinpfull]
    if precomputed_pepeffects:
        inps.append(pepeffect)
    confmod = keras.Model(inputs = inps, outputs = outputs)

    return confmod


def bind_mats_extract(bind_mats, odir):
    nitems = bind_mats['mats'].shape[0]
    for i in range(nitems):
        hla = bind_mats['hla'][i]
        pep = bind_mats['pep'][i]
        mat = np.round(bind_mats['mats'][i], 4)
        pkt = bind_mats['hla_pkts'][hla]
        rname = [f"hla{ind+1:02d}_{aa}" for ind, aa in enumerate(list(pkt))]
        cname = [f"pep{ind+1:02d}_{aa}" for ind, aa in enumerate(list(pep))]
        datf = pd.DataFrame(mat, columns = cname, index = rname)
        key = bind_mats['keys'][i]
        fbase = re.sub('[-*: ]', '-', key)
        fname = fbase + ".tsv"
        fpath = os.path.join(odir, fname)
        datf.to_csv(fpath, sep = '\t')

#------------------------------------------------#
#              Mutation generation               #
#------------------------------------------------#
def single_mut_gen(peptide, allele = None, aminos = None):
    if aminos is None:
        aminos = standard_aminos
    npos = len(peptide)
    mut_peps = []
    mut_pams = []
    for i in range(npos):
        orig_aa = peptide[i]
        pre_pep  = peptide[:i]
        post_pep = peptide[(i+1):]
        for aa in aminos:
            if aa != orig_aa:
                new_pep = pre_pep + aa + post_pep
                new_pam = PAM[orig_aa][aa]
                mut_peps.append(new_pep)
                mut_pams.append(new_pam)
            res = pd.DataFrame({
                'peptide'      : mut_peps,
                'pam'          : mut_pams,
                'orig_peptide' : peptide
            })
            if allele is not None:
                res.insert(0, 'hla', allele)
    return res
