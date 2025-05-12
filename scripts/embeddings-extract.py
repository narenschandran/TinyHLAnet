import re
import os
script_dir = os.path.dirname(__file__)
projroot   = os.path.join(script_dir, '..')

import sys
sys.path.append(projroot)

import argparse
import time
import numpy as np
import pandas as pd

import tensorflow as tf

from nbslpy.aacodes import _aa1lst as aminos

from nbslpy._utils import savePickle
from deephlaffylib.Utils import layernames, layer_get
from deephlaffylib.RunUtils import (
    load_and_check_input,
    hla_inputs, pep_inputs,
    configure_model, hla_pkts, bind_mats_extract
)

def emb_get(emblayer, prefix = 'd'):
    mat   = np.round(emblayer.embeddings.numpy(), 3)
    datf  = pd.DataFrame(mat, index = aminos)
    ncols = datf.shape[1]
    datf.columns = [f"{prefix}{i+1:03d}" for i in range(ncols)]
    return datf

mod = configure_model("default")

emblayers = {
    'hla'             : layer_get(mod, 'hla'),
    'peptide'         : layer_get(mod, 'pep'),
    'pepeffects'      : layer_get(layer_get(mod, 'pepeffect'),
                                  'pep_effects'),
    'hla-context'     : layer_get(mod, 'posimp').hla_embmaker,
    'peptide-context' : layer_get(mod, 'posimp').pep_embmaker,
}

embs = {feat: emb_get(emblayer) for feat, emblayer in emblayers.items()}

base_odir = os.path.join(projroot, 'results', '04-dissection')
emb_odir = os.path.join(base_odir, 'embeddings')

for feat, emb in embs.items():
    emb_o = os.path.join(emb_odir, feat)
    emb_f = os.path.join(emb_o, f"{feat}-embeddings.tsv")
    if not os.path.exists(emb_o):
        os.makedirs(emb_o)
    emb.to_csv(emb_f, sep = '\t')


hemb = embs['hla'].values
pemb = embs['peptide'].values
dotp_np = np.round(hemb @ np.transpose(pemb, (1, 0)), 3)
dotp = pd.DataFrame(dotp_np, index = aminos, columns = aminos)

pp_odir = os.path.join(base_odir, 'pair-potential')
if not os.path.exists(pp_odir):
    os.makedirs(pp_odir)

pp_f = os.path.join(pp_odir, 'pair-potential.tsv')
dotp.to_csv(pp_f, sep = '\t')
