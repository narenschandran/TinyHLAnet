# This script does some precomputation so that future runs are easier.
# 1. Since we have the list of all HLA alleles and their pockets, we
#    preprocess them into the model input format so that future runs
#    can just use a dictionary call rather than do the processing again.
# 2. To aid processing peptide sequences during runs, we compute
#    the input for triple amino acids so that future runs can piece together
#    the input for full peptides using the triplet processing
import os
import time
import pandas as pd
import sys

script_dir = os.path.dirname(__file__)
projroot   = os.path.join(script_dir, '..')

sys.path.append(projroot)
from nbslpy.aacodes import _aa1lst as AA
from nbslpy._utils import savePickle
from deephlaffylib.RunUtils  import *

mod = ModelLoad()

# We first carry out the processing for HLA pocket sequences
pkt_dct = {}
al_dct  = {}
for i, dat in hla_pkts.iterrows():
    al  = dat['allele']
    pkt = dat['mhcpocket']
    if pkt not in pkt_dct:
        pkt_dct[pkt] = seq_process(pkt)
    al_dct[al] = pkt_dct[pkt]

alleles = sorted(list(al_dct.keys()))

odir = os.path.join(projroot, 'run-data')
if not os.path.exists(odir):
    os.makedirs(odir, exist_ok = True)
savePickle(al_dct, os.path.join(odir, 'hlainp.pkl'))

# Triplet processing
triplets = [aa1 + aa2 + aa3 for aa1 in AA for aa2 in AA for aa3 in AA]
triplets_dct = {triplet: seq_process(triplet) for triplet in triplets}
savePickle(triplets_dct, os.path.join(odir, 'triplets.pkl'))
