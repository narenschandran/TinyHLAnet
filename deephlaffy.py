import re
import os
script_dir = os.path.dirname(__file__)
projroot   = os.path.join(script_dir)

import sys
sys.path.append(projroot)

import argparse
import time
import numpy as np
import pandas as pd
from nbslpy._utils import savePickle
from deephlaffylib.Utils import md5sum_compute

parser = argparse.ArgumentParser(
    prog="DeepHLAffy", description="Prediction of pHLA-I complex formation"
)

parser.add_argument("filename")
parser.add_argument(
    "-o",
    "--output-dir",
    type=str,
    action="store",
    default=None,
    help="Output directory",
)

parser.add_argument(
    "-f",
    "--file-base",
    type=str,
    action="store",
    default=None,
    help="Basename for file."
)

parser.add_argument(
    "-n",
    "--no-overwrite",
    action="store_true",
    default=False,
    help="Don't rerun if a prior run is detectable. This will not run any version of the program (baseline, -x, -X), if any version is detected"
)

parser.add_argument(
    "-X",
    "--full_output",
    action="store_true",
    default=False,
    help="Output all computations, including the binding matrix & binding-independet HLA and peptide effects."
)

parser.add_argument(
    "-E",
    "--extract",
    action="store_true",
    default=False,
    help="Extract binding affinity from .pkl file into matrices. Only works if -X/--full_output is specified."
)

parser.add_argument(
    "-x",
    "--extended_output",
    action="store_true",
    default=False,
    help="Output binding affinity, presentation probability, HLA and peptide effects."
)
args = parser.parse_args()

if args.output_dir is None:
    odir = os.path.dirname(f)
    if odir == '':
        odir = '.'
else:
    odir = args.output_dir

if args.full_output:
    conf = 'full'
elif args.extended_output:
    conf = 'extended'
else:
    conf = 'default'

f     = args.filename
# We do two splits in order to remove any compression extensions, if any
if args.file_base is None:
    fbase = os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] + '.deephlaffy'
else:
    fbase = args.file_base

fname = fbase + ".tsv.gz"
outf  = os.path.join(odir, fname)

fname_lock = fbase + ".input.md5"
outf_lock  = os.path.join(odir, fname_lock)
chksum     = md5sum_compute(f)

if args.no_overwrite:
    if os.path.exists(outf_lock):
        fcon = open(outf_lock, 'r')
        msum = fcon.read()
        print(msum)
        fcon.close()
        if msum == chksum:
            print(f"The MD5sum specified in: [{outf_lock}] matches with the current input file: [{f}]. Skipping run...", file = sys.stderr)
            sys.exit(0)

from deephlaffylib.Utils import layernames, layer_get, md5sum_compute
from deephlaffylib.RunUtils import (
    load_and_check_input,
    hla_inputs, pep_inputs,
    configure_model, hla_pkts, bind_mats_extract
)

datf  = load_and_check_input(f)


hlainp = hla_inputs(datf)
pepinp = pep_inputs(datf)

start = time.time()
mod   = configure_model(conf)
pred  = mod.predict([hlainp, pepinp])

end   = time.time()
elaps = end - start
print(f"Time elapsed: {elaps:0.2f}s")

if not os.path.exists(odir):
    os.makedirs(odir, exist_ok = True)

datf.loc[:,"pred_regressand"] = np.round(pred[0], 3)
datf.loc[:,"pred_binder"]     = np.round(pred[1], 3)

# We get the weigts for extended/full computation
rwts = layer_get(mod, 'reg').weights
rsc  = rwts[0].value.numpy()[0][0]

bwts = layer_get(mod, 'bind').weights
bsc1 = bwts[0].value.numpy()[0][0]
bsc2 = bwts[0].value.numpy()[1][0]


bbias = bwts[1].value.numpy()[0]


if (conf == 'full') or (conf == 'extended'):
    bstren_reg  = pred[3] * rsc
    bstren_bind = pred[3] * bsc1
    peffect = pred[2] * bsc2
    datf.loc[:,"OrientedBindingStrengthReg"]  = np.round(bstren_reg, 3)
    datf.loc[:,"OrientedBindingStrengthBind"] = np.round(bstren_bind, 3)
    datf.loc[:,"OrientedPepEffect"]  = np.round(peffect, 3)

datf.to_csv(outf, sep = '\t', index = False)

lock_fcon = open(outf_lock, 'w')
lock_fcon.write(chksum)
lock_fcon.close()

# We always ensure that the binding matrix we write
# out is sign corrected to have attractive interactions
# being assigned positive values
if conf == 'full':
    pred[4] = pred[4] * np.sign(rsc)
    bind_mats = {
        'hla' : datf.hla.values.tolist(),
        'pep' : datf.peptide.values.tolist(),
        'keys' : (datf.hla + " " + datf.peptide).values.tolist(),
        'mats' : pred[4]
    }
    pkts_tmp = hla_pkts[hla_pkts.allele.isin(bind_mats['hla'])]
    bind_mats['hla_pkts'] = {a[1]['allele'] : a[1]['mhcpocket'] for a in pkts_tmp.iterrows()}
    pkl_f = os.path.join(odir, fbase + '-binding-mats.pkl')
    savePickle(bind_mats, pkl_f)

    if args.extract:
        bmats_odir = os.path.join(odir, 'binding-mats')
        os.makedirs(bmats_odir, exist_ok = True)
        bind_mats_extract(bind_mats, bmats_odir)
