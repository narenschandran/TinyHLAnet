import os
script_dir = os.path.dirname(__file__)
projroot   = os.path.join(script_dir, '..')

import sys
sys.path.append(projroot)

import argparse
import time
import numpy as np
import pandas as pd

from nbslpy._utils import savePickle
from deephlaffylib.Utils import layernames, layer_get
from deephlaffylib.RunUtils import (
    hla_inputs, pep_inputs, ModelLoad
)

parser = argparse.ArgumentParser(
    prog="DeepHLAffy Peptide Effect", description="Prediction of binding-independent peptide effect."
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
    help="Base string used for file name"
)

parser.add_argument(
    "-r",
    "--round",
    action="store_true",
    default=False,
    help="Round to 3 decimal points"
)

args = parser.parse_args()

f = args.filename
if args.output_dir is None:
    odir = os.path.dirname(f)
    if odir == '':
        odir = '.'
else:
    odir = args.output_dir

if args.file_base is None:
    # We do two splits in order to remove any compression
    # extensions, if any.
    fbase = os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] + '.deephlaffy.pepeffect'
else:
    fbase = args.file_base

tsv_fname = fbase + ".tsv.xz"
tsv_fpath = os.path.join(odir, tsv_fname)

pkl_fname = fbase + ".pkl"
pkl_fpath = os.path.join(odir, pkl_fname)

pepinp_datf = pd.read_csv(f, sep = '\t')
if 'peptide' not in pepinp_datf.columns:
    # If there is no "peptide" column, assume that the input
    # is a list of peptide, with one peptide present in each
    # line. We read this by subsetting the first "column" using
    # pandas. We use pandas because it can handle pretty
    # much any file compression automagically.
    pepinp_data = pd.read_csv(f, header = None).iloc[:,0]
else:
    pepinp_data = pepinp_datf.loc[:,"peptide"]


peptides = pepinp_data.unique()

start      = time.time()
mod        = ModelLoad()
pepeffectl = layer_get(mod, 'pepeffect')
castl      = layer_get(mod, 'i32_2')
pepinp     = castl(pep_inputs(peptides)[:,:,0])
pred       = pepeffectl.predict(pepinp)
end        = time.time()

if not os.path.exists(odir):
    os.makedirs(odir, exist_ok = True)

if args.round:
    val = np.round(pred[:,0], 3)
else:
    val = pred[:,0]

out_datf = pd.DataFrame({
    'peptide'   : peptides,
    'pepeffect' : val
})
out_datf.to_csv(tsv_fpath, sep = '\t', index = False)

out_dct = pd.Series(out_datf.pepeffect.values, index = out_datf.peptide).to_dict()
savePickle(out_dct, pkl_fpath)
