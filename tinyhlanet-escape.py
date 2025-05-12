import os
import sys
script_dir = os.path.dirname(__file__)
projroot   = os.path.join(script_dir)
sys.path.append(projroot)

import argparse
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from deephlaffylib.Utils import md5sum_compute

parser = argparse.ArgumentParser(prog="Escape mutant scan")

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
    "-a",
    "--affinity-diff",
    type=float,
    action="store",
    default=0.2
)


parser.add_argument(
    "-b",
    "--bind-diff",
    type=float,
    action="store",
    default=0.3
)

parser.add_argument(
    "-p",
    "--pam50-cutoff",
    type=int,
    action="store",
    default=-1
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

args = parser.parse_args()

f     = args.filename
if args.output_dir is None:
    odir = 'immune-escape'
else:
    odir = args.output_dir

# We do two splits in order to remove any compression extensions, if any
if args.file_base is None:
    fbase = os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0] + '.deephlaffy.escape'
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

from deephlaffylib.RunUtils import (
        configure_model, single_mut_gen,
        hla_inputs, pep_inputs, check_input
)
mod = configure_model('default')

def predfn(datf, mod):
    datf = datf.copy()
    datf.reset_index(drop = True, inplace = True)
    hlainp   = hla_inputs(datf)
    pepinp   = pep_inputs(datf)
    pred     = mod.predict([hlainp, pepinp])
    datf.loc[:, "bindaff"]  = np.round(pred[0], 4)
    datf.loc[:, "bindprob"] = np.round(pred[1], 4)
    return datf


pdiff_cutoff = args.bind_diff
adiff_cutoff = args.affinity_diff
pam_cutoff   = args.pam50_cutoff

datf = check_input(pd.read_csv(f, sep = '\t'))
datf = datf.loc[:,['hla', 'peptide']].drop_duplicates().sort_values(['hla', 'peptide'])
datf.loc[:, "key"] = [v['hla'] + ' '+ v['peptide'] for n, v in datf.iterrows()]

datf = predfn(datf, mod)

mut_datf_lst = []
for num, row in datf.iterrows():
    tmp_mut_datf = single_mut_gen(row['peptide'],
                                  allele = row['hla'])
    tmp_mut_datf.loc[:,"orig_aff"]  = np.round(row["bindaff"], 3)
    tmp_mut_datf.loc[:,"orig_prob"] = np.round(row["bindprob"], 3)
    tmp_mut_datf.loc[:,"key"]       = row["key"]
    mut_datf_lst.append(tmp_mut_datf.copy())



mut_datf = predfn(pd.concat(mut_datf_lst, axis = 0), mod)
mut_datf.loc[:, "aff_diff"] = np.round(np.clip(mut_datf.loc[:,"orig_aff"] - mut_datf.loc[:, "bindaff"], 0, 1), 3)

mut_datf.loc[:, "prob_diff"] = np.round(np.clip(mut_datf.loc[:,"orig_prob"] - mut_datf.loc[:, "bindprob"], 0, 1), 3)

mut_datf_l =  (mut_datf.pam >= pam_cutoff) & ((mut_datf.prob_diff >= pdiff_cutoff) | (mut_datf.aff_diff >= adiff_cutoff))

mut_datf0 = mut_datf.copy()

mut_datf = mut_datf[mut_datf_l]


def rnkfn(x, descending = False):
    if descending:
        return rankdata(-x, method = 'min')
    else:
        return rankdata(x, method = 'min')

a = mut_datf.sort_values(['key', 'pam', 'aff_diff', 'prob_diff'],
                         ascending = [True, False, True, True])

ds = [d for _, d in a.groupby('key')]

fin_datf_lst = []
for _, e in a.groupby('key'):
    rnk = rnkfn(rnkfn(e.loc[:,"aff_diff"], True) * rnkfn(e.loc[:,"prob_diff"], True))
    e.loc[:,("rank")] = rnk
    fin = e.sort_values(['rank'])
    fin_datf_lst.append(fin.copy())

fin_datf = pd.concat(fin_datf_lst, axis = 0)
fin_datf.reset_index(drop = True, inplace = True)

cnames = list(fin_datf.columns)
cnames[1] = 'escape_peptide'
fin_datf.columns = cnames

if not os.path.exists(odir):

    os.makedirs(odir, exist_ok = True)
fin_datf.to_csv(outf, sep = '\t', index = False)


lock_fcon = open(outf_lock, 'w')
lock_fcon.write(chksum)
lock_fcon.close()
