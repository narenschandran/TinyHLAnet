# This script is used to preprocess datasets to facilitate the model training.
# The chief purpose of preprocessing is to prevent resource wastage that would
# occur in converting input (human-sensible) datasets into a form accepted by
# the model. This would otherwise occur once per instance of model training,
# which adds up to an enormous amount of time.
#
# The input files are expected to be tab-separated and require the following
# four columns to be present:
# mhcpocket  - HLA-I pocket ordered from N-to-C end
# peptide    - Peptide sequence
# regressand - Numerical value representing binding affinity (scaled to
#              be between 0 and 1, both inclusive)
# binder     - Class label representing non-binder (0) or binder (1)
#
# To reiterate, the input file must have all these columns. In case regressand
# or binder information are not available for some row, the value (-1) is used
# the same.

import os
import argparse
import pandas as pd
import numpy as np

script_path = os.path.realpath(__file__)
script_dir = os.path.dirname(script_path)
projroot = os.path.join(script_dir, "..")
os.sys.path.append(projroot)

from nbslpy.seq import (
    seq2ind,
    seq2ohe,
    seq2blosum50,
    seq2oheblosum50,
    seq2aaradii_pmhc_ca,
)
from nbslpy._utils import savePickle


def proc(datf, seq_fn):
    """
    Helper function to process data. Converts peptide and amino
    acid sequences into forms that are usable by downstream models.

    Args:
        datf: DataFrame containing MHC pocket pseudosequence, peptide sequence, binder and regressand outputs.
        seq_fn: Function used to process input sequences.

    Output:
        Dictionary with data in suitable format.
    """

    # Since each allele is expected to have multiple data points, the specific
    # processing done for the MHC pocket sequence in these entries need not
    # be repeated, and can be memoized.
    mhc_dict = {pkt: seq_fn(pkt) for pkt in list(datf.mhcpocket.unique())}
    out = {
        "pep": np.asarray([seq_fn(pep) for pep in datf.peptide.values]),
        "hla": np.asarray([mhc_dict[mhc] for mhc in datf.mhcpocket.values]),
        "binder": np.asarray([bind for bind in datf.binder.values]),
        "regressand": np.asarray([bind for bind in datf.regressand.values]),
    }
    out["regressand"] = out["regressand"].reshape(out["regressand"].shape[0], 1)
    out["binder"] = out["binder"].reshape(out["binder"].shape[0], 1)
    return out


def seq2indradii(s):
    indrepr = np.expand_dims(seq2ind(s), -1)
    radrepr = np.expand_dims(seq2aaradii_pmhc_ca(s), -1)
    reprvec = np.concatenate([indrepr, radrepr], -1)
    return reprvec


def ind_proc(datf):
    """
    Processing function that converts amino acids to positive integers. Can be used for embeddings.
    """
    return proc(datf, seq2indradii)


def ohe_proc(datf):
    """
    Processing function that converts amino acids to one hot encodings.
    """
    return proc(datf, seq2ohe)


def bls50_proc(datf):
    """
    Processing function that converts amino acids to BLOSUM50 encodings.
    """
    return proc(datf, seq2blosum50)


def ohebls50_proc(datf):
    """
    Processing function that converts amino acids to concatenated one hot and BLOSUM50 encodings.
    """
    return proc(datf, seq2oheblosum50)


parser = argparse.ArgumentParser(
    prog="Dataset assemble", description="Assemble dataset for deep learning"
)

parser.add_argument("filename")
parser.add_argument(
    "-t",
    "--type",
    type=str,
    action="store",
    default="index",
    help="Information on how the sequence should be processed",
)
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
    "--fbase",
    type=str,
    action="store",
    default=None,
    help="Base name for output file without extensions",
)
parser.add_argument(
    "-r",
    "--randomize",
    type=int,
    action="store",
    default=None,
    help="Seed to randomize data",
)

args = parser.parse_args()

input_file = args.filename
if not os.path.exists(input_file):
    raise ValueError(f"The input file [{input_file}] does not exist.")

output_dir = args.output_dir
if output_dir is None:
    input_dir = os.path.dirname(os.path.realpath(input_file))
    output_dir = input_dir


proc_type = args.type.lower()
allowed_types = ["bls50", "index", "ohe", "ohebls50"]

if proc_type not in allowed_types:
    allowed_types_str = "\n".join(allowed_types)
    raise ValueError(f"Accepted type options are:\n{allowed_types_str}")

out_fbase = args.fbase
if out_fbase is None:
    # We remove the last extension from the input file
    input_fname = os.path.basename(input_file)
    rm_exts = ["txt", "tsv", "gz", "xz", "zip", "tar"]
    keep_comps = [comp for comp in input_fname.split(".") if comp not in rm_exts]
    input_fbase = ".".join(keep_comps)
    out_fbase = input_fbase

output_fname = f"{out_fbase}.{proc_type}.pkl"

x = pd.read_csv(input_file, sep="\t")

# Ensure that the input has peptide and MHC pocket information.
req_cols = ["peptide", "mhcpocket", "regressand", "binder"]

missing_cols = []
for req_col in req_cols:
    if req_col not in x.columns:
        missing_cols.append(req_col)

if len(missing_cols) > 0:
    mis_str = "\n".join(missing_cols)
    msg = f"The following columns are missing:\n{mis_str}"
    raise ValueError(msg)

# Select the function to be used for processing the input.
proc_fn_dict = {
    "bls50": bls50_proc,
    "ohe": ohe_proc,
    "ohebls50": ohebls50_proc,
    "index": ind_proc,
}
proc_fn = proc_fn_dict[proc_type]

# Randomize data if necessary, and then create data suitable
# for downstream models.
rand = args.randomize
if rand is not None:
    y = x.sample(frac=1, random_state=rand)
else:
    y = x
z = proc_fn(y)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Save processed data
output_file = os.path.join(output_dir, output_fname)
savePickle(z, output_file)

# Save original data
output_file2 = os.path.join(output_dir, f"{output_fname}source.gz")
y.to_csv(output_file2, sep="\t", index=False)
