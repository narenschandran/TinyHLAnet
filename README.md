TinyHLAnet: A light-weight 3D structure-aware architecture for rapid and explainable identification of CD8+ T-cell antigens
====

## Before running any of the scripts in this repository, please check the sub-section titled: `Setting up a Python/Conda environment`. It details the easiest way to install the required packages (with matched versions). Even if you're using your custom environment to execute the code, that subsection will provide hints as to the specific packages that are absolutely essential.

## Quick start guide

The TinyHLAnet toolkit comes with three use interfaces. This section will showcase quick examples for all three.


### `tinyhlanet.py`

The first and simplest of these is the `tinyhlanet.py` script, which can be used to predict the binding affinity (S score) and immunopeptidome presentation probability (S\*) scores for pairs of peptide & HLA-I molecules. The input to this script is a tab-separated file which has two columns: `allele` and `peptide`. An example file can be found in `example/data.tsv`. The following is a typical command that you might be interested in:

```
python tinyhlanet.py -o example/quick-run example/data.tsv
```
The file `example/quick-run/data.tinyhlanet.tsv.gz` should have been generated at this point with four columns: `hla`, `peptide`, `pred_regressand` (S score), and `pred_binder` (S* score). As a general rule of thumb, any epitope that has S score >= 0.2 and S* score >= 0.5 can be considered as a probable epitope. It is not recommended to reduce these cutoffs beyond these values.

The `tinyhlanet.py` script can also provide detailed information about the specific reisude pairs that contribute most to pHLA-I complex formation. These can be accessed by running the following command:
```
python tinyhlanet.py -o example/quick-run-with-interactions -X -E example/data.tsv
```
This will generate binding matrices with the strength of interaction between every intermolecular residue pair in the pHLA-I complex in the `example/quick-run-with-interactions/binding-mats/` directory. For instance, the file 'HLA-A-02-01-HLMLQLVRV.tsv' will have the interaction between the `HLA-A*02:01` allele with the `HLMLQLVRV` peptide.

You can compare the results for both these commands with a set of pre-generated results placed in the `examples/pre-run` directory if you want to carry out a sanity check.

These are the two most common use cases for the `tinyhlanet.py` interface. The full set of options can be accessed using the command:

```
python tinyhlanet.py -h
```

### Repository organization
This repository contains the code and data required to run `TinyHLAnet` model for prediction of CD8+ T-cell epitopes. The files/directories that come with this repository include:

```
TinyHLAnet
|
|-bin             : Contains wrapper scripts used during model benchmarking
|-deephlaffylib   : Python library that contains the code used to built TinyHLAnet
|-nbslpy          : Python library that houses general-purpose and protein processing functions
|-prereq          : Contains preprocessed data that will be used to kickstart the manuscript-related analysis
|-scripts         : Contains the scripts required to reproduce the manuscript-related analysis & data from scratch.
|-tools           : Contains wrapper scripts used for TinyHLAnet-scan
|
|
|-reproduce.sh    : One shot script for reproducing the manuscript-related analysis.
|-tinyhlanet.py: Prediction of pHLA-I binding affinity (S) and immunopeptidome presentation probability (S\*). 
|-tinyhlanet-scan: Comprehensive proteome-wide scanning using multiple HLA-I alleles.
|-tinyhlanet-escape: Prediction of probable immune escape mutants.
```

The `datasets` and `results` directories will be created automatically if the `reproduce.sh` script is executed.
