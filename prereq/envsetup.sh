#! /usr/bin/env bash

#' This script generates the Python virtual environemnts used during
#' training and benchmarking DeepHLAffy in a Linux environemnt

SCRIPT_PATH=${BASH_SOURCE[0]}
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"
PREREQ_DIR="${PROJROOT}/prereq"
ENVBASEDIR="${PROJROOT}/.env"
mkdir -p "${ENVBASEDIR}"

# We fix a specific version of Python by hand-installing it.
bash "${BIN_DIR}/install-python3.11.11" "${BIN_DIR}"

PY="${BIN_DIR}/python-3.11.11/bin/python3"
ENVDIR="${ENVBASEDIR}/deephlaffy"
ENVACTV="${ENVDIR}/bin/activate"
ENVPY="${ENVDIR}/bin/python3"
"$PY" -m venv "${ENVDIR}"

# By using a subshell, we the enrivonment activation is contained
# and does not interfere when we set up the other environments
(
    source "${ENVACTV}"
    "$ENVPY" -m pip install -r "${PREREQ_DIR}/requirements.txt"
)

TPHLA_ENVDIR="${ENVBASEDIR}/transphla"
TPHLA_ENVACTV="${TPHLA_ENVDIR}/bin/activate"
TPHLA_ENVPY="${TPHLA_ENVDIR}/bin/python3"
"$PY" -m venv "$TPHLA_ENVDIR"

# The subshell here serves the same purpose as earlier and prevents
# spillover of the environment onto the next steps of the script.
# 
# We use the specific commit of TransPHLA that was most recent when
# we were benchmarking.
(
    source "$TPHLA_ENVACTV"
    cd "${TPHLA_ENVDIR}/bin"
    git clone https://github.com/a96123155/TransPHLA-AOMP.git
    cd TransPHLA-AOMP/TransPHLA-AOMP

    # Fix indentation error
    sed -i 's/        log = Logger(errLogPath)/    log = Logger(errLogPath)/' pHLAIformer.py

    # Fix bug in iteration variable (hla->hla_seq)
    sed -i 's/    if not (pep.isalpha() and hla.isalpha()):/    if not (pep.isalpha() and hla_seq.isalpha()):/' pHLAIformer.py
)


MHCFL_ENVDIR="${ENVBASEDIR}/mhcflurry"
MHCFL_ENVACTV="${MHCFL_ENVDIR}/bin/activate"
MHCFL_ENVPY="${MHCFL_ENVDIR}/bin/python3"
"$PY" -m venv "$MHCFL_ENVDIR"
(
    source "$MHCFL_ENVACTV"
    "$MHCFL_ENVPY" -m pip install mhcflurry==2.1.2
    "$MHCFL_ENVPY" -m pip install numpy==1.23.5
    mhcflurry-downloads fetch
    deactivate
)
