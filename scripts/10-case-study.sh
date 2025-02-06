# This script sets up the data for the SARS-CoV-2 case study
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

DEEPHLAFFY="${BIN_DIR}/deephlaffy-run.sh"

case_dir="${PROJROOT}/results/03-case-study"
prereq_dir="${PROJROOT}/prereq"
f="${prereq_dir}/case-study.tsv"
odir="${case_dir}"
bash "$DEEPHLAFFY" 'complete-extracted' "$f" "${case_dir}"
