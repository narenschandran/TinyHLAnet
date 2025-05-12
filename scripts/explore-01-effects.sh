SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"


env_dir="${PROJROOT}/.env/deephlaffy"
if [ -d "$env_dir" ]; then
    source "${env_dir}/bin/activate"
    PY="${env_dir}/bin/python3"
else
    PY="python"
fi

PEPEFFECT="${BIN_DIR}/deephlaffy-pepeffect.py"
RUNDATA_DIR="${PROJROOT}/run-data"

res_dir="${PROJROOT}/results/04-dissection"
mkdir -p "$res_dir"
input_f="${res_dir}/dataset.tsv"
input_xz="${input_f}.xz"

odir="${RUNDATA_DIR}/pepeffects"
mkdir -p "$odir"

datasets_dir="${PROJROOT}/datasets/raw"
xzcat "${datasets_dir}/train.tsv.xz" | head -1 > "$input_f"
for f in "${datasets_dir}/"*.tsv.xz; do
    awk 'NR > 1' <(xzcat "$f")
done >> "$input_f"
xz -9 "$input_f"

DEEPHLAFFY="${BIN_DIR}/deephlaffy-run.sh"

bash "${DEEPHLAFFY}" 'extended' "$input_xz" "${res_dir}/pepeffect"

"$PY" "${SCRIPT_DIR}/embeddings-extract.py"
