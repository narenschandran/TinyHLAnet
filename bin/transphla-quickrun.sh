# This is a wrapper script to run TransPHLA for only the binder scores
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

TDIR="${BIN_DIR}/TransPHLA-AOMP/TransPHLA-AOMP"
tenv_dir="${PROJROOT}/.env/transphla"
if [ -d "$tenv_dir" ]; then
    source "${tenv_dir}/bin/activate"
    PY=$(realpath -s "${tenv_dir}/bin/python3")
else
    PY=$(which python3)
fi
TPY="${TDIR}/pHLAIformer.py"

f="$1"
odir="$2"

if [ -z "$odir" ]; then
    odir=$(pwd)/transphla-run
fi
mkdir -p "$odir"

tmpdir="$3"
if [ -z "$tmpdir" ]; then
    tmpdir="${odir}/tmp"
fi

TPROC="${SCRIPT_DIR}/transphla-dataproc.R"
if [ ! -f "$TPROC" ]; then
    echo "Unable to find the TransPHLA input generation script at:"
    echo "[$TPROC]"
    kill -9 $$
fi

Rscript "$TPROC" "$f" "$tmpdir"
basef1=${f%.*} # Remove any compression extension
basef=$(basename "$basef1" '.tsv')
hfile="${tmpdir}/${basef}.transphla.hla.fasta"
pfile="${tmpdir}/${basef}.transphla.pep.fasta"

lfile="$odir/${basef}.transphla.time"
/usr/bin/time -o "$lfile" --format="%e %M" "$PY" "$TPY" --peptide_file "$pfile" --HLA_file "$hfile" --output_dir "$tmpdir"

cat "${tmpdir}/predict_results.csv" | tr ',' '\t' > "$odir/${basef}.transphla.tsv"
