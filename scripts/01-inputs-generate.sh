# This script preprocesses the different dataset splits so that input processing
# is not repeated each time a specific model instance is trained later on
script_path=${BASH_SOURCE[0]}
script_dir=$(dirname "$script_path")
projroot="${script_dir}/.."
bin_dir="${projroot}/bin"
datasets_dir="${projroot}/datasets"

# The next chunk decides which Python executable to call
env_dir="${projroot}/.env/deephlaffy"
if [ -d "$env_dir" ]; then
    source "${env_dir}/bin/activate"
    PY="${env_dir}/bin/python3"
else
    PY="python"
fi
ASSEMBLER="${bin_dir}/data-assemble.py"

# Our analysis converts sequence data into indices that are used for embedding.
# This is specified by the "index" proc_type
proc_types="index"
odir="${datasets_dir}/proc"
mkdir -p "$odir"

while read f; do
    for proc_type in $(echo "$proc_types"); do
        echo "Processing [$f] data into [$proc_type] type..."
        "${PY}" "$ASSEMBLER" -t "$proc_type" -o "$odir" "$f"
    done
done << EOF
${datasets_dir}/raw/train.tsv.xz
${datasets_dir}/raw/val.tsv.xz
${datasets_dir}/raw/test.tsv.xz
EOF
