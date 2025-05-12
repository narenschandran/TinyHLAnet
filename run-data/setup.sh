script_path=${BASH_SOURCE[0]}
script_dir=$(dirname "$script_path")
projroot="${script_dir}/.."
bin_dir="${projroot}/bin"
prereq_dir="${projroot}/prereq"

env_dir="${projroot}/.env/deephlaffy"
if [ -d "$env_dir" ]; then
    source "${env_dir}/bin/activate"
    PY="${env_dir}/bin/python3"
else
    PY="python"
fi

fbase=$(basename "$script_path" '.sh')
PRECOMPUTE="${script_dir}/setup.py"
"$PY" "$PRECOMPUTE"
