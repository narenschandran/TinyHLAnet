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

COLLECT="${script_dir}/prov-collect.py"
PLOT="${script_dir}/epoch-tracking-plots.R"

"$PY" "$COLLECT"
Rscript "$PLOT"
