script_path=${BASH_SOURCE[0]}
script_dir=$(dirname "$script_path")
projroot="${script_dir}/.."

env_dir="${projroot}/.env/deephlaffy"
if [ -d "$env_dir" ]; then
    source "${env_dir}/bin/activate"
    PY="${env_dir}/bin/python3"
else
    PY="python"
fi

TINYHLANET="${projroot}/tinyhlanet.py"

f="$1"
if [ -z "$f" ]; then
    echo "Usage: ${script_path} <input-tsv> <output-dir:optional>"
    kill -9 $$
fi


if [ ! -f "$f" ]; then
    echo "The input file <${f}> does not exist... Please check and try again."
    kill -9 $$
fi

odir="$2"
if [ -z "$odir" ]; then
    odir=$(pwd)
fi
mkdir -p "$odir"

ftmp="${f%.*}" # Remove compression extension, if any
fbase=$(basename "$ftmp" '.tsv')
log_f="${odir}/${fbase}.tinyhlanet.time"

/usr/bin/time -o "$log_f" --format "%e %M" "$PY" "$TINYHLANET" -o "$odir" "$f"
