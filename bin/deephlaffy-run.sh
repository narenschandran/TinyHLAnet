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

DEEPHLAFFY="${projroot}/deephlaffy.py"

runtype="$1"
if [ -z "$runtype" ]; then
    runtype='quick'
fi

f="$2"
if [ -z "$f" ]; then
    echo "No input file detected..."
    echo "Usage: ${script_path} <runtype> <input-tsv> <output-dir:optional>"
    echo "Runtype can be one of: quick, extended or complete."
    kill -9 $$
fi

if [ ! -f "$f" ]; then
    echo "The input file <${f}> does not exist... Please check and try again."
    kill -9 $$
fi

odir="$3"
if [ -z "$odir" ]; then
    odir=$(pwd)
fi
mkdir -p "$odir"


ftmp="${f%.*}" # Remove compression extension, if any
fbase=$(basename "$ftmp" '.tsv')
log_f="${odir}/${fbase}.deephlaffy.time"
if [ "$runtype" = 'quick' ]; then
    /usr/bin/time -o "$log_f" --format "%e %M" "$PY" "$DEEPHLAFFY" -o "$odir" "$f"
elif [ "$runtype" = 'extended' ]; then
    /usr/bin/time -o "$log_f" --format "%e %M" "$PY" "$DEEPHLAFFY" -x -o "$odir" "$f"
elif [ "$runtype" = 'complete' ]; then
    /usr/bin/time -o "$log_f" --format "%e %M" "$PY" "$DEEPHLAFFY" -X -o "$odir" "$f"
elif [ "$runtype" = 'complete-extracted' ]; then
    /usr/bin/time -o "$log_f" --format "%e %M" "$PY" "$DEEPHLAFFY" -X -E -o "$odir" "$f"
else
    echo "Unknown runtype [$runtype]"
fi

