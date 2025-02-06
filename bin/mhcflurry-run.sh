script_path=${BASH_SOURCE[0]}
script_dir=$(dirname "$script_path")
projroot="${script_dir}/.."
bin_dir="${projroot}/bin"
datasets_dir="${projroot}/datasets"

env_dir="${projroot}/.env/mhcflurry"
if [ -d "$env_dir" ]; then
    source "${env_dir}/bin/activate"
    PY="${env_dir}/bin/python3"
else
    PY="python"
fi

f="$1"
odir="$2"
if [ -z "$odir" ]; then
    odir=$(pwd)
fi

fbase=$(basename "$f" '.tsv')'.mhcflurry'
proc_f="${odir}/${fbase}.csv"
cat "$f" | tr '\t' ',' | grep -v "X" > "$proc_f"

log_f="${odir}/${fbase}.time"

pred_f="${odir}/${fbase}.csv"
/usr/bin/time -o "$log_f" --format="%e %M" mhcflurry-predict "$proc_f" --out "$pred_f"
