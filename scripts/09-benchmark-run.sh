# This script runs the benchmarking for all the methods
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

MHCFLURRY="${BIN_DIR}/mhcflurry-run.sh"
DEEPHLAFFY="${BIN_DIR}/deephlaffy-quickrun.sh"
TPHLA="${BIN_DIR}/transphla-quickrun.sh"
NPAN="${BIN_DIR}/netmhcpan-run.sh"

bench_basedir="${PROJROOT}/results"
bench_dir="${bench_basedir}/02-benchmark"
bench_datadir="${bench_dir}/data"
bench_transphla_datadir="${bench_dir}/transphla-data"

bench_resdir="$bench_dir/results"
mkdir -p "$bench_resdir"

bench_tmpdir="$bench_dir/tmp"
mkdir -p "$bench_tmpdir"

# Benchmarks for DeepHLAffy, netMHCpan (single-threaded), MHCflurry & TransPHLA
for d in "${bench_datadir}/"*; do
    dbase=$(basename "$d")
    echo "$dbase"
    for f in "${d}/"*.tsv; do
        basef=$(basename "$f" '.tsv')
        odir="${bench_resdir}/${dbase}"
        mkdir -p "$odir"
        bash "$DEEPHLAFFY" "$f" "${odir}"
        bash "$MHCFLURRY" "$f" "${odir}"
        bash "$NPAN" "$f" "${odir}" "${bench_tmpdir}/netmhcpan/${dbase}/${basef}"

        # TransPHLA requires a lot of resources for 1 million datapoints,
        # so we split it into two 500k files.
        if [ "$basef" = 'peps-6' ]; then
            echo "$basef"
            tphla_tmp="${bench_tmpdir}/transphla/${dbase}/${basef}"
            tphla_tmpa="${tphla_tmp}a"
            mkdir -p "$tphla_tmpa"
            tphla_f1="${tphla_tmpa}/peps-6a.tsv"
            head -1 "$f" > "$tphla_f1"
            awk 'NR > 1' "$f" | head -500000 >> "$tphla_f1"
            bash "$TPHLA" "$tphla_f1" "${tphla_tmpa}" "${tphla_tmpa}"

            tphla_tmpb="${tphla_tmp}b"
            mkdir -p "$tphla_tmpb"
            tphla_f2="${tphla_tmpb}/peps-6b.tsv"
            head -1 "$f" > "$tphla_f2"
            awk 'NR > 1' "$f" | tail -500000 >> "$tphla_f2"
            bash "$TPHLA" "$tphla_f2" "${tphla_tmpb}" "${tphla_tmpb}"
            cat "${tphla_tmpa}/peps-6a.transphla.time" \
                "${tphla_tmpb}/peps-6b.transphla.time" |
            awk 'BEGIN {
                tot_time = 0
                ram_usage = 0
            } {
                tot_time = tot_time + $1
                if ($2 > ram_usage) {
                    ram_usage =$2
                }
            } END {
                print tot_time" "ram_usage
            }' > "${odir}/peps-6.transphla.time"

            head -1 "${tphla_tmpa}/peps-6a.transphla.tsv" > "${odir}/peps-6.tsv"
            awk 'NR > 1' "${tphla_tmpa}/peps-6a.transphla.tsv" >> "${odir}/peps-6.tsv"
            awk 'NR > 1' "${tphla_tmpb}/peps-6b.transphla.tsv" >> "${odir}/peps-6.tsv"
        else
            bash "$TPHLA" "$f" "${odir}" "${bench_tmpdir}/transphla/${dbase}/${basef}"
        fi
    done
done


par_datadir="${bench_dir}/netmhcpan-parallel/data"
par_odir="${bench_dir}/netmhcpan-parallel/output"
mkdir -p "$par_odir"
par_base_idir="${bench_dir}/netmhcpan-parallel/intermediate"
mkdir -p "$par_base_idir"


par_cmd_dir="${bench_dir}/netmhcpan-parallel/cmds"
mkdir -p "$par_cmd_dir"

run_par() {
    local f="$1"
    cat "$f" | parallel -j16
}

# Speed benchmarking for netMHCpan (multi-threaded)
for d in "${par_datadir}/"*; do
    dbase=$(basename "$d")
    odir="${bench_tmpdir}/netmhcpan-parallel/${dbase}"

    # Has the list of files that will be used for the
    # parallelization script.
    tmpfile="${par_base_idir}/${dbase}.txt"
    lfile="${par_odir}/${dbase}.time"

    # Is the parallelization script.
    tmpscr="${par_cmd_dir}/${dbase}.sh"

    echo "$d"
    for f in "${d}/"*.txt; do
        fbase=$(basename "$f" '.txt')
        tmpdir="${bench_tmpdir}/netmhcpan-parallel/${fbase}"
        echo bash "$NPAN" "$f" "${odir} ${tmpdir}"
    done > "$tmpfile"
    echo cat "$tmpfile" '| parallel -j16' > "$tmpscr"
    mkdir -p "$odir"
    mkdir -p "$tmpdir"
    /usr/bin/time -o "${lfile}" --format="%e %M" bash "$tmpscr"
done

PLOTTING="${SCRIPT_DIR}/benchmark-plots.R"
Rscript "$PLOTTING"
