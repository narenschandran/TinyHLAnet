# This script sets up the data for the benchmarking
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

MHCFLURRY="${BIN_DIR}/mhcflurry-run.sh"
DEEPHLAFFY="${BIN_DIR}/deephlaffy-quickrun.sh"
TPHLA="${BIN_DIR}/transphla-quickrun.sh"
MMPRED="${BIN_DIR}/mixmhcpred-run.sh"
NPAN="${BIN_DIR}/netmhcpan-run.sh"

bench_basedir="${PROJROOT}/results"
bench_dir="${bench_basedir}/02-benchmark"
bench_datadir="${bench_dir}/data"
bench_transphla_datadir="${bench_dir}/transphla-data"

bench_resdir="$bench_dir/results"
mkdir -p "$bench_resdir"

bench_tmpdir="$bench_dir/tmp"
mkdir -p "$bench_tmpdir"


#------------------------------------------------#
#                   DeepHLAffy                   #
#------------------------------------------------#
# DeepHLAffy ~ 2 minutes 45 seconds
for d in "${bench_datadir}/"*; do
    echo "$d"
    dbase=$(basename "$d")
    for f in "${d}/"*.tsv; do
        basef=$(basename "$f" '.tsv')
        odir="${bench_resdir}/${dbase}"
        mkdir -p "$odir"
        bash "$DEEPHLAFFY" "$f" "${odir}"
        bash "$MHCFLURRY" "$f" "${odir}"
        bash "$NPAN" "$f" "${odir}" "${bench_tmpdir}/netmhcpan/${dbase}/${basef}"
        bash "$MMPRED" "$f" "${odir}" "${bench_tmpdir}/mixmhcpred/${dbase}/${basef}"

        # TransPHLA requires a lot of resources for 1 million datapoints,
        # so we split it into two 500k files.
        if [ "$basef" = 'peps-6' ]; then
            tphla_tmp="${bench_tmpdir}/transphla/${dbase}/${basef}"
            tphla_tmpa="${tphla_tmp}a"
            mkdir -p "$tphla_tmpa"
            tphla_f1="${tphla_tmpa}/peps-6a.tsv"
            head -1 "$f" > "$tphla_f1"
            awk 'NR > 1' "$f" | head -500000 > "$tphla_f1"
            bash "$TPHLA" "$tphla_f1" "${odir}" "${tphla_tmpa}"

            tphla_tmpb="${tphla_tmp}b"
            mkdir -p "$tphla_tmpb"
            tphla_f2="${tphla_tmpa}/peps-6b.tsv"
            head -1 "$f" > "$tphla_f2"
            awk 'NR > 1' "$f" | tail -500000 > "$tphla_f2"
            bash "$TPHLA" "$tphla_f2" "${odir}" "${tphla_tmpb}"
        else
            bash "$TPHLA" "$f" "${odir}" "${bench_tmpdir}/transphla/${dbase}"
        fi
    done
done
