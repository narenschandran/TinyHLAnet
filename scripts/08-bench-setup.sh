# This script sets up the data for the benchmarking
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

NPROC="${BIN_DIR}/netmhcpan-dataproc.sh"
TPROC="${BIN_DIR}/transphla-dataproc.R"

bench_basedir="${PROJROOT}/results"

bench_dir="${bench_basedir}/02-benchmark"

bench_datadir="${bench_dir}/data"
mkdir -p "$bench_datadir"


#--------------------------------------------------------------------#
#                         Benchmarking setup                         #
#--------------------------------------------------------------------#

#------------------------------------------------#
#               Train, Val & Test                #
#------------------------------------------------#
datasets_dir="${PROJROOT}/datasets/raw"
for f in "${datasets_dir}/"*.tsv.xz; do
    fbase=$(basename "$f" '.tsv.xz')
    odir="${bench_datadir}/${fbase}"
    mkdir -p "$odir"
    outfbase="${odir}/${fbase}"
    echo "$fbase"
    outf="${outfbase}.tsv"
    xzcat "$f" | awk '{print $1"\t"$2"\t"$4"\t"$5}' > "${outf}"
done


#------------------------------------------------#
#                     Speed                      #
#------------------------------------------------#
# This creates the peptides for the speed benchmark
SPEED_PEPS_CREATE="${SCRIPT_DIR}/speed-peps-create.R"
Rscript "$SPEED_PEPS_CREATE"
