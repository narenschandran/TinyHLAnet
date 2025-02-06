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

bench_netmhc_datadir="${bench_dir}/netmhcpan-data"
bench_transphla_datadir="${bench_dir}/transphla-data"

#--------------------------------------------------------------------#
#                         Benchmarking setup                         #
#--------------------------------------------------------------------#

#------------------------------------------------#
#               Train, Val & Test                #
#------------------------------------------------#
datasets_dir="${PROJROOT}/datasets/raw"
# for f in "${datasets_dir}/"*.tsv.xz; do
#     fbase=$(basename "$f" '.tsv.xz')
#     odir="${bench_datadir}/${fbase}"
#     mkdir -p "$odir"
#     outfbase="${odir}/${fbase}"
#     echo "$fbase"
#     outf="${outfbase}.tsv"
#     xzcat "$f" | awk '{print $1"\t"$2"\t"$4"\t"$5}' > "${outf}"
# done
# 
# 
# #------------------------------------------------#
# #                     Speed                      #
# #------------------------------------------------#
# # This creates the peptides for the speed benchmark
# SPEED_PEPS_CREATE="${SCRIPT_DIR}/speed-peps-create.R"
# Rscript "$SPEED_PEPS_CREATE"


#------------------------------------------------#
#                   SARS-CoV-2                   #
#------------------------------------------------#
sars_cov2_bench_url='https://ars.els-cdn.com/content/image/1-s2.0-S2666379121001555-mmc2.xlsx'
wget -P "$bench_dir" -N "$sars_cov2_bench_url"
sars_cov2_inpf="${bench_dir}/"$(basename "$sars_cov2_bench_url")

sars_cov2_bench_datadir="${bench_datadir}/sars-cov-2"
mkdir -p "$sars_cov2_bench_datadir"
sars_cov2_f="${sars_cov2_bench_datadir}/sars-cov2.tsv"
excl_f="${bench_dir}/sars-cov2-exclude.txt"

Rscript -e "
library(readxl)

f    <- '${sars_cov2_inpf}'
outf <- '${sars_cov2_f}'
excl <- readLines('${excl_f}')

x <- read_excel(f)

hla_pat <- 'HLA[-][ABC][*][0-9][0-9]*:[0-9][0-9]*$'
y <- local({
    # Data starts from second row
    tmp0 <- setNames(x[-1,c(2, 1)], c('allele', 'peptide'))
    tmp1 <- as.data.frame(tmp0[complete.cases(tmp0),])
    tmp2 <- subset(tmp1, nchar(peptide) == 9)
    tmp3 <- subset(tmp2, grepl(hla_pat, allele))
    tmp3[['pmhc_key']] <- with(tmp3, paste(allele, peptide))
    tmp4 <- tmp3[!(tmp3[['pmhc_key']] %in% excl),]
    tmp4 <- tmp4[with(tmp4, order(pmhc_key)),]
})

write.table(y, outf, sep = '\t', row.names = F, quote = F)
"
