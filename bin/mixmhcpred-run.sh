# This script processes the standard data frames that stores our
# inputs into the MixMHCpred peptide input format.
# [1] Input File
# [2] Output directory (Optional)
# [3] Temp directory (Optional)
# The input file. Should be tab-separated. First column is expected to
# have the HLA allele in the full format "HLA-<gene>*<num1>:<num2>".
# Second column should have the peptide sequence. The first line will be
# considered to be a header and skipped, so be mindful of that.

# This is a wrapper script to run TransPHLA for only the binder scores
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

MDIR="${BIN_DIR}/MixMHCpred"
menv_dir="${PROJROOT}/.env/mixmhcpred"
if [ -d "$menv_dir" ]; then
    source "${menv_dir}/bin/activate"
    PY=$(realpath -s "${tenv_dir}/bin/python3")
else
    PY=$(which python3)
fi

MMPRED="${MDIR}/MixMHCpred"
if [ ! -x "$MMPRED" ]; then
    echo "Unable to find the MixMHCpred executable at:"
    echo "[$MMPRED]"
    kill -9 $$
fi


MMP_PARSE="${BIN_DIR}/mixmhcpred-parse.sh"
if [ ! -x "$MMPRED" ]; then
    echo "Unable to find the MixMHCpred results parsing script at:"
    echo "[$MMP_PARSE]"
    kill -9 $$
fi

f="$1"
fbase=$(basename "$f")

odir="$2"
if [ -z "$odir" ]; then
    odir="$(pwd)/mixmhcpred-proc/$(echo $fbase | cut -d '.' -f1)"
fi

tmpdir="$3"
if [ -z "$tmpdir" ]; then
    tmpdir="${odir}/mixmhcpred-tmp-"$(mktemp -u XXXXXXXXX)
fi

mkdir -p "$odir"
mkdir -p "$tmpdir"

awk -vodir="$tmpdir" '

# This will remove the "*", which is not used in the NetMHCpan allele format
function alproc(s) {
    sub("[*]", "", s)
    return s
}

# This will be used to remove bad characters for the file path. The "_" can
# be converted to ":" to get the NetMHCpan allele format
function clean_fpath(s) {
    sub("[*]", "", s)
    sub("[:]", "_", s)
    return s
}

# We skip the header, and then process every line
NR > 1 {
    allele_str = $1
    allele     = alproc(allele_str)
    peptide    = $2
    fpath      = odir"/"clean_fpath(allele)".mixmhcpred-input.txt"
    print peptide > fpath
}' "$f"

tmp_main_outf="${tmpdir}/${fbase}"
bash "$MMP_PARSE" 'header' > "$tmp_main_outf"
 
for tmpf in "${tmpdir}/"*.mixmhcpred-input.txt; do
    basef=$(basename "$tmpf" '.mixmhcpred-input.txt')
    alname=$(echo "$basef" | sed 's/_/:/g')
    outf="${tmpdir}/${basef}.mixmhcpred-output.txt"
    outlf="${tmpdir}/${basef}.mixmhcpred-output.time"
    /usr/bin/time -o "$outlf" --format="%e %M" "$MMPRED"  -i "$tmpf" -o "$outf" -a "$alname"
    gene=$(echo "$alname" | sed 's/HLA[-]//' | sed 's/[0-9]*:[0-9]*$//')
    allele="HLA-${gene}"$(echo "$alname" | sed 's/^HLA[-][A-Z]*/*/')
    echo "$allele"
    bash "$MMP_PARSE" "$outf" "$allele" >> "$tmp_main_outf"
done


time_outf="${odir}/"$(basename "${fbase}" '.tsv')".mixmhcpred.time"
cat "${tmpdir}/"*.mixmhcpred-output.time | awk '
BEGIN {
    time = 0
    maxmem = 0
}
{
    time = time + $1
    if ($2 > maxmem) {
        maxmem = $2
    }
}
END {
    print time" "maxmem
}' > "$time_outf"


main_outf="${odir}/"$(basename "${fbase}" '.tsv')'.mixmhcpred.tsv'
Rscript -e "
f <- '$f'
procf <- '$tmp_main_outf'
outf  <- '$main_outf'

readf <- function(f) {
    x <- read.table(f, sep = '\t', header = T)
    x[['pmhc_key']] <- with(x, paste(allele, peptide))
    x
}

x <- readf(f)
y <- readf(procf)

ind <- match(x[['pmhc_key']], y[['pmhc_key']])
z <- y[ind,]

stopifnot(all.equal(x[['pmhc_key']], z[['pmhc_key']]))

x[,c('bindscore')] <- z[,c('bindscore')]

write.table(x, outf, sep = '\t', row.names = F, quote = F)
"

