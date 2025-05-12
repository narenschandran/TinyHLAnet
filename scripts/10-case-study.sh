# This script sets up the data for the SARS-CoV-2 case study
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

# The next chunk decides which Python executable to call
env_dir="${PROJROOT}/.env/deephlaffy"
if [ -d "$env_dir" ]; then
    source "${env_dir}/bin/activate"
    PY="${env_dir}/bin/python3"
else
    PY="python"
fi

TINYHLANET="${PROJROOT}/tinyhlanet.py"
PROTEOME_SCAN="${PROJROOT}/tinyhlanet-scan"
ESCAPE_MUTANTS="${PROJROOT}/tinyhlanet-escape.py"
case_dir="${PROJROOT}/results/03-case-study"

sars_cov2_dir="${case_dir}"
mkdir -p "$sars_cov2_dir"

sars_cov2_data_dir="${sars_cov2_dir}/data"
mkdir -p "$sars_cov2_data_dir"

wget -P "$sars_cov2_dir" -N 'https://ars.els-cdn.com/content/image/1-s2.0-S2211124723000062-mmc2.xlsx'

SARS_CLEAN="${SCRIPT_DIR}/sars-cov-2-clean.R"
Rscript "$SARS_CLEAN"

wget -P "$sars_cov2_dir" -N 'http://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/pep/Sars_cov_2.ASM985889v3.pep.all.fa.gz'
seq_f="${sars_cov2_dir}/Sars_cov_2.ASM985889v3.pep.all.fa.gz"
gawk '
 /^>/ { print $1 }
!/^>/ { print $0 }
' <(zcat "$seq_f") > "${sars_cov2_data_dir}/sars-cov-2.proteome.fa"

gawk '
BEGIN { print "protein\tgene\ttranscript\tsymbol" }

/>/ {
    gene="UNSET"
    transcript="UNSET"
    sym="UNSET"

    prot=$1
    sub(/^>/, "", prot)
    for (i = 2; i <= NF; i++) {
        curr_word = $i
        if (curr_word ~ /gene:/) {
            sub(/gene:/, "", curr_word)
            gene = curr_word
        } else if (curr_word ~ /gene_symbol:/) {
            sub(/gene_symbol:/, "", curr_word)
            sym = curr_word
        } else if (curr_word ~ /transcript:/) {
            sub(/transcript:/, "", curr_word)
            transcript = curr_word
        }
    }
    print prot"\t"gene"\t"transcript"\t"sym
}' <(zcat "$seq_f") > "${sars_cov2_data_dir}/sars-cov-2.proteome.idmap"

bash "$PROTEOME_SCAN" \
    -l "${sars_cov2_data_dir}/allele-list.txt" \
    -a 0 -b 0 \
    -o "${sars_cov2_dir}/epitope-scan" \
    "${sars_cov2_data_dir}/sars-cov-2.proteome.fa"

"$PY" "$ESCAPE_MUTANTS" \
    -o "${sars_cov2_dir}/immune-escape" \
    -f 'escape-mutants' \
    "${sars_cov2_data_dir}/escape-study.tsv"


intstudy_f="${sars_cov2_data_dir}/interaction-study.tsv"

cat "${sars_cov2_data_dir}/escape-study.tsv" <(zcat "${sars_cov2_dir}/immune-escape/escape-mutants.tsv.gz" | awk 'NR > 1 {print $1"\t"$2}') > "$intstudy_f"

"$PY" "$TINYHLANET" -X -E -o "${sars_cov2_dir}/mutant-interactions" -f "mutant-interactions" "$intstudy_f"
