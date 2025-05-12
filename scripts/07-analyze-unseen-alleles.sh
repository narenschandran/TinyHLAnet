# Script setup to activate environment and to select the correct
# Python version.
SCRIPT_PATH="${BASH_SOURCE[0]}"
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
PROJROOT="${SCRIPT_DIR}/.."
BIN_DIR="${PROJROOT}/bin"

TINYHLANET="${BIN_DIR}/tinyhlanet-quickrun.sh"

d="${PROJROOT}/results/01-model-tuning/04-unseen-alleles"
f="${d}/excluded.tsv.xz"

bash "$TINYHLANET" "$f" "$d"

outf="${d}/excluded.tinyhlanet.tsv.gz"
echo "$outf"

out_tiff="${d}/excluded-performance.tiff"

Rscript -e "
x <- read.table('$outf', sep = '\t', header = T)

pred <- ifelse(x[['pred_binder']] > 0.5,
               'True Positive', 'False Negative')
tbl <- local({
    tmp0 <- table(pred)
    tmp <- round((tmp0 / sum(tmp0)) * 100, 2)
    nms <- paste0(names(tmp), '\n(', unname(tmp), '%)')
    setNames(tmp, nms)
})

tiff('$out_tiff', width = 2500, height = 2500, res = 300)
par(family = 'symbol')
pie(tbl, col = c('indianred2', 'cornflowerblue'), radius = 0.3)
dev.off()
"

mogrify -compress lzw -trim "$out_tiff"
