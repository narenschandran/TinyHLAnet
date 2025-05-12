# This is the set of scripts required to reproduce
# the entire model tuning & analysis for the
# manuscript.

Rscript scripts/00-data-split.R
bash scripts/01-inputs-generate.sh
bash scripts/02-regmodel-optimize-embdim-contacts.sh
bash scripts/03-regmodel-optimize-posconf.sh
bash scripts/04-effects.sh
Rscript scripts/05-model-tuning-stats.R
Rscript scripts/06-analyze-model-tuning.R 
bash scripts/07-analyze-unseen-alleles.sh
bash run-data/setup.sh
bash scripts/08-bench-setup.sh
bash scripts/09-benchmark-run.sh
bash scripts/10-case-study.sh
Rscript scripts/11-case-study-analyze.R
bash scripts/12-dissect-embeddings-and-effects.sh
Rscript scripts/13-plot-embeddings-and-effects.R
