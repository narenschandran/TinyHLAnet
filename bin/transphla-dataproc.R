library(scriptloc)
library(seqinr)

script_file   <- scriptloc()
script_dir    <- script_dir_get()
projroot      <- file.path(script_dir, "..")
transphla_dir <- file.path(projroot, "bin",
                           "TransPHLA-AOMP")
seq_file      <- file.path(transphla_dir, "TransPHLA-AOMP",
                           "common_hla_sequence.csv")
if (!file.exists(seq_file)) stop("Unable to find file with HLA pseudo-sequences")
seq_df        <- read.table(seq_file, sep = ',', header = T,
                            stringsAsFactors = F)

args        <- commandArgs(trailingOnly = T)

if (length(args) != 2) {
    msg <- sprintf("<%s> <input-tsv> <output-dir>", script_file)
    writeLines(msg)
    stop("Exactly two inputs required")
}

# This input file must have the columns:
# <allele> and <peptide>
# <allele> must be in the normal HLA format: (eg) HLA-A*02:01
f   <- args[1]
if (!file.exists(f)) {
    msg <- sprintf("Unable to find input file: [%s]", f)
    writeLines(msg)
    stop("Check input 1")
}

# Read in the data and skip any cases where the amino acid is unknown0
dat <- subset(read.table(f, sep = '\t', header = T,
                         stringsAsFactors = F),
              !grepl("X", peptide),
              select = c("allele", "peptide"))

# Get the pseudosequence
seq_i <- match(dat$allele, seq_df$HLA)
dat$hla_pseudo <- seq_df$HLA_sequence[seq_i]

na_l <- is.na(seq_i)
n_na <- sum(na_l)
fin_dat <- if (n_na > 0) {
    na_hla   <- sort(unique(dat$allele[na_l]))
    n_na_hla <- length(na_hla)
    writeLines(sprintf("%d HLA IDs corresponding to %d entries have no available psuedosequence. Skipping them...", n_na_hla, n_na))
    writeLines("")
    writeLines(paste(na_hla, collapse = '\n'))
    dat[!na_l,]
} else {
    dat
}

# This is a key that we can use later
fin_dat$nm <- sprintf("phla%06d", seq_len(nrow(fin_dat)))


out_dir <- args[2]
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = T)

bname <- basename(tools::file_path_sans_ext(f, compression = T))
outhla   <- file.path(out_dir, paste0(bname, '.transphla.hla.fasta'))
outpep   <- file.path(out_dir, paste0(bname, '.transphla.pep.fasta'))

# These will be the inputs to TransPHLA
hla <- as.list(fin_dat$hla_pseudo)
pep <- as.list(fin_dat$peptide)

write.fasta(hla, fin_dat$allele, file.out = outhla)
write.fasta(pep, fin_dat$peptide, file.out = outpep)
