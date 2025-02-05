library(scriptloc)
library(seqinr)

script_dir <- script_dir_get()
projroot   <- file.path(script_dir, '..')
res_dir    <- file.path(projroot, 'results')
bench_dir  <- file.path(res_dir, '02-benchmark')

# Amino acids used to compose the peptides
aminos <- c("A", "C", "D", "E", "F", 
            "G", "H", "I", "K", "L", 
            "M", "N", "P", "Q", "R",
            "S", "T", "V", "W", "Y")


#' Inputs:
#' [1]: Seed
#' [2]: Return as string?
pepmaker <- function(seed, as_str) {
    set.seed(seed)
    pep <- sample(aminos, 9, replace = T)
    if (as_str) paste(pep, collapse = "") else pep
}


# We want a million peptides. We make a few more than that
# to account for duplicate creations
N <- 1500000
system.time({
z = sapply(seq_len(N), pepmaker, as_str = T)
})

# We remove duplicates
v = z[!duplicated(z)]

# This is our actual goal: a million peptides
M <- 1000000
stopifnot(M <= length(v))

# Our peptide set
peps <- v[seq_len(M)]

odir <- file.path(bench_dir, 'data', 'speed')
if (!dir.exists(odir)) dir.create(odir, recursive = T)

allele        <- "HLA-A*02:01"
transphla_pkt <- 'YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY'

powers <- 1:6
ms <- 10 ^ (powers)
for (p in powers) {
    m <- ms[p]

    # Used for netMHCpan & MixMHCpred
    netmhcpan_odir <- file.path(odir, '..', '..', 'netmhcpan-data', 'speed')
    if (!dir.exists(netmhcpan_odir)) dir.create(netmhcpan_odir, recursive = T)
    fname <- sprintf("peps-%d.txt", p)
    fpath <- file.path(netmhcpan_odir, fname)
    pepset <- peps[seq_len(m)]
    writeLines(pepset, fpath)

    # Used for DeepHLAffy & MHCflurry
    fname2 <- sprintf("peps-%d.tsv", p)
    fpath2 <- file.path(odir, fname2)
    datf <- data.frame(
        allele = allele,
        peptide = pepset
    )
    write.table(datf, fpath2, sep = '\t', row.names = F, quote = F)

    # Used for TransPHLA
    tphla_odir <- file.path(odir, '..', '..', 'transphla-data', 'speed')

    if (!dir.exists(tphla_odir)) dir.create(tphla_odir, recursive = T)
    p_fname <- sprintf("peps-%d.fa", p)
    p_fpath <- file.path(tphla_odir, p_fname)
    write.fasta(as.list(pepset), peps, p_fpath)


    h_fname <- sprintf("hla-%d.fa", p)
    h_fpath <- file.path(tphla_odir, h_fname)
    hvec  <- rep(transphla_pkt, length(pepset))
    hname <- rep(allele, length(pepset))
    write.fasta(as.list(hvec), hname, h_fpath)
}
