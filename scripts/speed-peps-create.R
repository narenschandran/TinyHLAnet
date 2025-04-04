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

NTHREADS <- 16
powers <- 1:6
ms <- 10 ^ (powers)
netmhc_par_base_odir <- file.path(bench_dir, 'nethmhcpan-parallel', 'data')
for (p in powers) {
    m <- ms[p]
    pepset <- peps[seq_len(m)]

    fname <- sprintf("peps-%d.tsv", p)
    fpath <- file.path(odir, fname)
    datf <- data.frame(
        allele = allele,
        peptide = pepset
    )
    write.table(datf, fpath, sep = '\t', row.names = F, quote = F)


    # This is for the netMHCpan parallel run
    sp_size <- ceiling(nrow(datf) / NTHREADS)
    sp_ind <- rep(seq_len(NTHREADS), each = sp_size)
    sp <- Filter(function(x) nrow(x) > 0, suppressWarnings(split(datf, sp_ind)))
    netmhc_par_odir <- file.path(netmhc_par_base_odir, sprintf('peps-%d', p))
    if (!dir.exists(netmhc_par_odir)) dir.create(netmhc_par_odir, recursive = T)

    for (i in seq_along(sp)) {
        ofile <- file.path(netmhc_par_odir,
                           sprintf("peps-%d%s.txt", i, letters[i]))
        write.table(sp[[i]], ofile, sep = '\t', row.names = F,
                    quote = F)
    }
}
