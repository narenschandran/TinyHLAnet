library(scriptloc)
library(readxl)

SCRIPT_DIR   <- script_dir_get()
PROJROOT     <- file.path(SCRIPT_DIR, "..")

case_dir     <- file.path(PROJROOT, 'results', '03-case-study')
f <- file.path(case_dir, '1-s2.0-S2211124723000062-mmc2.xlsx')


case_data_dir <- file.path(case_dir, 'data')
if (!dir.exists(case_data_dir)) dir.create(case_data_dir, recursive = T)

dat_f <- file.path(case_data_dir, 'sars-cov-2.epitopes')
al_f  <- file.path(case_data_dir, 'allele-list.txt')
escape_f <- file.path(case_data_dir, 'escape-study.tsv')


x    <- read_excel(f, sheet = 1)
dat0 <- local({
    hd             <- as.character(x[5,])
    tmp0           <- as.data.frame(x[6:nrow(x),])
    colnames(tmp0) <- hd
    sel_cols       <- c(
        "allele"   = "HLA allele",
        "epitope"  = "Peptide seq",
        "method"   = "Selection method",
        "protname" = "Protein"
    )
    tmp1           <- tmp0[, unname(sel_cols)]
    tmp2           <- setNames(tmp1, names(sel_cols))
    transform(tmp2, allele = paste0("HLA-", allele))
})

dat <- subset(dat0,
               grepl("MS", method) &
              (nchar(epitope) == 9))
dat <- dat[with(dat, order(allele, epitope, method)),]
writeLines(sort(unique(dat$allele)), al_f)
write.table(dat, dat_f, sep = '\t', row.names = F, quote = F)


mut_study <- data.frame(
    allele  = "HLA-A*02:01",
    peptide = c("YLQPRTFLL", "GLMWLSYFI",
                "LLFNKVTLA", "SIIAYTMSL")
)
write.table(mut_study, escape_f, sep = '\t',row.names = F, quote = F)
