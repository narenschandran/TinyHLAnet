# https://www.hiv.lanl.gov/content/immunology/variants/ctl_variant.csv
library(scriptloc)
script_dir <- script_dir_get()
projroot   <- file.path(script_dir, '..')
bench_dir  <- file.path(projroot, 'results', '02-benchmark')

f  <- file.path(bench_dir, 'lanl-hiv.csv')

# The second line is a timestamp, so we skip it.
hd <- gsub('"', "", strsplit(readLines(f, n = 1), ",")[[1]])
x  <- read.table(f, sep = ',', header = F, skip = 2,
                 stringsAsFactors = F)
colnames(x) <- hd

codes <- Reduce(union, strsplit(x$`Mutation Type Code`, ', '))

# These are the specific codes from the Mutation Type Codes
# from the dataset that will be used to
# https://www.hiv.lanl.gov/mojo/immunology/variant/mutation_types
ancodes <- local({
    sp1 <- strsplit(x$`Mutation Type Code`, ', ')
    sp2 <- strsplit(x$`Mutation Type Description`, ', ')

    stopifnot(all.equal(sapply(sp1, length), sapply(sp2, length)))

    an_df0 <- do.call("rbind.data.frame", lapply(seq_along(sp1), function(ind) {
        data.frame(
            code = sp1[[ind]],
            desc = sp2[[ind]]
        )
    }))
    an_df <- unique(an_df0)
    #setNames(an_df[[2]], an_df[[1]])
})

# Should be used for binding affinity check
# DHB  = Diminished HLA binding,
# 
# # Should check binder score for these
# E    = Escape,
# EL   = Epitope loss,
# LE   = Literature escape,
# 
# # Should check peptide-specific effects for these
# P    = Processing,
# 
# # Some evidence that diminished HLA binding, but not
# # always a functional assay.
# CHB  = Calculated diminished HLA binding,
# A    = HLA association,


# | Some evidence that escape is present in at least
# | some alleles, but not always a direct functional   |
# | assay.                                             |
# | CE   | Calculated escape                           |
# | IE   | Inferred escape                             |
# |----------------------------------------------------|
# | # Others                                           |
# |----------------------------------------------------|
# | AHE  | Altered HLA expression                      |
# | AKB  | Altered KIR binding                         |
# | HBOK | HLA binding unchanged                       |
# | R    | Reversion                                   |
# | I    | Insertion                                   |
# | ?    | Unclear sequence                            |
# | DI   | Drug induced                                |
# |----------------------------------------------------|
# |                                                    |
# |----------------------------------------------------|
# | CM   | Compensatory mutation                       |
# | F    | Fitness enhanced                            |
# | RCOK | Replicative capacity is not abrogated       |
# | RCR  | Replicative capacity reduced                |
# |----------------------------------------------------|
# |                                                    |
# |----------------------------------------------------|
# | OV   | Observed variant                            |
# | TCR  | TCR-related mutation                        |
# | DR   | Diminished response                         |
# |----------------------------------------------------|
# |                                                    |
# |----------------------------------------------------|
# | # Others                                           |
# | SF   | Susceptible form                            |
# | NSF  | Non-susceptible form                        |
# | SNSF | Subtype-specific non-susceptible form       |
# | SSF  | Subtype-specific susceptible form           |
# |----------------------------------------------------|

y <- subset(x, (HLA != ""))

valid_hla <- function(x) grepl("^[ABC][*][0-9][0-9]*:[0-9][0-9]*$", x)
al_l <- sapply(strsplit(y$HLA, ", "), function(x) any(valid_hla(x)))

y$alleles <- sapply(strsplit(y$HLA, ", "), function(al) {
    l <- valid_hla(al)
    if (sum(l) == 0) return("")
    if (sum(l) == 1) return(al[l])
    return(paste(al[l], collapse = ","))
})

valid_9mer <- function(x) grepl("^[A-Za-z]{9}$", x)
z <- subset(y, (alleles != "")      &
               valid_9mer(Epitope)  &
               valid_9mer(`Variant Epitope`))

single_z <- subset(z, !grepl(",", alleles))
multi_z  <- local({
    tmp <- subset(z, grepl(",", alleles))
    datf0 <- lapply(seq_len(nrow(tmp)), function(ind) {
        zz <- tmp[ind,]
        als <- strsplit(zz$alleles, ",")[[1]]
        tmp2 <- do.call("rbind.data.frame", lapply(seq_along(als), function(x) {
            zz
        }))
        tmp2$alleles <- als
        tmp2
    })
    datf <- do.call("rbind.data.frame", datf0)
    rownames(datf) <- NULL
    datf
})


clean_df <- rbind.data.frame(single_z, multi_z)


code_split <- list(
    'hla-binding'    = c("DHB"),
    'pep-processing' = c("P"),
    'epitope-escape' = c("E", "EL", "LE")
)

datf0 <- data.frame(
    allele      = toupper(paste0("HLA-", clean_df$alleles)),
    peptide     = toupper(clean_df$Epitope),
    mut_peptide = toupper(clean_df$`Variant Epitope`),
    mutcode     = clean_df$`Mutation Type Code`
)

datf_split <- lapply(code_split, function(code_sp) {
    sp <- strsplit(datf0$mutcode, ", ")
    l  <- sapply(sp, function(x) any(x %in% code_sp))
    datf0[l,]
})


bench_dir <- file.path(projroot, 'results', '02-benchmark')
hiv_bench_dir <- file.path(bench_dir, 'data', 'hiv')

if (!dir.exists(hiv_bench_dir)) dir.create(hiv_bench_dir, recursive = T)

for (nm in names(datf_split)) {
    datf_sp <- datf_split[[nm]]
    a <- do.call("rbind.data.frame", apply(datf_sp, 1, function(x) {
        y <- c(x, "label" = 1)
        tmp <- y
        tmp[["peptide"]] <- tmp[["mut_peptide"]]
        tmp[["label"]] <- 0
        rbind(y, tmp)[,-c(3, 4)]
    }, simplify = F))
    rownames(a) <- NULL
    a$set <- sprintf('%s-%03d', nm, rep(seq_len(nrow(a) / 2), each = 2))
    b <- unique(a)
    if (nm == "hla-binding") {
        b[["regressand"]] <- b[["label"]]
    } else if (nm == "epitope-escape") {
        b[["binder"]] <- b[["label"]]
    } else if (nm == "pep-processing") {
        b[["pep_processing"]] <- b[["label"]]
        b[["binder"]] <- b[["label"]]
    } else {
        stop("Unknown data split type")
    }
    fname <- paste0("hiv-", nm, ".tsv")
    fpath <- file.path(hiv_bench_dir, fname)
    write.table(b, fpath, sep = '\t', row.names = F, quote = F)
}
