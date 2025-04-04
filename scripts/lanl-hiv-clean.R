# https://www.hiv.lanl.gov/content/immunology/variants/ctl_variant.csv
library(scriptloc)
script_dir <- script_dir_get()
projroot   <- file.path(script_dir, '..')
case_dir  <- file.path(projroot, 'results', '03-case-study')
prereq_dir <- file.path(projroot, 'prereq')

excl_f <- file.path(case_dir, '..', '02-benchmark', 'exclude-keys.txt')
excl <- readLines(excl_f)

f  <- file.path(case_dir, 'lanl-hiv.csv')

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

valid_9mer <- function(x) grepl("^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]{9}$", x)
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
    'epitope-escape' = c("E", "EL", "LE")
)

datf0 <- data.frame(
    allele      = toupper(paste0("HLA-", clean_df$alleles)),
    peptide     = toupper(clean_df$Epitope),
    mut_peptide = toupper(clean_df$`Variant Epitope`),
    mutcode     = clean_df$`Mutation Type Code`
)


excl_l <- with(datf0,
    (paste(allele, peptide) %in% excl)|
    (paste(allele, mut_peptide) %in% excl))

datf1 <- datf0[!excl_l,]
datf_split <- lapply(code_split, function(code_sp) {
    sp <- strsplit(datf1$mutcode, ", ")
    l  <- sapply(sp, function(x) any(x %in% code_sp))
    datf1[l,]
})

hiv_case_dir <- file.path(case_dir, 'data', 'hiv')

if (!dir.exists(hiv_case_dir)) dir.create(hiv_case_dir, recursive = T)

nm <- 'epitope-escape'
datf_sp <- datf_split[[nm]]

dup_l   <- duplicated(apply(datf_sp[,1:3], 1, paste, collapse = " "))

datf_out0 <- datf_sp[!dup_l,]
a1 <- do.call("rbind", strsplit(datf_out0$peptide, ""))
a2 <- do.call("rbind", strsplit(datf_out0$mut_peptide, ""))
single_l <- rowSums(a1 != a2) == 1
datf_out <- datf_out0[single_l,]
datf_out <- datf_out[!grepl("RCR", datf_out$mutcode),]
datf_out <- datf_out[!grepl("^A, |, A", datf_out$mutcode),]
datf_out <- datf_out[!grepl("^SF, |, SF", datf_out$mutcode),]
datf_out <- datf_out[!grepl("^DI, |, DI", datf_out$mutcode),]
datf_out <- datf_out[!grepl("^IE, |, IE", datf_out$mutcode),]

fname <- paste0("hiv-", nm, ".tsv")
fpath <- file.path(hiv_case_dir, fname)
write.table(datf_out, fpath, sep = '\t', row.names = F, quote = F)
