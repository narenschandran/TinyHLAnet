# This script generates train-val splits from the training datasets
# and structures these splits along with the test split in a location
# that is convenient to access for the subsequent scripts.

library(scriptloc)
library(tools)
script_dir    <- scriptloc::script_dir_get()
projroot      <- file.path(script_dir, "..")
prereq_dir    <- file.path(projroot, "prereq")
datasets_dir  <- file.path(prereq_dir, "datasets")


#-----------------------------------------------------------------#
#                      Function definitions                       #
#-----------------------------------------------------------------#

#' Order datapoints by allele and peptide
#'
#' Input : Table with pHLA-I data points
#' Output: Table with pHLA-I data points ordered by allele
#'         and peptide
order_data <- function(x) {
    # We order the data by allele and peptide key
    # to ensure that even if the input data is
    # shuffled, we end up with the same order before
    # we do our own shuffling.
    stopifnot(all(c("pmhc_key", "regressand", "binder") %in%
                  colnames(x)))
    ind <- order(x$pmhc_key, x$regressand, x$binder)
    tmp <- x[ind,]
}

#' Convenient function to generate train-test splits
#'
#' Input: Table with pHLA-I data points
#' Output: A list with two tables where the first table has
#'         a proprotion of datapoints (specified by `prop`), with
#'         the second table having the rest. The data points
#'         are split randomly.
split_data <- function(x, prop, seed, split_names = NULL) {
    stopifnot((prop >= 0) & (prop <= 1))
    y        <- order_data(x)

    stopifnot("pmhc_key" %in% colnames(y))
    keys <- sort(unique(y$pmhc_key))

    N        <- length(keys)
    sp1_size <- ceiling(prop * N)

    set.seed(seed)
    key_sp <- within(
        list(sp1 = sort(sample(keys, sp1_size, replace = F))), {
        sp2 <- sort(setdiff(keys, sp1))
    })
    intr <- Reduce(intersect, key_sp)
    stopifnot(length(intr) == 0)
    sp <- lapply(key_sp, function(k) subset(y, pmhc_key %in% k))
    if (!is.null(split_names)) names(sp) <- split_names

    return(sp)
}


#' Same as `split_data`, but splits occur separately according to
#' data type (qualitative and quantitative) to ensure uniformity
split_within_data_type <- function(x, prop, seed, split_names = NULL) {
    dat_lst <- split(x, ifelse(x$regressand > (-1),
                               "reg", "bind"))
    sp_lst <- lapply(dat_lst, split_data, prop = prop, seed = seed)

    sp <- list(
        sp1 = do.call('rbind.data.frame', lapply(sp_lst, `[[`, 1)),
        sp2 = do.call('rbind.data.frame', lapply(sp_lst, `[[`, 2))
    )
    if (!is.null(split_names)) names(sp) <- split_names
    lapply(sp, order_data)
}


#' Simple wrapper function for reading files
readf <- function(f) read.table(f, sep = '\t', header = T,
                                stringsAsFactors = F)


#' Simple wrapper function for writing files
writef <- function(dat, f) {
    cols <- c("allele", "peptide", "mhcpocket", "regressand",
              "binder", "pmhc_key")
    stopifnot(all(cols %in% colnames(dat)))
    fcon <- xzfile(f, 'wb', compression = 9)
    out_dat <- order_data(dat[,cols])
    write.table(out_dat, fcon, sep = '\t', row.names = F, quote = F)
    close(fcon)
    f
}

#-----------------------------------------------------------------#

# Get list of dataset files
fs <- local({
    tmp <- list.files(datasets_dir, full.names = T)
    setNames(tmp, basename(file_path_sans_ext(tmp, compression = T)))
})

odir <- file.path(projroot, "datasets", "raw")
if (!dir.exists(odir)) dir.create(odir, recursive = T)

# Take the test data verbatim
test_data <- readf(fs[["test"]])

# Split the training data into training and validation splits
x <- readf(fs[["train"]])
y <- split_within_data_type(x, 0.9, 1, c('train', 'val'))
train_data0 <- y$train

# Some alleles have only binders and no non-binders (possibly)
# because of removing peptide linkage across datasets. We're
# removing this from the training data
bx <- subset(train_data0, binder > (-1))

bsp <- split(bx$binder, bx$allele)
bmat <- do.call("rbind", lapply(bsp, function(x) {
    l <- x == 1
    c("nonbinders" = sum(!l), "binders" = sum(l))
}))

bdatf0 <- data.frame(
    allele = names(bsp),
    bmat, row.names = NULL
)

ncutoff <- 5
selected_alleles <- bdatf0[apply(bdatf0[,-1], 1, min) >= ncutoff,1]
excluded_alleles <- bdatf0[apply(bdatf0[,-1], 1, min) < ncutoff,1]

# We keep all regression data, but remove biased alleles
train_data <- local({
    tmp <- train_data0
    l  <- (tmp$regressand > (-1)) |
          ((tmp$binder > (-1)) & (tmp$allele %in% selected_alleles))
    tmp[l,]
})

val_data   <- y$val
writef(train_data, file.path(odir, 'train.tsv.xz'))
writef(val_data, file.path(odir, 'val.tsv.xz'))
writef(test_data, file.path(odir, 'test.tsv.xz'))

exl_data <- local({
    subset(train_data0,
           (binder == 1) &
           (allele %in% excluded_alleles))
})

odir2 <- file.path(projroot, 'results', '01-model-tuning', '04-unseen-alleles')

if (!dir.exists(odir2)) dir.create(odir2, recursive = T)

writef(exl_data, file.path(odir2, 'excluded.tsv.xz'))
