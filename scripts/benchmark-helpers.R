library(scriptloc)
library(ROCR)
source(file.path(script_dir_get(), "plotting-helpers.R"))
fn_gen <- function(fn) {
    function(x) {
        stopifnot(ncol(x) == 2)
        a <- x[,1]
        b <- x[,2]
        signif(fn(a, b), 3)
    }
}

MSE    <- fn_gen(function(a, b) sqrt(mean((a - b) ^ 2)))
SRCC   <- fn_gen(function(a, b) cor(a, b, method = 'spearman'))
AUC    <- fn_gen(function(a, b) performance(prediction(b, a), "auc")@y.values[[1]])
PRAUC  <- fn_gen(function(a, b) performance(prediction(b, a), "aucpr")@y.values[[1]])
BCE    <- fn_gen(function(a, b) {
    lab  <- a
    nlab <- 1 - lab
    pred <- b
    pred[pred >= 1] <- 0.9999
    pred[pred <= 0] <- 0.0001
    lpred  <- log(pred)
    neg_lpred <-  log(1 - pred)
    v <- (lab * lpred) + (nlab * neg_lpred)
    mean(v * (-1))
})

perf_fngen <- function(things) {
    function(x) {
        template <- "%s(x)"
        setNames(sapply(things, function(thing) {
            eval(parse(text = sprintf(template, thing)))
        }), things)
    }
}
regr_perf <- perf_fngen(c("MSE", "SRCC"))
bind_perf <- perf_fngen(c("BCE", "AUC", "PRAUC"))



fs_datf_fetch <- function(d, f_types = NULL) {
    fs    <- list.files(d, full.names = T)
    fbase <- basename(fs)
    sp    <- strsplit(fbase, "[.]")
    fs_datf <- data.frame(
        benchmark = sapply(sp, `[[`, 1),
        method    = sapply(sp, `[[`, 2),
        f_ext     = sapply(sp, `[[`, 3),
        file      = fs
    )
    fs_datf$f_type <- sapply(fs_datf$f_ext, function(xt) {
        if (xt %in% c("csv", "tsv")) {
            "results"
        } else if (xt == "time") {
            "time"
        } else {
            stop("Unknown extension")
        }
    })
    if (is.null(f_types)) {
        fs_datf
    } else {
        subset(fs_datf, f_type %in% f_types)
    }
}

bench_results_load <- function(fs_datf) {
    avail_methods <- c("deephlaffy", "mhcflurry", "netmhcpan",
                       "transphla", "mixmhcpred", "tinyhlanet")
    stopifnot(c("method", "file", "benchmark") %in% colnames(fs_datf))
    stopifnot(all(fs_datf$method %in% avail_methods))
    stopifnot(!anyDuplicated(fs_datf$method))
    res <- apply(fs_datf, 1, function(x) {
        readfn <- eval(parse(text = paste0(x[["method"]], "_read")))
        y           <- readfn(x[["file"]])
        y$method    <- x[["method"]]
        y$benchmark <- x[["benchmark"]]
        y
    }, simplify = F)
    names(res) <- fs_datf$method
    res
}


# The functions to read pHLA-I predictions are expected to have the following
# columns:
# 1. allele     - HLA allele
# 2. peptide    - Peptide sequence
# 3. regressand - Binding affinity
# 4. binder     - Presentation probability
#
# If [3] & [4] do not exist, they will be substituted based on other data.

deephlaffy_read <- function(f) read.table(f, sep = '\t', header = T, stringsAsFactors = F)
tinyhlanet_read <- deephlaffy_read


mhcflurry_read <- function(f) {
    y <- read.table(f, sep = ',', header = T, stringsAsFactors = F)
    y$tmp <- y$mhcflurry_affinity
    y$tmp[y$mhcflurry_affinity > 20000] <- 20000
    y$pred_regressand <- 1 - (log(y$tmp) / log(20000))
    y$pred_binder <- y$mhcflurry_presentation_score
    y$pred_bind_perc <- y$mhcflurry_presentation_percentile
    return(y)
}

netmhcpan_read <- function(f) {
    y <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)
    y$tmp <- y$affinity
    y$tmp[y$affinity > 20000] <- 20000
    y$pred_regressand <- 1 - (log(y$tmp) / log(20000))
    y$pred_binder <- y$bindscore
    y$pred_bind_perc   <- y$bindrank
    return(y)
}

transphla_read <- function(f) {
    y <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)

    # We don't have regression output, so we're using the binder probability
    y$pred_regressand <- y$y_prob
    y$pred_binder     <- y$y_prob
    colnames(y)[colnames(y) == "HLA"] <- "allele"
    return(y)
}

mixmhcpred_read <- function(f) {
    y <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)

    # We don't have regression output, so we're using the binder probability
    y$pred_regressand <- y$bindscore
    y$pred_binder     <- y$bindscore
    return(y)
}
