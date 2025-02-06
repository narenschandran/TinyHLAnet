# This script retrieves the performance information from the entireity of
# the model tuning process.

#------------------------------------------------------------------------------#
#                               Helper functions                               #
#------------------------------------------------------------------------------#

#' Simple reader function
readf <- function(f) read.table(f, sep = '\t', header = T, check.names = T,
                                stringsAsFactors = F)

#' Wrapper to read performance table
perf_read <- function(f) {
    if (is.null(f)) {
        print("Encountered NULL file")
        return(NULL)
    }
    data.frame(readf(f), file = f)
}

#' Wrapper to read model configuration file
conf_read <- function(f) {
    if (is.null(f)) {
        print("Encountered NULL file")
        return(NULL)
    }
    tmp0 <- readLines(f)
    sp   <- strsplit(tmp0, "\t")
    tmp  <- setNames(sapply(sp, `[[`, 2),
                     sapply(sp, `[[`, 1))
    as.list(tmp)
}

#' Wrapper to report model training time in seconds
time_read <- function(f) {
    rm_pat <- "^Model training time: | seconds$"
    as.numeric(gsub(rm_pat, "", readLines(f)))
}

#' Generate a data.frame chronicling the improvement of total model loss over time
chkpt_evolution <- function(chkpt_d) {
    fs <- list.files(chkpt_d, full.names = T)
    evol_dat <- local({
        tmp  <- sub("[.]keras$", "", basename(fs))
        sp   <- strsplit(tmp, "[-]")
        datf <- data.frame(
            epoch = as.numeric(sapply(sp, `[[`, 1)),
            loss  = as.numeric(sapply(sp, `[[`, 2)),
            stringsAsFactors = F
        )
        datf[order(datf$epoch),]
    })
}

#' Retrieve model information (if any present) by recursively descending
#' the directory containing model runs
dat_get <- function(d) {
    # Create a record of the current directory along with an easy-to-refer name
    # based on the filename
    curr_fs <- local({
        tmp <- list.files(d, full.names = T)
        setNames(tmp, tools::file_path_sans_ext(basename(tmp)))
    })

    # The next chunk scans any recursively scans any directories found in the
    # current run. If it detects model configuration or performance files,
    # it loads them. Otherwise, it returns the file as is.
    dat <- lapply(curr_fs, function(f) {
        if (dir.exists(f) && (basename(f) == "chkpt")) {
            chkpt_evolution(f)
        } else if (dir.exists(f)) {
            dat_get(f)
        } else if (basename(f) == "perf.tsv") {
            perf_read(f)
        } else if (basename(f) == "conf.dct") {
            conf_read(f)
        } else if (basename(f) == "time.txt") {
            time_read(f)
        } else {
            f
        }
    })
}

#' Accumulate results in a recursive fashion from the output of `dat_get`.
#'
#' This function descends the output of `dat_get` recursively and retrives
#' them in a flat list.
acc_res <- function(dlst, acc = c()) {
    nms <- names(dlst)
    res <- list()
    # Each time we descend a level, we accumulate that level's name. We
    # expect the output to be organized in a specific way, where the innermost
    # directory has seed information, the next level has model name, the
    # next level has data type, etc. Depending on the number of descents,
    # we dynamically add these other labels.
    if (all(c("perf", "conf", "time", "chkpt") %in% nms)) {
        hierarchy <- c("expt", "data_type", "model_key", "seed_id")
        acc_names <- rev(rev(hierarchy)[seq_len(length(acc))])
        acc <- setNames(acc, acc_names)

        perf <- dlst$perf
        perf$time <- dlst$time
        for (acc_nm in acc_names) {
            perf[[acc_nm]] <- acc[[acc_nm]]
        }
        conf <- dlst$conf
        if ("model_key" %in% acc_names) {
            conf$model_key <- acc[["model_key"]]
        }

        chkpt <- dlst$chkpt

        summ <- local({
            hyperp_items <- c("pp_embdim", "contacts", "pos_conf", "effects_conf")
            info_items   <- c("model_key", "data_type", "time", "seed_id")
            tmp <- unique(perf[,info_items])
            stopifnot(rownames(tmp) == 1)
            perf_stats <- as.list(setNames(perf$test, perf$measure))
            model_info <- as.list(tmp)
            hyperp     <- conf[hyperp_items]
            c(model_info, perf_stats, hyperp)
        })
        summ$model_id <- paste0(summ$model_key, "-", summ$data_type, "-", sub("^seed[-]", "", summ$seed_id))
        summ$epoch    <- max(chkpt$epoch)

        return(list(list(perf = perf, conf = conf, chkpt = chkpt, summ = summ)))
    } else {
        for (nm in nms) {
            thing <- dlst[[nm]]
            if (is.list(thing)) {
                res <- c(res, acc_res(thing, c(acc, nm)))
            }
        }
    }
    return(res)
}

#------------------------------------------------------------------------------#

library(scriptloc)
script_dir <- script_dir_get()
projroot   <- file.path(script_dir, '..')
base_d     <- file.path(projroot, 'models', 'deephlaffy')
res_dir    <- file.path(projroot, 'results')
tune_dir   <- file.path(res_dir, '01-model-tuning')

dat_tree <- dat_get(base_d)
dat      <- local({
    tmp  <- acc_res(dat_tree)
    ids  <- lapply(tmp, function(x) x$summ$model_id)
    setNames(tmp, ids)
})

if (!dir.exists(tune_dir)) dir.create(tune_dir, recursive = T)
saveRDS(dat, file.path(tune_dir, 'data.rds'))
