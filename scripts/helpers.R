library(RColorBrewer)
library(pheatmap)

#------------------------------------------------------------------------------#
#                                Image handling                                #
#------------------------------------------------------------------------------#

#' In-place image autocropping
#' Inputs: Image file path
cropfn <- function(f) {
    fmt <- "mogrify -compress lzw -trim %s"
    cmd <- sprintf(fmt, f)
    system(cmd)
    f
}

#' Wrapper function with sane defaults to enform conformity
#' in the images I generate.
tiff_open <- function(fpath, width = 2500, height = 2500, res = 600,
                      compression = 'lzw') {
    tiff(fpath, compression = compression,
         width = width, height = height, res = res)
}

#' Pair to the "tiff_open" function. Comes with easy option to toggle cropping.
tiff_close <- function(fpath, crop = T) {
    dev.off()
    if (crop) cropfn(fpath)
    fpath
}
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#             Parse and plot for hyperparameter choice exploration             #
#------------------------------------------------------------------------------#

# Functions and associated kit to help clean up names
actv_map <- c(
    "None"    = "None",
    "sigmoid" = "Sigmoid",
    "linear"  = "Linear",
    "relu"    = "ReLU",
    "gelu"    = "GELU",
    "tanh"    = "Tanh"
)
contacts_map <- c(
    "allpairs" = "All",
    "simple"   = "Observed",
    "expdecay" = "Observed\n(Decay)"
)

#' Baseline function to map names from one vector based on a guide
annotate_fn <- function(x, guide) {
    stopifnot(all(x %in% names(guide)))
    guide[x]
}

#' Clean contact hyperpameter names
annotate_contacts <- function(x) annotate_fn(x, contacts_map)

#' Clean activation function hypermarater
annotate_actv     <- function(x) annotate_fn(x, actv_map)

#' Clean embedding dimensions and ANN node numbers
dim2str <- function(x, digits = 3) {
    fmt <- paste0("%0", digits, "d")
    sprintf(fmt, as.numeric(x))
}

#' Get only the numberic characters from string values
num_get  <- function(x) gsub("[A-Za-z]", "", x)

#' Get only non-numeric characters from string values
char_get <- function(x) gsub("[0-9]", "", x)


#' Parse ANN string configuration into the individual hyperparameters
conf_parse <- function(x, sufx) {
    cols  <- c("EmbDim", "Node", "ActvIn", "ActvOut")
    n     <- length(x)
    dummy <- rep(NA, n)

    datf <- local({
        tmp <- do.call("cbind.data.frame",
                       lapply(cols, function(x) dummy))
        setNames(tmp, cols)
    })

    procfn <- function(sp, none_val, comp, parse = "all") {
        parsefn <- if (parse == "num") num_get else if (parse == "char") char_get else identity
        parsefn(sapply(sp, function(x) {
            if ((length(x) == 1) && (x == "None")) {
                none_val
            } else if (length(x) == 1) {
                stop("Should be None or multisplit")
            } else {
                x[[comp]]
            }
        }))
    }

    sp <- strsplit(x, "[-_]")
    datf$EmbDim  <- dim2str(procfn(sp, 0, 1, "num"))
    datf$Node    <- dim2str(procfn(sp, 0, 2, "num"))
    datf$ActvIn  <- annotate_actv(procfn(sp, "None", 2, "char"))
    datf$ActvOut <- annotate_actv(procfn(sp, "None", 3, "char"))

    if (!is.null(sufx)) {
        colnames(datf) <- paste0(sufx, colnames(datf))
    }
    datf
}

pos_conf_parse  <- function(x) conf_parse(x, "Env")
fx_conf_parse   <- function(x) conf_parse(x, "Fx")

#' In case of missing values, we put the worst-possible value that
#' can be achieved for each measure.
deal_with_nan <- function(perf) {
    worst_case <-  list(
        # Quantitative
        "SRCC"  = 0  , "MSE"   = Inf,

        # Qualitative
        "BCE"   = Inf, "AUC"   = 0.5, "PRAUC" = 0
    )
    for (meas in names(worst_case)) {
        perf[[meas]][is.nan(perf[[meas]])] <- worst_case[[meas]]
    }
    perf
}


#' Rank the models by their average rank across all performance measures
rank_perf_ <- function(x, measures, select_best) {
    nms <- paste0(measures, "_rank")
    rank_lst <- setNames(lapply(measures, function(meas) {

        # We ensure that SRCC, AUC & PRAUC are given lower ranks
        # if the values are higher.
        rankfn <- if (meas %in% c("SRCC", "AUC", "PRAUC")) function(x) rank(-x) else rank
        rankfn(x[[meas]])
    }), nms)

    # Get the average rank across all the measures.
    datf <- local({
        tmp <- do.call("cbind.data.frame", rank_lst)
        tmp$rank <- rowMeans(tmp)
        tmp
    })

    # Order by ranks
    y <- local({
        tmp <- x
        tmp[,c(nms, "rank")] <- datf
        tmp[order(tmp$rank),]
    })

    # If required, select only the best model instance for each
    # hyperparameter set.
    z <- if (select_best) {
        stopifnot(all(c("model_key", "rank") %in% colnames(y)))
        y[!duplicated(y$model_key),]
    } else {
        y
    }
}

reg_meas  <- c("SRCC", "MSE")
bind_meas <- c("BCE", "AUC", "PRAUC")
all_meas  <- c(reg_meas, bind_meas)

rank_regperf  <- function(x) rank_perf_(x, reg_meas , T)
rank_bindperf <- function(x) rank_perf_(x, bind_meas, T)
rank_perf     <- function(x) rank_perf_(x, all_meas , T)


#' Generate summary matrix
mat_gen <- function(datf, split_cols, data_col = "SRCC", fn = "identity") {
    splitter <- if (is.character(split_cols)) {
        datf[,split_cols]
    } else {
        split_cols
    }
    out  <- tapply(datf[[data_col]], splitter, fn)
    return(out)
}

gen_fn <- function(data_col) {
    tmp <- mat_gen
    formals(tmp)$data_col <- data_col
    return(tmp)
}

srcc_gen  <- gen_fn("SRCC")
mse_gen   <- gen_fn("MSE")
auc_gen   <- gen_fn("AUC")
bce_gen   <- gen_fn("BCE")
prauc_gen <- gen_fn("PRAUC")
rank_gen  <- gen_fn("rank")

ph <- local({
    tmp <- pheatmap
    formals(tmp)$cluster_rows    <- F
    formals(tmp)$cluster_cols    <- F
    formals(tmp)$display_numbers <- T
    formals(tmp)$number_format   <- "%0.3f"
    formals(tmp)$cellheight      <- 30
    formals(tmp)$cellwidth       <- 30
    formals(tmp)$fontsize_number <- 10
    formals(tmp)$border_color    <- "black"
    formals(tmp)$font_family     <- "symbol"
    formals(tmp)$angle_col       <- 45
    tmp
})

reg_pal1  <- colorRampPalette(c("coral1", "lightyellow"))
reg_pal2  <- colorRampPalette(c("cornflowerblue", "lightyellow"))

bind_pal1  <- colorRampPalette(c("turquoise4", "lightyellow"))
bind_pal2 <- colorRampPalette(c("seagreen4", "lightyellow"))
bind_pal3 <- colorRampPalette(c("#d986ec", "lightyellow"))


hyperp_ph <- function(x, pal, best_fn, breaks = NA, revpal = F, best_ncol = "gold", def_ncol = "black") {
     # This prevents issues if there is only one value
    uvals <- sort(unique(c(x)))
    best_ncol <- if (length(uvals) == 1) "black" else best_ncol
    breaks <- local({
        tmpbr0 <- uvals
        if (length(tmpbr0) == 1) {
            tmpbr <- c(min(tmpbr0) - c(0.001), tmpbr0)
            #tmpbr[-length(tmpbr)]
        } else {
            breaks
        }
    })
    num_col <- matrix("black", nrow = nrow(x), ncol = ncol(x))
    num_col[x == best_fn(x)] <- best_ncol



    cl <- if ((length(breaks) == 1) && is.na(breaks)) pal(100) else pal(length(breaks))
    if (revpal) cl <- rev(cl)
    ph(x, color = cl, breaks = breaks, legend = F, number_col = num_col)
}

hyperp_phgen <- function(best_fn, revpal = F) {
    tmp <- hyperp_ph
    formals(tmp)$best_fn <- best_fn
    formals(tmp)$revpal  <- revpal
    tmp
}

mxfn <- function(x) max(x, na.rm = T)
mnfn <- function(x) min(x, na.rm = T)

srcc_hyperp_ph   <- hyperp_phgen(mxfn, T)
mse_hyperp_ph    <- hyperp_phgen(mnfn)
bce_hyperp_ph    <- hyperp_phgen(mnfn)
auc_hyperp_ph    <- hyperp_phgen(mxfn, T)
prauc_hyperp_ph  <- hyperp_phgen(mxfn, T)

#------------------------------------------------------------------------------#
