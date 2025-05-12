library(scriptloc)
script_dir <- script_dir_get()
source(file.path(script_dir, "benchmark-helpers.R"))
projroot   <- file.path(script_dir, '..')
res_dir    <- file.path(projroot, 'results')
bench_dir  <- file.path(res_dir, '02-benchmark')
plots_dir  <- file.path(bench_dir, 'plots')

base_d <- file.path(bench_dir, 'results')

ds <- file.path(base_d, c("train", "test", "val"))

methods <- c("deephlaffy", "mhcflurry", "netmhcpan", "transphla")

fs_lst <- local({
    tmp <- list.files(ds, full.names = T, pattern = '[.]tsv$|[.]csv$|[.]tsv.gz$')
    method <- sapply(strsplit(basename(tmp), "[.]"), `[[`, 2)
    split(tmp, method)
})
fs_lst <- fs_lst[names(fs_lst) %in% methods]


perf_lst <- lapply(names(fs_lst), function(method) {
    print(method)
    readfn <- eval(parse(text = sprintf("%s_read", method)))
    fs <- local({
        tmp <- fs_lst[[method]]
        spl <- sapply(strsplit(basename(tmp), "[.]"), `[[`, 1)
        substr(spl, 1, 1) <- toupper(substr(spl, 1, 1))
        setNames(tmp, spl)[c("Train", "Val", "Test")]
    })

    dat_lst <- lapply(fs, function(f) {
        x <- readfn(f)
        regr <- subset(x, regressand > (-1), select = c("regressand", "pred_regressand"))
        bind <- subset(x, binder > (-1), select = c("binder", "pred_binder"))
        bind$binder <- bind$binder

        # TransPHLA has missing predictions.
        bind <- bind[complete.cases(bind),]
        list(regr = regr, bind = bind)
    })

    regr_lst <- lapply(dat_lst, `[[`, "regr")
    bind_lst <- lapply(dat_lst, `[[`, "bind")

    regr_mat <- do.call("rbind", lapply(regr_lst, regr_perf))
    bind_mat <- do.call("rbind", lapply(bind_lst, bind_perf))

    mat0 <- cbind(regr_mat, bind_mat)

    worst_case <- setNames(sapply(colnames(mat0), function(meas) {
        summ_fn <- if (meas %in% c("MSE", "BCE")) {
            max
        } else if (meas %in% c("SRCC", "AUC", "PRAUC")) {
            min
        } else {
            stop()
        }
        summ_fn(mat0[,meas])
    }), colnames(mat0))

    mat1 <- rbind(mat0, WorstCase = worst_case)


    regr_full <- do.call("rbind.data.frame", lapply(dat_lst, `[[`, "regr"))
    bind_full <- do.call("rbind.data.frame", lapply(dat_lst, `[[`, "bind"))

    mat <- rbind(mat1,
                 All = c(regr_perf(regr_full), bind_perf(bind_full)))

    summ <- data.frame(
        method     = method,
        data_split = rownames(mat),
        mat,
        row.names = NULL
    )

    list(summ = summ, regr_full = regr_full, bind_full = bind_full)
})
names(perf_lst) <- names(fs_lst)

summ_lst <- lapply(perf_lst, `[[`, "summ")

summ_df <- do.call("rbind.data.frame", summ_lst)
noreg_methods <- c("transphla")
summ_df$SRCC[summ_df$method %in% noreg_methods] <- NA
summ_df$MSE[summ_df$method %in% noreg_methods] <- NA

# In our case, the "Test" set is guaranteed to be an
# independent test, but for the other methods, they
# might have some training data here. We give ourselves
# this handicap to show that our algorithm works well
# despite this.
bind_lst <- lapply(perf_lst, function(x) {
    bind <- x$bind_full
    bind[grepl("Test", rownames(bind)),]
})

plot_data <- lapply(bind_lst, function(x) {
    to_datf <- function(z, nms = NULL) {
        datf <- data.frame(
            x = z@x.values[[1]],
            y = z@y.values[[1]]
        )
        if (!is.null(nms)) colnames(datf) <- nms
        return(datf)
    }

    get_val <- function(z) z@y.values[[1]]
    pred <- prediction(x[[2]], x[[1]])
    s <- list(
        f1      = to_datf(performance(pred, "f"),
                          c("Cutoff", "F1 score")),
        roc     = to_datf(performance(pred, "tpr", "fpr"),
                          c("FPR", "TPR")),
        auc     = get_val(performance(pred, "auc")),
        pr_roc  = to_datf(performance(pred, "prec", "rec"),
                          c("Recall", "Precision")),
        prauc   = get_val(performance(pred, "aucpr"))
    )
})

if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = T)

#------------------------------------------------#
#                      ROC                       #
#------------------------------------------------#
{
roc_file <- file.path(plots_dir, 'roc-comparison.tiff')
tiff_open(roc_file, res = 550)
par(family = 'symbol')
for (method in methods) {
    plot_fn <- if (method == methods[1]) plot else points
    plot_fn(plot_data[[method]]$roc, type = 'l',
            lty = which(methods == method),
            col = method_colors(method), lwd = 2,
            xlab = "False Positive Rate",
            ylab = "True Positive Rate"
    )
}

roc_legend <- local({
    tmp    <- signif(sapply(plot_data, `[[`, "auc"), 3)
    mnames <- method_names(names(tmp))
    sprintf("%s (AUC: %0.3f)", mnames, unname(tmp))
})

legend("bottomright", bty = 'n', lty = seq_along(methods), lwd = 2,
       legend = roc_legend, col = method_colors(methods), cex = 0.75)
tiff_close(roc_file)
}


{
roc_inset_file <- file.path(plots_dir, 'roc-comparison-inset.tiff')
tiff_open(roc_inset_file, res = 550)
par(family = 'symbol')
for (method in methods) {
    plot_fn <- if (method == methods[1]) plot else points
    inset_data <- subset(plot_data[[method]]$roc,
                         (FPR > 0.01) & (FPR < 0.1))
    plot_fn(inset_data, type = 'l',
            lty = which(methods == method),
            col = method_colors(method), lwd = 2,
            xlim = c(0.01, 0.1),
            xlab = NA, ylab = NA
    )

}

tiff_close(roc_inset_file)
}

#------------------------------------------------#
#             Precision-Recall plot              #
#------------------------------------------------#
{
pr_roc_file <- file.path(plots_dir, 'pr-roc-comparison.tiff')
tiff_open(pr_roc_file, res = 550)
par(family = 'symbol')
for (method in methods) {
    plot_fn <- if (method == methods[1]) plot else points
    plot_fn(plot_data[[method]]$pr_roc, type = 'l',
            lty = which(methods == method),
            col = method_colors(method), lwd = 2,
            xlab = "Recall",
            ylab = "Precision"
    )
}

pr_roc_legend <- local({
    tmp    <- signif(sapply(plot_data, `[[`, "prauc"), 3)
    mnames <- method_names(names(tmp))
    sprintf("%s (PRAUC: %0.3f)", mnames, unname(tmp))
})

legend("bottomleft", bty = 'n', lty = seq_along(methods), lwd = 2,
       legend = pr_roc_legend, col = method_colors(methods), cex = 0.75)
tiff_close(pr_roc_file)
}

#------------------------------------------------#
#                 F1 score plot                  #
#------------------------------------------------#
{
f1_file <- file.path(plots_dir, 'f1-comparison.tiff')
tiff_open(f1_file, res = 550)
par(family = 'symbol')
for (method in methods) {
    plot_fn <- if (method == methods[1]) plot else points
    plot_fn(plot_data[[method]]$f1, type = 'l',
            lty = which(methods == method),
            col = method_colors(method), lwd = 2,
            xlab = "Cutoff",
            ylab = expression("F"[1]*" score"),
            ylim = c(0, 1)
    )
}

f1_legend <- local({
    tmp    <- signif(sapply(plot_data, function(x) max(x$f1[,"F1 score"], na.rm = T)) , 3)
    mnames <- method_names(names(tmp))
    sprintf("%s (Max F1 score: %0.3f)", mnames, unname(tmp))
})

legend("bottomleft", bty = 'n', lty = seq_along(methods), lwd = 2,
       legend = f1_legend, col = method_colors(methods), cex = 0.75)
tiff_close(f1_file)
}

speed_d <- file.path(base_d, 'speed')

speed_d2 <- file.path(speed_d, '..', '..', 'netmhcpan-parallel',
                      'output')
speed_dat2 <- local({
    tmp <- list.files(speed_d2, pattern = '[.]time$', full.names = T)
    names(tmp) <- sub("[.]time", "", basename(tmp))
    s <- strsplit(sapply(tmp, readLines), " ")
    tm <- as.numeric(sapply(s, `[[`, 1))
})

speed_mem_dat <- local({
    speed_get <- function(f) strsplit(readLines(f), "[ ]")[[1]][1]
    mem_get <- function(f) strsplit(readLines(f), "[ ]")[[1]][2]
    tmp <- list.files(speed_d, pattern = ".time$", full.names = T)
    tmp <- tmp[!grepl("csv", tmp)]
    sp <- strsplit(basename(tmp), "[.]")
    Power <- as.numeric(sub("peps[-]", "", sapply(sp, `[[`, 1)))
    Set   <- 10 ^ Power
    Method <- sapply(sp, `[[`, 2)
    speed_vec <- as.numeric(sapply(tmp, speed_get))
    mem_vec <- signif(as.numeric(sapply(tmp, mem_get)) / 1000000, 3)
    list(
        speed = tapply(speed_vec, list(Set = Set, Method = Method), identity),
        mem = tapply(mem_vec, list(Set = Set, Method = Method), identity)
    )
})

speed_dat <- cbind(speed_mem_dat$speed[,methods],
                   'netmhcpan_p' = speed_dat2)
mem_dat   <- speed_mem_dat$mem[,methods]

l10 <- function(x) log10(as.numeric(x))
l2 <- function(x) log2(as.numeric(x))
l10p1 <- function(x) log10(as.numeric(x)+1)

{
speed_file <- file.path(plots_dir, 'speed-comparison.tiff')
tiff_open(speed_file, res = 550)
par(family = 'symbol')
speed_methods <- c('deephlaffy', 'mhcflurry', 'netmhcpan',
                   'netmhcpan_p', 'transphla')

for (method in speed_methods) {
    plot_fn <- if (method == methods[1]) plot else points
    plot_fn(l10(rownames(speed_dat)),
            speed_dat[,method], type = 'b',
            lty = which(speed_methods == method),
            pch = which(speed_methods == method),
            col = method_colors(method), lwd = 2,
            xlab = "# peptides", ylab = "Time to run (s)",
            ylim = c(0, 900), xaxt = 'n')
}
axis(1, at = seq_len(nrow(speed_dat)),
     c("10", "100", "1k", "10k", "100k", "1M"))
speed_legend <- local({
    tmp    <- speed_dat["1e+06",speed_methods]
    mnames <- method_names(names(tmp))
    sprintf("%s (Max time: %0.1f s)", mnames, unname(tmp))
})
legend("topleft", lwd = 2,
       bty = 'n', lty = seq_along(speed_methods),
       pch = seq_along(speed_methods),
       legend = speed_legend,
       col = method_colors(speed_methods), cex = 0.75)
tiff_close(speed_file)
}

{
speed_file2 <- file.path(plots_dir, 'speed-comparison-2.tiff')
tiff_open(speed_file2, width = 1800, height = 2000, res = 300)
par(family = 'symbol', oma = c(6, 3, 3, 3))
speed_vec <- local({
    per_sec <- 1000000 / speed_dat["1e+06",]
    per_sec[order(per_sec)]
})
bp <- barplot(speed_vec, xaxt = 'n',
        ylab = "# predictions per second",
        col = method_colors(names(speed_vec)),
        ylim = c(0, 25000))
tmpnms <- method_names(names(speed_vec))
tmpnms[tmpnms == "netMHCpan 4.1 (par)"] <- "netMHCpan 4.1\n(par)"
text(bp[,1] - 0.6, -4300, labels = tmpnms, xpd = NA, srt = 45)
in_mille <- sprintf("~%0.1fk", unname(speed_vec) / 1000)
in_mille[in_mille == "~18.0k"] <- "~18k"
text(bp[,1], speed_vec + 750, labels = in_mille)
tiff_close(speed_file2)
}


{
mem_file <- file.path(plots_dir, 'mem-comparison.tiff')
tiff_open(mem_file, res = 550)
par(family = 'symbol')
for (method in methods) {
    plot_fn <- if (method == methods[1]) plot else points
    plot_fn(l10(rownames(mem_dat)),
            mem_dat[,method], type = 'b',
            lty = which(methods == method),
            pch = which(methods == method),
            col = method_colors(method), lwd = 2,
            xlab = "# peptides",
            ylab = expression("Max. resident memory (GB)"),
            ylim = c(0, 45), xaxt = 'n')
}
axis(1, at = seq_len(nrow(speed_dat)),
     c("10", "100", "1k", "10k", "100k", "1M"))
mem_legend <- local({
    tmp    <- mem_dat["1e+06",]
    mnames <- method_names(names(tmp))
    sprintf("%s (Max memory: %0.1f GB)", mnames, unname(tmp))
})
legend("topleft", lwd = 2,
       bty = 'n', lty = seq_along(methods),
       pch = seq_along(methods),
       legend = mem_legend,
       col = method_colors(methods), cex = 0.75)
tiff_close(mem_file)
}


test_f_get <- function(f) {
    fname <- grep("^test", basename(f), value = T)
    file.path(dirname(f[1]), fname)
}

test_bind_lst <- list(
    deephlaffy = subset(deephlaffy_read(test_f_get(fs_lst$deephlaffy)),
                        binder > (-1), select = c("hla", "binder", "pred_binder")),
    mhcflurry  = subset(mhcflurry_read(test_f_get(fs_lst$mhcflurry)),
                        binder > (-1), select = c("allele", "binder", "pred_binder"))
)

auc_datf <- as.data.frame(lapply(test_bind_lst, function(datf) {
    #### Check later
    sp <- Filter(function(x) length(unique(x[[1]])) == 2, split(datf[,c(2,3)], datf[[1]]))
    s <- sapply(sp, AUC)
}))

prauc_datf <- as.data.frame(lapply(test_bind_lst, function(datf) {
    sp <- Filter(function(x) length(unique(x[[1]])) == 2, split(datf[,c(2,3)], datf[[1]]))
    s <- sapply(sp, PRAUC)
}))

f1_datf <- as.data.frame(lapply(test_bind_lst, function(datf) {
    sp <- Filter(function(x) length(unique(x[[1]])) == 2, split(datf[,c(2,3)], datf[[1]]))
    s <- sapply(sp, function(x) {
        pred <- prediction(x[[2]], x[[1]])
        max(performance(pred, "f")@y.values[[1]], na.rm = T)
    })
    signif(s, 3)
}))

{
alauc_file <- file.path(plots_dir, 'allelewise-auc-comparison.tiff')
tiff_open(alauc_file, res = 550)
plot(deephlaffy ~ mhcflurry, auc_datf,
     xlab = "TinyHLAnet AUC", ylab = "MHCflurry 2 AUC", pch = 21,
     bg = adjustcolor('#6eaf5e', alpha.f = 0.5),
     xlim = c(0.95, 1), ylim = c(0.95, 1)
)
abline(coef = c(0, 1))
tiff_close(alauc_file)
}

{
alprauc_file <- file.path(plots_dir, 'allelewise-prauc-comparison.tiff')
tiff_open(alprauc_file, res = 550)
plot(deephlaffy ~ mhcflurry, prauc_datf,
     xlab = "TinyHLAnet PRAUC", ylab = "MHCflurry 2 PRAUC", pch = 21,
     bg = adjustcolor("#8981a7", alpha.f = 0.5),
     xlim = c(0.9, 1), ylim = c(0.9, 1)
)
abline(coef = c(0, 1))
tiff_close(alprauc_file)
}

{
f1_file <- file.path(plots_dir, 'allelewise-f1-comparison.tiff')
tiff_open(f1_file, res = 550)
plot(deephlaffy ~ mhcflurry, f1_datf,
     xlab = expression("TinyHLAnet max "*"F"[1]*" score"),
     ylab = expression("MHCflurry 2 max "*"F"[1]*" score"),
     pch = 21, bg = adjustcolor('#bf8502', alpha.f = 0.5),
     xlim = c(0.85, 1), ylim = c(0.85, 1))
abline(coef = c(0, 1))
tiff_close(f1_file)
}
