cropfn <- function(f) {
    cmd <- sprintf("mogrify -trim %s", f)
    system(cmd)
    f
}

library(scriptloc)
SCRIPT_DIR   <- script_dir_get()
PROJROOT     <- file.path(SCRIPT_DIR, "..")

case_dir     <- file.path(PROJROOT, 'results', '03-case-study')
f <- file.path(case_dir, '1-s2.0-S2211124723000062-mmc2.xlsx')


case_data_dir <- file.path(case_dir, 'data')
epi_f <- file.path(case_data_dir, 'sars-cov-2.epitopes')

epi <- subset(read.table(epi_f, sep = '\t', header = T),
              nchar(epitope) == 9)

escan_dir <- file.path(case_dir, 'epitope-scan')
d  <- file.path(escan_dir, 'epitopes', 'aff0_bind0')
fs <- local({
    tmp <- list.files(d, full.names = T, pattern = '.tsv.gz')
    atmp0 <- sub("[.]tsv[.]gz", "", basename(tmp))
    al    <- sapply(strsplit(atmp0, "[-]"), function(x)paste0(x[1], "-", x[2], "*", x[3], ":", x[4]))
    setNames(tmp, al)
})

pred <- do.call("rbind.data.frame", lapply(names(fs), function(al) {
    x <- read.table(fs[[al]], sep = '\t', header = T, stringsAsFactors = F)
    x$allele <- al
    x
}))


alleles <- unique(epi$allele)

res_lst <- lapply(alleles, function(al) {
    print(al)
    al_epi  <- subset(epi, allele == al)

    al_pred <- local({
        tmp0 <- subset(pred, allele == al)
        tmp1 <- tmp0[order(tmp0$rank),]
        tmp2 <- tmp1[!duplicated(tmp1$epitope),]
        stopifnot(all(al_epi$epitope %in% tmp2$epitope))
        tmp3 <- tmp2[tmp2$epitope %in% al_epi$epitope,]
    })
})

res <- do.call("rbind.data.frame", res_lst)

res$detected <- res$bindprob > 0.5

detvec <- local({
    tmp <- c(sapply(split(res$detected, res$allele), mean),
             "Combined" = mean(res$detected))
    round(tmp, 2)
})


split_vec <- function(v, breaks, last_inclusive = F, fmt = "%0.2f") {
    stopifnot(length(breaks) >= 2)

    iter <- cbind(
        start = breaks[seq_len(length(breaks) - 1)],
        end   = breaks[seq_len(length(breaks))[-1]]
    )

    start_nm <- sprintf(fmt, iter[,"start"])
    end_nm   <- sprintf(fmt, iter[,"end"])
    last_nm  <- sprintf(fmt, breaks[length(breaks)])
    bnames <- c(paste0(start_nm, "-", end_nm),
                paste0(last_nm, "+"))

    res <- setNames(rep(NA, length(breaks)), bnames)
    for (k in seq_len(nrow(iter))) {
        startv <- iter[k, "start"]
        endv   <- iter[k, "end"]
        l      <- (v >= startv) & (v < endv)
        res[k] <- sum(l)
    }
    last_start <- iter[k, "end"]
    if (last_inclusive) {
        res <- res[1:k]
        res[k] <- res[k] + sum(v >= last_start)
    } else {
        res[k+1] <- sum(v >= last_start)
        if (res[k+1] == 0) res <- res[1:k]
    }
    res
}

prob_bar <- split_vec(res$bindprob, seq(0, 1, 0.05), T)

rank_pie <- split_vec(res$rank, c(1, 5, 10, 30, 50), F, "%d")
names(rank_pie)[1] <- "Top 5"
names(rank_pie) <- paste0(names(rank_pie), " (", rank_pie, ")")

rankp_pie <- split_vec(res$rankp, c(0, 1, 2, 3, 4, 5), F, "%d")
names(rankp_pie) <- paste0(names(rankp_pie), "%")
names(rankp_pie)[1] <- "Top 1%"
names(rankp_pie) <- paste0(names(rankp_pie), " (", rankp_pie, ")")


cols <- c("#FB7979",
          "#F7F4A6",
          "#74CE98",
          "#00BFFF",
          "white")

plot_f <- file.path(escan_dir, 'epitope-scan.png')
{
png(plot_f, width = 2000, height = 2000, res = 300)
par(family = 'symbol', oma = c(0, 0, 0, 0), mar = c(1, 4, 5, 2))
layout(rbind(c(1, 1), c(2, 3)))
bp <- barplot(
    prob_bar, las = 2, xaxt = 'n',
    ylab = "# epitopes"
)
legend("topleft", bty = 'n', cex = 1.25,
       legend = c("Predicted binding probability distribution",
                  sprintf("Median: %0.3f", median(res$bindprob))))
text(x = bp[,1] - 0.4, y = -1.5, labels = names(prob_bar),
     srt = 45, xpd = NA, cex = 0.8)
pie(rank_pie,  col = cols)
mtext("Intraprotein rank of epitopes", side = 3)
pie(rankp_pie, col = cols)
mtext("Intraprotein rank (%) of epitopes", side = 3)
dev.off()
cropfn(plot_f)
}

int_d <- file.path(case_dir, 'mutant-interactions')
mat_d <- file.path(int_d, 'binding-mats')

mat_fs <- local({
    tmp <- list.files(mat_d, full.names = T)
    setNames(tmp, tools::file_path_sans_ext(basename(tmp)))
})

readmat <- function(f) {
    AA <- function(x) sapply(strsplit(x, "_"), `[[`, 2)
    tmp <- read.table(f, sep = '\t', row.names = 1, header = T,
                      check.names = F)

    m   <- as.matrix(tmp)
    colnames(m) <- AA(colnames(m))
    rownames(m) <- AA(rownames(m))
    t(m)
}


cullmat <- function(m, cutoff = 0.01) {
    m[abs(m) < cutoff] <- NA
    m
}

gmat_f <- mat_fs[['HLA-A-02-01-YLQPRTFLL']]
bmat_f <- mat_fs[['HLA-A-02-01-YFQPRTFLL']]

library(pheatmap)
library(RColorBrewer)


rdylbu_r <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
ph <- local({
    br <- seq(-0.1, 0.1, 0.005)
    tmp <- pheatmap
    formals(tmp)$cluster_rows <- F
    formals(tmp)$cluster_cols <- F
    formals(tmp)$cellwidth    <- 10
    formals(tmp)$cellheight   <- 10
    formals(tmp)$fontfamily   <- "symbol"
    formals(tmp)$legend       <- F
    formals(tmp)$na_col       <- "white"
    formals(tmp)$breaks       <- br
    formals(tmp)$color        <- rdylbu_r(length(br))
    formals(tmp)$angle_col    <- 0
    formals(tmp)$border_color <- "black"
    tmp
})
g0 <- readmat(gmat_f)
b0 <- readmat(bmat_f)
dif0 <- b0 - g0

cutoff <- 0.01
g <- cullmat(g0, cutoff)
b <- cullmat(b0, cutoff)
matl <- (!is.na(g)) & (!is.na(b))
dif <- local({
    tmp <- dif0
    tmp[!matl] <- NA
    ind <- rownames(g0) != rownames(b0)
    stopifnot(sum(ind) > 0)
    nm <- rownames(g0)
    nm[ind] <- paste0(rownames(g0)[ind], "→", rownames(b0)[ind])
    rownames(tmp) <- nm
    tmp
})


gplot_f  <- file.path(int_d, 'good.png')
bplot_f  <- file.path(int_d, 'bad.png')
dplot_f  <- file.path(int_d, 'diff.png')
dplot_f2 <- file.path(int_d, 'diff2.png')

{
png(gplot_f, width = 2000, height = 2000, res = 300)
ph(g)
dev.off()
cropfn(gplot_f)
}

{
png(bplot_f, width = 2000, height = 2000, res = 300)
ph(b)
dev.off()
cropfn(bplot_f)
}

{
png(dplot_f, width = 2000, height = 2000, res = 300)
ph(dif)
dev.off()
cropfn(dplot_f)
}

{
png(dplot_f2, width = 3000, height = 2000, res = 300)
ph(dif, legend = T)
dev.off()
cropfn(dplot_f2)
}
