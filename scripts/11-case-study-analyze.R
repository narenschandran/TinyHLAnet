library(scriptloc)
library(pheatmap)
library(RColorBrewer)

rdylbu_r <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))

get_aa <- function(x) sapply(strsplit(x, "_"), `[[`, 2)
ph <- local({
    br <- seq(-0.2, 0.2, 0.05)
    tmp <- pheatmap
    formals(tmp)$cluster_rows    <- F
    formals(tmp)$cluster_cols    <- F
    formals(tmp)$display_numbers <- F
    formals(tmp)$breaks <- br
    formals(tmp)$color  <- rdylbu_r(length(br))
    formals(tmp)$border_color  <- "black"
    formals(tmp)$na_col <- "white"
    formals(tmp)$fontfamily <- "symbol"
    formals(tmp)$legend <- T
    formals(tmp)$cellwidth  <- 10
    formals(tmp)$cellheight <- 10
    formals(tmp)$angle_col  <- 0
    function(mat, main = NA, ...) {
        num_col <- matrix("black", nrow = nrow(mat),
                          ncol = ncol(mat))
        num_col[abs(mat) > 0.1] <- "white"
        num_col[is.na(mat)] <- "white"
        tmp(mat, main = main, number_color = num_col, ...)
    }
})

read_mat <- function(f) {
    x <- read.table(f, sep = '\t', header = T, row.names = 1,
                    stringsAsFactors = F)
    m <- as.matrix(x)
    colnames(m) <- get_aa(colnames(m))
    rownames(m) <- get_aa(rownames(m))
    return(m)
}


script_dir     <- script_dir_get()
projroot       <- file.path(script_dir, '..')
case_study_dir <- file.path(projroot, 'results', '03-case-study')


f <- file.path(case_study_dir, 'case-study.deephlaffy.tsv')
x <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)
x$pmhc_key <- with(x, paste(hla, peptide))
x$fname <- paste0(gsub("[*: ]", "-", x$pmhc_key), '.tsv')
x$fpath <- file.path(case_study_dir, 'binding-mats', x$fname)

xsp <- split(x, x$set)

dif <- sapply(xsp, function(y) {
    diff(y$pred_regressand)
})

ysp <- xsp[order(dif)]

get_difmat <- function(datf) {
    zlst <- lapply(datf$fpath, read_mat)
    stopifnot(all.equal(rownames(zlst[[1]]), rownames(zlst[[2]])))

    al <- unique(datf$hla)
    stopifnot(length(al) == 1)

    m <- as.matrix(zlst[[2]] - zlst[[1]])
    rownames(m) <- rownames(zlst[[1]])

    peps <- sapply(zlst, colnames, simplify = F)
    pep1 <- peps[[1]]
    pep2 <- peps[[2]]

    newpep <- pep1
    ch_ind <- which(pep1 != pep2)
    newpep[ch_ind] <- paste0(pep1[ch_ind], '→', pep2[ch_ind])
    colnames(m) <- newpep

    mats <- list(
        difmat = m,
        mat1   = zlst[[1]],
        mat2   = zlst[[2]]
    )
    res <- lapply(mats, function(mat) {
        mat[abs(mat) < 0.05] <- NA
        mat
    })
    res$allele <- al
    res
}

difmat_lst <- lapply(ysp, get_difmat)



formals(tiff)$width  <- 2500
formals(tiff)$height <- 1500
formals(tiff)$res    <- 300
formals(tiff)$compression <- "lzw"
formals(ph)$main <- "HLA-A*02:01"


base_odir <- file.path(case_study_dir, 'interactions')
for (nm in names(difmat_lst)) {
    odir <- file.path(base_odir, nm)
    if (!dir.exists(odir)) dir.create(odir, recursive = T)

    fbase <- file.path(odir, nm)
    dat <- difmat_lst[[nm]]

    dif_file <- paste0(fbase, '-diff.tiff')
    tiff(dif_file)
    ph(t(dat$difmat))
    dev.off()

    good_file <- paste0(fbase, '-good.tiff')
    tiff(good_file)
    ph(t(dat$mat1))
    dev.off()

    bad_file <- paste0(fbase, '-bad.tiff')
    tiff(bad_file)
    ph(t(dat$mat2))
    dev.off()
}
