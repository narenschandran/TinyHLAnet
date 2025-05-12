library(Boruta)
library(RColorBrewer)
library(pheatmap)
library(scriptloc)
library(dendextend)
library(writexl)
script_dir <- script_dir_get()
projroot   <- file.path(script_dir, '..')
res_dir    <- file.path(projroot, 'results', '04-dissection')
prereq_dir <- file.path(projroot, 'prereq')
source(file.path(script_dir, 'plotting-helpers.R'))


#------------------------------------------------#
#            Peptide processing score            #
#------------------------------------------------#

tap_read <- function(f) {
    x <- read.table(f, sep = '\t', header = T, row.names = 1)
    as.matrix(x)
}

#' Compute TAP transport propensity score
#'
#' The formula we use is generalizable to different k-mers, but we
#' apply it to 9-mers. It takes a weighted sum of the scores from
#' Positions 1, 2, 3, and 9. For ease of plotting, we orient the
#' score to link higher scores with better transport.
tapscore_compute <- function(peptides, tapmat) {
    pmat <- do.call("rbind", strsplit(peptides, ""))
    colnames(pmat) <- paste0("P", 1:9)

    poss <- c("P1", "P2", "P3", "P9")

    sc <- matrix(NA, nrow = nrow(pmat), ncol = length(poss))
    rownames(sc) <- peptides
    colnames(sc) <- poss

    for (pos in poss) {
        aa <- pmat[,pos]
        sc[,pos] <- tapmat[aa, pos]
    }

    sc[,c("P1", "P2", "P3")] <- 0.2 * sc[,c("P1", "P2", "P3")]

    (-1.0) * rowMeans(sc)
}

tap_f  <- file.path(prereq_dir, 'tap-smmscan-mat.tsv')
tapmat <- tap_read(tap_f)

pepeffect_dir <- file.path(res_dir, 'pepeffect')
f             <- file.path(pepeffect_dir,
                           'dataset.deephlaffy.tsv.gz')
x             <- read.table(f, sep = '\t', header = T,
                            stringsAsFactors = F)

# We next compute the contribution of binding strength and peptide
# effect to the binding scores.
contrib_cols  <- c("hla", "peptide",
                   "OrientedBindingStrengthBind",
                   "OrientedPepEffect")
xcontrib  <- x[x$pred_binder > 0.5,contrib_cols]

# Contributions converted to proportions.
xcontribp <- t(apply(xcontrib[,-c(1:2)], 1, function(x) {
    absx <- abs(x)
    sumx <- sum(absx)
    round(absx/sumx, 3)
}))

# Plot the median contribution of the two components to
# the binder score.
contrib_file <- file.path(res_dir, 'contribution.tiff')
tiff_open(contrib_file, width = 2000)
par(oma = c(8, 3, 3, 3), mar = c(3, 8, 3, 3), family = 'symbol')
bp <- barplot(setNames(apply(xcontribp, 2, median),
                 c("Binding Strength", "Peptide Processing")),
        col = c("#E16173", "#729FCF"),
        ylab = "Contribution to\npresentation probability",
        xaxt = 'n', ylim = c(0, 1))
text(bp - 0.20, c(-0.15, -0.15), srt = 45, xpd = NA, cex = 0.7,
     labels = c("Binding\nStrength", "Peptide\nProcessing"))
tiff_close(contrib_file)

# We compute the TAP score for all the peptides.
pep_df <- setNames(unique(x[,c("peptide", "OrientedPepEffect")]),
                   c("peptide", "effect"))
pep_df$good <- as.numeric(pep_df$effect > 0)
pep_df$`TAP score` <- tapscore_compute(pep_df$peptide, tapmat)

# Plot the distribution of the peptide score.
scoredist_file <- file.path(pepeffect_dir, 'scoredist.tiff')
{
tiff_open(scoredist_file, height = 2000)
par(family = 'symbol')
hist(pep_df$effect, xlab = "Peptide processing score", main = NA,
     ylim = c(0, 150000))
abline(v = 0, lty = 5, lwd = 2, col = 'red')
text(x = 5, y = 100000,
     labels = sprintf("%d\npeptides", sum(pep_df$good)))
text(x = -7.5, y = 100000,
     labels = sprintf("%d\npeptides", sum(pep_df$good == 0)))
tiff_close(scoredist_file)
}

# We plot the distribution of the TAP scores for positive
# and non-positive peptide processing scores.
tap_file <- file.path(pepeffect_dir, 'tap-score.tiff')
{
tiff_open(tap_file)
par(oma = c(8, 3, 3, 3), mar = c(3, 8, 3, 3), family = 'symbol')
peptap_sp <- split(pep_df$`TAP score`, pep_df$effect > 0)
tmp_ppp <- wilcox.test(peptap_sp[[1]], peptap_sp[[2]])

# 2.2e-16 is reported by the summary when printing
# the results.
stopifnot(tmp_ppp$p.value < 2.2e-16)
boxplot(split(pep_df$`TAP score`, pep_df$effect > 0),
        pch = ".", col = c("#729FCF", "#E16173"),
        ylim = c(-1, 1), ylab = "TAP score", xaxt = 'n')
text(c(1, 2), y = -1.33, xpd = NA,
     labels = c("Less\nPresentable", "More\npresentable"))
legend("bottomright", legend = c("p-value < 2.2e-16"), bty = 'n',
       cex = 0.75)
tiff_close(tap_file)
}

#' Compute position-wise residue frequency across a set of peptides
freq_get <- function(peptides, K) {
    peptides <- unique(peptides)
    stopifnot(all(nchar(peptides) == K))
    psplit <- strsplit(peptides, "")
    aminos     <- "ACDEFGHIKLMNPQRSTVWY"
    aminos_vec <- strsplit(aminos, "")[[1]]

    pmat <- do.call("cbind", lapply(seq_len(K), function(ind) {
        s <- sapply(psplit, `[[`, ind)
        stopifnot(all(s %in% aminos_vec))
        tbl <- setNames(sapply(aminos_vec, function(aa) {
            sum(s %in% aa)
        }), aminos_vec)
    }))
    colnames(pmat) <- paste0("P", seq_len(K))
    pmat
}

# Frequency in the full peptide set.
nmat1 <- freq_get(pep_df$peptide, 9)
n1    <- colSums(nmat1)[1]
pmat1 <- nmat1 / n1

# Frequency in the postive contribution peptide set.
nmat2 <- freq_get(pep_df$peptide[pep_df$effect > 0], 9)
n2    <- colSums(nmat2)[1]
pmat2 <- nmat2 / n2

# Ratio of residue occurence in the good set compared to
# the full peptide set.
rmat      <- round(pmat2 / pmat1, 3)
log2_rmat <- log2(rmat)

# Plot the position-wise residue preference for enhanced
# peptide processing.
res_f <- file.path(pepeffect_dir, 'peptide-processing-propensity.tiff')
br <- seq(-2, 2, 0.5)
rdylbu_r <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
l2mat <- log2_rmat
l2mat[abs(log2_rmat) < log2(1.5)] <- NA
tiff_open(res_f, width = 3000, height = 2500)
pheatmap(t(l2mat), cluster_rows = F, cluster_cols = F,
         breaks = br, col = rdylbu_r(length(br)),
         border_col = "black", display_numbers = F,
         na_col = "white", family = 'symbol', angle_col = 0,
         cellheight = 20, cellwidth = 20,
         legend_breaks = c(-2, -1, 0, 1, 2),
         legend_labels = c("0.25", "0.5", "1", "2", "4"))
tiff_close(res_f)
 

# Compute the position-wise importance using Boruta.
pep_bdf <- local({
    m <- do.call("rbind", strsplit(pep_df$peptide, ""))
    colnames(m) <- paste0("P", seq_len(ncol(m)))
    cbind.data.frame(m, good = pep_df$good)
})

set.seed(1)
print(system.time({
    boruta_res <- Boruta(good ~ ., pep_bdf)
}))

ccode <- c("lightgoldenrod1", "cornflowerblue", "indianred1", "grey80")

pepeffects_posimp_file <- file.path(pepeffect_dir,
                                    'pepeffect-posimp.tiff')
{
tiff_open(pepeffects_posimp_file,
    width = 2000, height = 2400, res = 400)
par(oma = c(6.5, 2, 2, 2), family = 'symbol')
plot(boruta_res, ccode, las = 2, xlab = NA, ylab = "Importance", xaxt = 'n', pch = "-")
nitems <- ncol(boruta_res$ImpHistory)
nms    <- names(sort(apply(boruta_res$ImpHistory, 2, median)))
#axis(side = 1, seq_len(nitems), labels = NA)
offset <- -c(-0.8, -0.8, -0.8, rep(0.1, 9))
text(seq_len(nitems) - offset, c(-32, -32, -32, rep(-18, 9)), labels = nms, xpd = NA, srt = 45, cex = 0.6)
tiff_close(pepeffects_posimp_file)
}

#------------------------------------------------#
#              Embeddings analysis               #
#------------------------------------------------#

emb_dir <- file.path(res_dir, 'embeddings')

emb_fs <- local({
    tmp <- list.files(emb_dir, full.names = T, recursive = T,
                      pattern = '[-]embeddings[.]tsv$')
    setNames(tmp, basename(dirname(tmp)))
})

# We read in the representative VHSE AAIndex1 data for
# comparison later
vhse_f   <- file.path(prereq_dir, 'vhse-repr-index.tsv.xz')
vhse     <- read.table(vhse_f, sep = '\t', header = T)
vhse_mat <- local({
    tmp  <- as.matrix(vhse[,-c(1,2)])
    rownames(tmp) <- vhse$aaindex
    tmp
})

emb_read <- function(f, vhse_mat = NULL) {

    # Read the matrix, without the unknown amino acid "X".
    embm <- local({
        tmp <- read.table(f, sep = '\t', header = T,
                          row.names = 1)
        tmp[rownames(tmp) != "X", ]
    })

    # Calculate distance matrix.
    embm_dist <- apply(embm, 1, function(x) {
        apply(embm, 1, function(y) {
            round(sqrt(mean((x - y) ^ 2)), 3)
        })
    })

    # Hierarchical clustering
    set.seed(1)
    embm_hc <- hclust(as.dist(embm_dist))

    res <- list(
        mat  = embm,
        dist = embm_dist,
        hc   = embm_hc
    )

    if (!is.null(vhse_mat)) {
        vhse_mat <- vhse_mat[,rownames(embm)]
        cr_datf <- local({
            cr <- cor(embm, t(vhse_mat), method = 'spearman')
            tmp <- which(abs(cr) > 0.5, arr.ind = T)
            datf <- data.frame(
                feature   = colnames(embm)[tmp[,1]],
                aaindex   = rownames(vhse_mat)[tmp[,2]],
                SRCC      = round(cr[tmp], 3),
                sign      = ifelse(cr[tmp] >= 0, 1, -1)

            )
            datf$desc <- vhse$description[match(datf$aaindex, vhse$aaindex)]
            datf
        })
        res[["vhse"]] <- cr_datf
    }

    res
}

emb_dat <- lapply(emb_fs, emb_read, vhse_mat = vhse_mat)

{
# We plot the clustering of the amino acids across
# the different embedding systems.
emb_file <- file.path(emb_dir, 'embeddings.tiff')
tiff_open(emb_file, width = 2850, height = 1800, res = 300)
par(mfrow = c(2, 3), family = 'symbol')
nm_ord <- c("peptide", "hla", "pepeffects",
            "peptide-context", "hla-context")
for (nm in nm_ord) {
    nmstr <- if (nm == "hla") {
        "Pair-potential\nHLA embedding"
    } else if (nm == "peptide") {
        "Pair-potential\nPeptide embedding"
    } else if (nm == "hla-context") {
        "Env. context\nHLA embedding"
    } else if (nm == "peptide-context") {
        "Env. context\nPeptide embedding"
    } else if (nm == "pepeffects") {
        "Processing score\nPeptide embedding"
    } else {
        stop("Unknown entity")
    }
    plot(emb_dat[[nm]]$hc, main = nmstr,
         sub = NA, xlab = NA)
}
tiff_close(emb_file)
}

vhse_dat <- lapply(emb_dat, `[[`, "vhse")

vhse_cats <- sort(Reduce(union, lapply(vhse_dat, `[[`, "desc")))

# We identify which VHSE AAIndex1 categories are correlated
# with in each embedding system and then plot them.
a <- lapply(vhse_dat, function(x) {
    vhse_cats %in% x$desc
})

aa <- local({
    tmp0 <- do.call("cbind", a) * 1
    rownames(tmp0) <- vhse_cats

    tmp <- tmp0

    colnames(tmp) <- sapply(colnames(tmp0), function(nm) {
        nmstr <- if (nm == "hla") {
            "INT: HLA"
        } else if (nm == "peptide") {
            "INT: Peptide"
        } else if (nm == "hla-context") {
            "ENV: HLA"
        } else if (nm == "peptide-context") {
            "ENV: Peptide"
        } else if (nm == "pepeffects") {
            "PEP"
        } else {
            stop("Unknown entity")
        }
    })
    tmp
})


{
vhse_hfile <- file.path(emb_dir, 'embeddings-to-vhse.tiff')
tiff_open(vhse_hfile, width = 3600, height = 4800, res = 300)
pheatmap(
    aa[,c(5, 2, 4, 1, 3)],
    cluster_rows = F, cluster_cols = F,
    display_numbers = F,
    col = c("white", "turquoise2"),
    border_col = 'black', angle = 45,
    fontfamily = 'symbol',
    cellwidth = 25, cellheight = 15,
    legend = F
)
tiff_close(vhse_hfile)
}

for (nm in names(vhse_dat)) vhse_dat[[nm]]$emb <- nm

vhse_emb_datf <- do.call("rbind.data.frame", vhse_dat)
row.names(vhse_emb_datf) <- NULL

write.table(vhse_emb_datf, file.path(emb_dir, 'embeddings-to-vhse.tsv'),
            sep = '\t', row.names = F, quote = F)

sp_odir <- file.path(emb_dir, 'embeddings-by-system')
if (!dir.exists(sp_odir)) dir.create(sp_odir, recursive = T)

emb_sp <- split(vhse_emb_datf, vhse_emb_datf$emb)

for (nm in names(emb_sp)) {
    sp <- emb_sp[[nm]]
    sp_edges <- local({
        tmp <- sp[,c("feature", "aaindex", "SRCC", "sign")]
        setNames(tmp, c("Source", "Target", "SRCC", "Sign"))
    })
    sp_edges$color <- ifelse(sp_edges$SRCC > 0,
                             "#FB6A4A", "#74A9CF")


    write.table(sp_edges, sep = '\t', row.names = F, quote = F,
                file.path(sp_odir, paste0(nm, '-edges.tsv')))

}


nodes <- local({
    feats <- sort(Reduce(union, lapply(emb_dat, function(x) colnames(x$mat))))
    inds  <- vhse$aaindex
    descs <- vhse$description
    data.frame(
        name = c(feats, inds),
        type = c(rep("Feature", length(feats)),
                 rep("AAIndex1", length(inds))),
        disp = c(feats, descs),
        color = c(rep("#84C9EC", length(feats)),
                  rep("#EFD277", length(inds))),
        border_width = c(rep(5, length(feats)),
                         rep(2, length(inds))),
        border_color = "black",
        height = c(rep(75, length(feats)),
                   rep(250, length(inds))),
        width  = c(rep(75, length(feats)),
                   rep(250, length(inds))),
        shape  = c(rep("Ellipse", length(feats)),
                   rep("Rounded Rectangle", length(inds)))
    )
})
write.table(nodes, sep = '\t', row.names = F, quote = F,
            file.path(sp_odir, 'node-table.tsv'))


out_dat <- lapply(emb_dat, function(dat) {
    ed <- dat$vhse
    ed <- ed[, colnames(ed) != "sign"]
    colnames(ed) <- sapply(colnames(ed), function(x) {
        switch(x,
               "feature" = "Feature",
               "aaindex" = "AAIndex1",
               "SRCC"    = "SRCC",
               "desc"    = "Description",
               x
        )
    })
    ed
})
names(out_dat) <- sapply(names(emb_dat), function(nm) {
    switch(nm,
        "hla-context"     = "ENV-HLA",
        "peptide-context" = "ENV-PEP",
        "hla"             = "INT-HLA",
        "peptide"         = "INT-PEP",
        "pepeffects"      = "PEP",
        nm
    )
})
out_dat <- out_dat[sort(names(out_dat))]
    

xlsx_f <- file.path(emb_dir, "embeddings-vhse.xlsx")
write_xlsx(out_dat, xlsx_f)
