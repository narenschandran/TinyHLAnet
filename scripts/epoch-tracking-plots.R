library(scriptloc)
script_dir <- script_dir_get()
source(file.path(script_dir, 'helpers.R'))
projroot   <- file.path(script_dir, '..')
res_dir    <- file.path(projroot, 'results')
tune_dir   <- file.path(res_dir, '01-model-tuning')
tune_f     <- file.path(tune_dir, 'data.rds')

regtune_dir <- file.path(tune_dir, '01-regression')
dualtune_dir <- file.path(tune_dir, '02-dualtask')
prov_dir <- file.path(tune_dir, '03-provenance')

prov_dat_f <- file.path(prov_dir, 'prov-models-data.tsv')

dat <- read.table(prov_dat_f, sep = '\t', header = T)

mdl_ord <- c("Baseline"        , "(+) Contacts",
             "(+) Env. context", "(+) Effects")

desc <- local({
    tmp0 <- dat[,c("Epoch", "Model", "ModelKey", "Seed")]
    tmp1 <- tmp0[order(-as.numeric(tmp0$Epoch)),]
    tmp2 <- tmp1[!duplicated(tmp1$ModelKey),]
    stopifnot(all(tmp2$Model %in% mdl_ord))
    tmp2[match(mdl_ord, tmp2$Model),]
})

desc$Color <- c("purple", "blue", "darkred", "black")
colvec <- with(desc, setNames(Color, Model))

datf <- setNames(dat[,c("Epoch", "loss", "val_loss", "Model")],
                  c("Epoch", "Train", "Val", "Model"))


uround_fn <- function(x) ceiling(y <- (x * 100) + 1)/100
dround_fn <- function(x) floor((x * 100) - 1)/100

drng <- function(x) {
    ymin <- dround_fn(min(unlist(x[,c("Train", "Val")])))
    ymax <- uround_fn(max(unlist(x[,c("Train", "Val")])))
    c(ymin, ymax)
}
erng <- function(x) c(0, max(x$Epoch))

ylim_full <- drng(datf)
xlim_full <- erng(datf)

datf_sp <- split(datf, datf$Model)[mdl_ord]

{
prov_plot1 <- file.path(prov_dir, '06-provenance-1.tiff')
tiff_open(prov_plot1, width = 3500, height = 3500)
par(mfrow = c(2, 2), family = "symbol",
    oma = c(5, 5, 0, 0), mar = c(2, 2, 2, 2))
for (nm in mdl_ord) {
    sp <- datf_sp[[nm]]

    ylim_sp <- drng(sp)
    xlim_sp <- erng(sp)
    with(sp, {
        plot(Epoch, Train, xlim = xlim_sp, ylim = ylim_sp,
             xlab = NA, ylab = NA, pch = 19,
             type = 'o', main = nm, col = colvec[[nm]])
        points(Epoch, Val, type = 'o', lty = 'dotted',
               col = colvec[[nm]])
        legend("topright", bty = 'n',
               lty = c("solid", "dotted"),
               pch = c(19, 1),
               legend = c("Train", "Val"),
               col = colvec[[nm]])
    })
}
mtext("# Epochs", side = 1, xpd = T, outer = T, line = 1)
mtext("Dual-task Loss (MSE + BCE)", side = 2, xpd = T,
      outer = T, line = 1)
tiff_close(prov_plot1)
}


{
lout <- matrix(rep(c(1, 1, 2), 2), ncol = 3, byrow = T)
prov_plot2 <- file.path(prov_dir, '07-provenance-2.tiff')
tiff_open(prov_plot2, width = 3000, height = 3500, res = 750)
nm <- mdl_ord[[1]]
sp <- datf_sp[[nm]]
layout(lout)
par(family = "symbol", mar = c(4, 4, 0, 0))
with(sp, {
    plot(Epoch, Train, xlim = xlim_full, ylim = ylim_full,
         xlab = "# Epochs", ylab = "Dual-task Loss (MSE + BCE)", pch = 19,
         type = 'o', main = NA, col = colvec[[nm]])
    points(Epoch, Val, type = 'o', lty = 'dotted',
           col = colvec[[nm]])
})


for (nm in mdl_ord[-1]) {
    sp <- datf_sp[[nm]]

    ylim_sp <- drng(sp)
    xlim_sp <- erng(sp)
    with(sp, {
        points(Epoch, Train, xlim = xlim_full, ylim = ylim_full,
             xlab = NA, ylab = NA, pch = 19,
             type = 'o', col = colvec[[nm]])
        points(Epoch, Val, type = 'o', lty = 'dotted',
               col = colvec[[nm]])
    })
}
par(family = "symbol", mar = c(0, 1, 0, 0))

plot.new()
legend("topleft", lty = 0, pch = 15,
       col = desc$Color, legend = desc$Model)

legend("left",
        lty = c("solid", "dotted"),
        pch = c(19, 1),
        legend = c("Train", "Val"))
tiff_close(prov_plot2)
}

