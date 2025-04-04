library(scriptloc)
script_dir <- script_dir_get()
source(file.path(script_dir, 'helpers.R'))
projroot   <- file.path(script_dir, '..')
res_dir    <- file.path(projroot, 'results')
tune_dir   <- file.path(res_dir, '01-model-tuning')
tune_f     <- file.path(tune_dir, 'data.rds')


dat  <- readRDS(tune_f)

# This data.frame should have all the performance statistics
# along with the hyperparameters already formatted in the
# correct way.
perf <- local({
    lst <- lapply(dat, function(x) as.data.frame(x$summ))
    tmp <- do.call("rbind.data.frame", lst)
    rownames(tmp) <- NULL
    tmp$BaseEmbDim <- dim2str(tmp$pp_embdim)
    tmp$Contacts   <- annotate_contacts(tmp$contacts)
    tmp1 <- cbind.data.frame(tmp, pos_conf_parse(tmp$pos_conf))
    tmp2 <- cbind.data.frame(tmp1, fx_conf_parse(tmp$effects_conf))
    #deal_with_nan(tmp2)
})

reg_perf  <- rank_regperf(subset(perf, data_type == "regressand"))
both_perf <- rank_perf(subset(perf, data_type == "both"))

perf_splits <- list(
    Base       = subset(reg_perf , pos_conf    == "None"),
    PosActv    = subset(reg_perf , (EnvEmbDim  == "016") &
                                   (EnvNode    == "016")),
    PosDimNode = subset(reg_perf , (EnvActvIn  == "Sigmoid") &
                                   (EnvActvOut == "Sigmoid")),
    FxActv     = subset(both_perf, (FxEmbDim   == "016") &
                                   (FxNode     == "016")),
    FxDimNode  = subset(both_perf, (FxActvIn   == "GELU") &
                                   (FxActvOut  == "Linear"))
)

perf_cols <- list(
    "Base"       = c("Contacts" , "BaseEmbDim"),
    "PosActv"    = c("EnvActvOut", "EnvActvIn"),
    "PosDimNode" = c("EnvEmbDim", "EnvNode"),
    "FxActv"     = c("FxActvOut" , "FxActvIn"),
    "FxDimNode"  = c("FxEmbDim" , "FxNode")
)

perf_mats <- local({
    tmp <- list()
    for (nm in names(perf_splits)) {
        cls <- perf_cols[[nm]]
        sp  <- perf_splits[[nm]]
        tmp[[nm]][["SRCC"]]      <- srcc_gen(sp, cls)
        tmp[[nm]][["MSE"]]       <- mse_gen(sp, cls)
        if (nm %in% c("FxActv", "FxDimNode")) {
            tmp[[nm]][["BCE"]]   <- bce_gen(sp, cls)
            tmp[[nm]][["AUC"]]   <- auc_gen(sp, cls)
            tmp[[nm]][["PRAUC"]] <- prauc_gen(sp, cls)
        }
    }
    tmp
})

regtune_dir <- file.path(tune_dir, '01-regression')
if (!dir.exists(regtune_dir)) dir.create(regtune_dir, recursive = T)

dualtune_dir <- file.path(tune_dir, '02-dualtask')
if (!dir.exists(dualtune_dir)) dir.create(dualtune_dir, recursive = T)


local({
for (nm in names(perf_mats)) {

    pmats <- perf_mats[[nm]]

    base_odir <- if ("BCE" %in% names(pmats)) dualtune_dir else regtune_dir
    srcc_fpath <- file.path(base_odir, paste0(nm, 'SRCC.tiff'))
    tiff_open(srcc_fpath)
    srcc_hyperp_ph(pmats$SRCC, reg_pal1)
    tiff_close(srcc_fpath)

    mse_fpath <- file.path(base_odir, paste0(nm, 'MSE.tiff'))
    tiff_open(mse_fpath)
    mse_hyperp_ph(pmats$MSE, reg_pal2)
    tiff_close(mse_fpath)

    if ("BCE" %in% names(pmats)) {
        bce_fpath <- file.path(base_odir, paste0(nm, 'BCE.tiff'))
        tiff_open(bce_fpath)
        bce_hyperp_ph(pmats$BCE, bind_pal1)
        tiff_close(bce_fpath)

        # This prevents issues if there is only one value
        auc_fpath <- file.path(base_odir, paste0(nm, 'AUC.tiff'))
        tiff_open(auc_fpath)
        auc_hyperp_ph(pmats$AUC, bind_pal2)
        tiff_close(auc_fpath)


        prauc_fpath <- file.path(base_odir, paste0(nm, 'PRAUC.tiff'))
        tiff_open(prauc_fpath)
        prauc_hyperp_ph(pmats$PRAUC, bind_pal3)
        tiff_close(prauc_fpath)
    }

}
})


# Provenance
prov_lst <- list(
    'Baseline' = subset(both_perf, (Contacts == "All") &
                                   (pos_conf == "None") &
                                   (effects_conf == "None")),

    '(+) Contacts' = subset(both_perf, (Contacts == "Observed") &
                                       (pos_conf == "None") &
                                       (effects_conf == "None")),

    '(+) Env. context' = subset(both_perf,
                                (Contacts == "Observed") &
                                (pos_conf == "8-sigmoid128_sigmoid53") &
                                (effects_conf == "None")),

    '(+) Effects' = both_perf[1,]
)


subsetm <- function(x, select_cols) as.matrix(x[,select_cols, drop = F])

prov <- do.call("rbind", prov_lst)[,all_meas]

prov_dir <- file.path(tune_dir, '03-provenance')
if (!dir.exists(prov_dir)) dir.create(prov_dir, recursive = T)

mse_file <- file.path(prov_dir, '01-MSE.tiff')
tiff_open(mse_file)
mse_hyperp_ph(subsetm(prov, "MSE"), pal = reg_pal2,
              best_ncol = "white")
tiff_close(mse_file)

srcc_file <- file.path(prov_dir, '02-SRCC.tiff')
tiff_open(srcc_file)
srcc_hyperp_ph(subsetm(prov, "SRCC"), pal = reg_pal1,
               best_ncol = "white")
tiff_close(srcc_file)

bce_file <- file.path(prov_dir, '03-BCE.tiff')
tiff_open(bce_file)
bce_hyperp_ph(subsetm(prov, "BCE"), pal = bind_pal1, best_ncol = "white")
tiff_close(bce_file)

auc_file <- file.path(prov_dir, '04-AUC.tiff')
tiff_open(auc_file)
auc_hyperp_ph(subsetm(prov, "AUC"), pal = bind_pal2, best_ncol = "white")
tiff_close(auc_file)

prauc_file <- file.path(prov_dir, '05-PRAUC.tiff')
tiff_open(prauc_file)
prauc_hyperp_ph(subsetm(prov, "PRAUC"), pal = bind_pal3, best_ncol = "white")
tiff_close(prauc_file)

library(Boruta)

set.seed(0)
reg_boruta <- Boruta(MSE ~ BaseEmbDim + Contacts + EnvEmbDim + EnvNode + EnvActvIn + EnvActvOut, reg_perf, pValue = 1e-6, maxRuns = 2000)

set.seed(0)
both_boruta <- Boruta(BCE ~ FxEmbDim + FxNode + FxActvIn, both_perf, pValue = 1e-6, maxRuns = 2000)


ccode <- c("lightgoldenrod1", "cornflowerblue", "indianred1", "grey80")
reg_boruta_fpath <- file.path(regtune_dir,
                              'HyperparameterImportanceMSE.tiff')
tiff_open(reg_boruta_fpath, 2700, 3500, 600)
par(oma = c(4, 3, 3, 3), family = "serif")
plot(reg_boruta, ccode, xlab = NA, las = 2, horizontal = F)
tiff_close(reg_boruta_fpath)


both_boruta_fpath <- file.path(dualtune_dir,
                               'HyperparameterImportanceBCE.tiff')
tiff_open(both_boruta_fpath, 3500, 3500, 600)
par(oma = c(4, 3, 3, 3), family = "serif")
plot(both_boruta, ccode, xlab = NA, ylab = NA, horizontal = T, las = 2)
tiff_close(both_boruta_fpath)
