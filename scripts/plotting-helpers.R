library(RColorBrewer)
library(pheatmap)

tiff_open <- local({
    tmp <- tiff
    formals(tmp)$width       <- 2600
    formals(tmp)$height      <- 2400
    formals(tmp)$res         <- 450
    formals(tmp)$compression <- "lzw"
    tmp
})


tiff_close <- function(f, crop = T) {
    dev.off()
    if (crop) {
        fmt <- 'mogrify -compress lzw -trim %s'
        cmd <- sprintf(fmt, f)
        system(cmd)
    }
    f
}

method_names <- function(x) {
    mnames <- c(
        "tinyhlanet"  = "TinyHLAnet",
        "mhcflurry"   = "MHCflurry 2",
        "mixmhcpred"  = "MixMHCpred",
        "netmhcpan"   = "netMHCpan 4.1",
        "transphla"   = "TransPHLA",
        "netmhcpan_p" = "netMHCpan 4.1 (par)"
    )
    stopifnot(all(x %in% names(mnames)))
    mnames[x]
}

method_colors <- function(x) {
    mcols <- c(
      "tinyhlanet" = "#E16173",
      "mhcflurry"  = "#729FCF",
      "mixmhcpred" = "black",
      "netmhcpan"  = "#6B5E9B",
      "transphla"  = "#FF972F",
      "netmhcpan_p"  = "#6B5E9B"
    )
    stopifnot(all(x %in% names(mcols)))
    mcols[x]
}
