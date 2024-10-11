#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)
library(optparse)
library(ggplot2)
library(gplots)
library(ggExtra)
library(cowplot)
library(tidyr)
library(dplyr)
library(yaml)

theme_set(theme_bw(base_size = 12))


option_list <- list(
  make_option(
    c("--bestfit"), type = "character", help = "Joint qc metrics file"
  ),
  make_option(
    c("--model"), type = "character", help = "Joint qc metrics file"
  ),
  make_option(
    c("--droplet_utils_nuclei"), type = "character", help = "DropletUtils nuclei barcodes list"
  ),
  make_option(
    c("--outdir"), type = "character", help = "Output directory for plots"
  ),
  make_option(
    c("--library"), type = "character", help = "Sample name"
  )
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# ## Testing
# opts=list()
# #opts$qc = "results/joint_qc_gmm/10965-VD-2/joint_qc.csv"
# opts$bestfit = "results/joint_qc_gmm/10965-VD-2/gmm_rna_umis__rna_fraction_mitochondrial_bestfit.csv"
# opts$model = "gmm_rna_umis__rna_fraction_mitochondrial"
# opts$droplet_utils_nuclei = "results/droplet_utils/10965-VD-2/barcodes_nuclei.txt"
# opts$outdir = "results/joint_qc_gmm/10965-VD-2"
# opts$library = "10965-VD-2"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Reading files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Reading files", opts$bestfit, "\n")
d <- read.csv(opts$bestfit)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Infer which GMM component is "Pass nuclei"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d$Component = paste0("Component_", d$bestfit + 1)

sumTable = d %>%
  group_by(Component) %>%
  summarize(
    median_rna_umis = median(rna_umis, na.rm = TRUE),
    one_minus_median_rna_fraction_mitochondrial = 1 - median(rna_fraction_mitochondrial, na.rm = TRUE),
    #median_atac_hqaa = median(atac_hqaa, na.rm = TRUE),
    median_atac_tss_enrichment = median(atac_tss_enrichment, na.rm = TRUE)
  )
    #n = n()

sumTable = as.data.frame(sumTable)
winner_list = list()
for (i in 2:ncol(sumTable)) {
  colname = names(sumTable)[i]
  winner_list[[colname]] = sumTable[order(sumTable[, colname], decreasing = TRUE), "Component"][1]
}

winner = unique(unlist(winner_list))
if (!length(winner) == 1) {
  cat("CAUTION: the gaussian components cannot clearly be assigned as pass QC and fail QC components. Please review the QC plots manually.\n")
  d$Component_plotting = factor(d$Component, levels = paste0("Component_", sort(unique(d$bestfit + 1))))
} else {
  d$Component_plotting = factor(d$Component == winner, levels = c(TRUE, FALSE), labels = c("Pass QC barcodes", "Fail QC barcodes"))
}

# Add counts to factor labels
component_levels = levels(d$Component_plotting)
component_labels = NULL
for (val in component_levels) {
  n = sum(d$Component_plotting == val)
  label = paste0(val, " (n = ", n, ")")
  component_labels <- c(component_labels, label)
}

d$Component_plotting_with_n = factor(d$Component_plotting, levels = component_levels, labels = component_labels)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cross modality visualization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Visualizing qc metrics across modalities\n")
metrics = list(
  rna_reads = list(var = "rna_umis", label = "UMIs in GEX"),
  rna_mt = list(var = "rna_fraction_mitochondrial", label = "Fraction MT in GEX"),
  rna_nuc_frac = list(var = "rna_exon_to_full_gene_body_ratio", label = "RNA exonic fraction"),
  atac_reads = list(var = "atac_hqaa", label = "HQAA in ATAC"),
  atac_mt = list(var = "atac_fraction_mitochondrial", label = "Fraction MT in ATAC"),
  atac_tss = list(var = "atac_tss_enrichment", label = "TSS enrichment in ATAC")
  #atac_autosome = list(var = "atac_max_fraction_reads_from_single_autosome", label = "Max fraction HQAA from single autosome")
)

### Get rid of zero values for plotting
for (i in 1:length(metrics)) {
  m = metrics[[i]]$var
  if (sum(d[, m] == 0, na.rm = TRUE) > 0) {
    plus_val = 0.5 * min(d[!d[, m] == 0, m], na.rm = TRUE)
    cat(m, "=", m, "+", plus_val, "\n")
    d[, m] <- d[, m] + plus_val
  }
}

make2DDensityPlot = function(d, x, y) {
  p = ggplot(d, aes_string(x = x, y = y)) +
    geom_bin2d(bins = 100) +
    geom_point(alpha = 0) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    scale_fill_distiller(palette = 4, direction = 1)
  p2 = ggMarginal(p, type = "density")
  return(p2)
}

make2DDensityPlot_withFit = function(d, x, y) {
  p = ggplot(d, aes_string(x = x, y = y)) +
    #geom_bin2d(bins = 100) +
    geom_point(alpha = 0.1, shape = 16, size = 2, aes(colour = Component_plotting_with_n)) +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    scale_fill_distiller(palette = 4, direction = 1) +
    theme(legend.title = element_blank())
  p2 = ggMarginal(p, type = "density")
  return(p2)
}


addTitle = function(plots, titletext) {
  title <- ggdraw() + draw_label(titletext, fontface = 'bold', x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))
  p = plot_grid(title, plots, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control vertical title margins
  return(p)
}

makeDensityGrid = function(d, plotfile) {
  qcplots = list()
  qcplots[[1]] = make2DDensityPlot(d, x = metrics[["rna_reads"]]$var, y = metrics[["atac_reads"]]$var)
  qcplots[[2]] = make2DDensityPlot_withFit(d, x = metrics[["rna_reads"]]$var, y = metrics[["atac_reads"]]$var)
  qcplots[[3]] = make2DDensityPlot(d, x = metrics[["rna_reads"]]$var, y = metrics[["rna_mt"]]$var)
  qcplots[[4]] = make2DDensityPlot_withFit(d, x = metrics[["rna_reads"]]$var, y = metrics[["rna_mt"]]$var)
  qcplots[[5]] = make2DDensityPlot(d, x = metrics[["rna_reads"]]$var, y = metrics[["atac_tss"]]$var)
  qcplots[[6]] = make2DDensityPlot_withFit(d, x = metrics[["rna_reads"]]$var, y = metrics[["atac_tss"]]$var)
  qcplots[[7]] = make2DDensityPlot(d, x = metrics[["rna_reads"]]$var, y = metrics[["rna_nuc_frac"]]$var)
  qcplots[[8]] = make2DDensityPlot_withFit(d, x = metrics[["rna_reads"]]$var, y = metrics[["rna_nuc_frac"]]$var)
  pngfile = plotfile
  cat("Writing plot to", pngfile, "\n")
  plotsonly = plot_grid(plotlist = qcplots, ncol = 2, byrow = TRUE)
  plotsPlusTitle = addTitle(plotsonly, paste(opts$library))
  ggsave(pngfile, plotsPlusTitle, width = 15, height = 13, dpi = 300)
}

makeDensityGrid(d, file.path(opts$outdir, paste0(opts$model, "_bestfit_viz.png")))
