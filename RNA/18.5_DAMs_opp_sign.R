####
library(tidyverse)
library(Matrix)
## library(DESeq2)
## library(biobroom)
library(data.table)
library(motifmatchr)
library(GenomicRanges)
library(SummarizedExperiment)
##
library(cowplot)
library(RColorBrewer)
library(scales)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(ggrastr)
library(tidyr)
library(dplyr)
###
###
rm(list=ls())


#setwd("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/Example_show/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/ATAC/")


############################
### show examples
##########################


#outdir2 <- "./Examples/v4/"
#if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings = F, recursive = T)

####
### need provide cluster, variable, motif names  


cluster_rna <- c("C0", "C4")
cluster_atac <- c("C0", "C3")
celltypes <- c("CD4+ T cell", "Monocyte")
variable <-  c("PSS_all_mean", "ISEL_Mean")
variable_name <-c("Psychological Stress", "social Support")

cl_atac = "C3"
cl_rna <- "C0"
cl_atac = "C0"

#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz
res_motif <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz", header=T)
#motif2 <- sort(unique(res_motif$motif_name))
motif2 <- sort(unique(res_motif$gene))

## motif2 <- motif2[grepl("IRF|STAT|REL", motif2)]
## motif2 <- motif2[grepl("IRF", motif2)]

# Filter for the two variables of interest
sub <- res_motif %>%
  filter(psycho_variable %in% c("ISEL_Mean", "PSS_all_mean"), padj_t<0.1, Cluster == cl_atac) %>%
  dplyr::select(Cluster, psycho_variable, gene, estimate, padj_t)

# Pivot wider to have zval_diff for both variables side-by-side
wide <- sub %>%
  pivot_wider(names_from = psycho_variable, values_from = c(estimate, padj_t)) %>% as.data.frame() %>%
  filter(!is.na(estimate_ISEL_Mean), !is.na(estimate_PSS_all_mean))

# Identify motifs with opposite sign z-scores
opposite_sign_df <- wide %>%
  filter(sign(estimate_ISEL_Mean) != sign(estimate_PSS_all_mean))

fname = paste0("00_C0_DAMs_opposite_sign.txt")
  write.table(opposite_sign_df, fname, sep="\t", row.names=F, quote=F)


###################################################### monocyte: zero with opposite signs

