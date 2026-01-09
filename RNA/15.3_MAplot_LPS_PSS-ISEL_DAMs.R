library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DESeq2)
library(data.table)
library(qvalue)
library(annotables)

library(qqplotr)
library(tidyverse)
library(qqman)


#######################
#### scatter plot #####

rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/LPS_ATAC/v1/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

#platePrefix <- "ALL.0.15.13." 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq" 


#2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_summary/2_th20_plotData.comb.txt.gz
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz"
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_summary/2_th20_plotData.comb.txt.gz"
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz"
fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"

dams <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')


# motif activity matrix: 
#2_motif.outs/correct_excludeX/Cluster_res0.07/1.3_YtX_ave.th20.clean.rds
fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/option_nFeature15K_cluster_res0.12/1.3_YtX_ave.th20.clean.rds"
motifmat <- readRDS(fname)

# make metadata
meta <- tibble(sample = colnames(motifmat)) %>%
  separate(sample, into = c("cluster", "treat", "dbgap.ID"), sep = "_", remove = FALSE) %>% as.data.frame()

rownames(meta) <- meta$sample

# assign atac metadata: 
meta$celltype <- "NA"
meta$celltype[meta$Cluster == "C0"] <- "A0 T CD4+"
meta$celltype[meta$Cluster == "C1"] <- "A1 T CD8+"
meta$celltype[meta$Cluster == "C2"] <- "A2 NK"
meta$celltype[meta$Cluster == "C3"] <- "A3 T CD4+"
meta$celltype[meta$Cluster == "C4"] <- "A4 Monocyte"
meta$celltype[meta$Cluster == "C5"] <- "A5 T CD4+"
meta$celltype[meta$Cluster == "C6"] <- "A6 B"
meta$celltype[meta$Cluster == "C7"] <- "A7 T CD4+"
meta$celltype[meta$Cluster == "C8"] <- "A8 T CD4+"
meta$celltype[meta$Cluster == "C9"] <- "A9 T CD4+"
meta$celltype[meta$Cluster == "C10"] <- "A10 DC"

# load DAMs for variables
fname <-  "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
vardams <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
vardams <- vardams %>% filter(psycho_variable %in% c("PSS_all_mean", "ISEL_Mean")) %>% filter(padj_t < 0.1)
vardams$identifier <- vardams$gene 
vardams$gene_cluster <- paste0(vardams$identifier, "_", vardams$Cluster)

# assign atac metadata: 
vardams$celltype <- "NA"
vardams$celltype[vardams$Cluster == "C0"] <- "A0 T CD4+"
vardams$celltype[vardams$Cluster == "C1"] <- "A1 T CD8+"
vardams$celltype[vardams$Cluster == "C2"] <- "A2 NK"
vardams$celltype[vardams$Cluster == "C3"] <- "A3 T CD4+"
vardams$celltype[vardams$Cluster == "C4"] <- "A4 Monocyte"
vardams$celltype[vardams$Cluster == "C5"] <- "A5 T CD4+"
vardams$celltype[vardams$Cluster == "C6"] <- "A6 B"
vardams$celltype[vardams$Cluster == "C7"] <- "A7 T CD4+"
vardams$celltype[vardams$Cluster == "C8"] <- "A8 T CD4+"
vardams$celltype[vardams$Cluster == "C9"] <- "A9 T CD4+"
vardams$celltype[vardams$Cluster == "C10"] <- "A10 DC"
trl_deseqres_sig <- vardams %>% filter(psycho_variable == "PSS_all_mean")
iselctrl_deseqres_sig <- vardams %>% filter(psycho_variable == "ISEL_Mean")


##############################################33
# prep for scatter plot
treat="LPS"
control="CTRL"


clusters <- unique(meta$cluster)
all_cluster_means <- lapply(clusters, function(clust) {
  
  # Control group
  samples_c <- meta %>%
    filter(treat == control, cluster == clust) %>%
    pull(sample)
  
  mat_c <- motifmat[, colnames(motifmat) %in% samples_c, drop = FALSE]
  ave_c <- rowMeans(mat_c, na.rm = TRUE)

  # Treatment group
  samples_t <- meta %>%
    filter(treat == treat, cluster == clust) %>%
    pull(sample)
  
  mat_t <- motifmat[, colnames(motifmat) %in% samples_t, drop = FALSE]
  ave_t <- rowMeans(mat_t, na.rm = TRUE)

  # Combine into data.frame
  df <- data.frame(
    identifier = rownames(motifmat),
    cluster = clust,
    motif_activity_CTRL = ave_c,
    motif_activity_LPS = ave_t  )

  return(df)
})

# Combine all clusters
cpmb <- do.call(rbind, all_cluster_means)
cpmb$gene_cluster <- paste0(cpmb$identifier, "_", cpmb$cluster)

# Annotate with significance
cpmb$bothsig <- ifelse(
  cpmb$gene_cluster %in% PSSctrl_deseqres_sig$gene_cluster & cpmb$gene_cluster %in% iselctrl_deseqres_sig$gene_cluster, "bothsig",
  ifelse(cpmb$gene_cluster %in% PSSctrl_deseqres_sig$gene_cluster, "PSS_sig",
         ifelse(cpmb$gene_cluster %in% iselctrl_deseqres_sig$gene_cluster, "ISEL_sig", "not_var_sig"))
)

cpmb$bothsig <- factor(cpmb$bothsig, levels = c("not_var_sig", "ISEL_sig", "PSS_sig", "bothsig"))


# assign atac metadata: 
cpmb$celltype <- "NA"
cpmb$celltype[cpmb$Cluster == "C0"] <- "A0 T CD4+"
cpmb$celltype[cpmb$Cluster == "C1"] <- "A1 T CD8+"
cpmb$celltype[cpmb$Cluster == "C2"] <- "A2 NK"
cpmb$celltype[cpmb$Cluster == "C3"] <- "A3 T CD4+"
cpmb$celltype[cpmb$Cluster == "C4"] <- "A4 Monocyte"
cpmb$celltype[cpmb$Cluster == "C5"] <- "A5 T CD4+"
cpmb$celltype[cpmb$Cluster == "C6"] <- "A6 B"
cpmb$celltype[cpmb$Cluster == "C7"] <- "A7 T CD4+"
cpmb$celltype[cpmb$Cluster == "C8"] <- "A8 T CD4+"
cpmb$celltype[cpmb$Cluster == "C9"] <- "A9 T CD4+"
cpmb$celltype[cpmb$Cluster == "C10"] <- "A10 DC"

cpmb_filt <- cpmb %>% filter(!cluster %in% c("C8", "C9"))
cpmb_filt$bothsig <- factor(cpmb_filt$bothsig, levels = c("not_var_sig", "ISEL_sig", "PSS_sig", "bothsig"))
cpmb_filt$bothsig <- factor(cpmb_filt$bothsig, levels = c("bothsig", "ISEL_sig", "PSS_sig", "not_var_sig"))

########33

p <- ggplot(cpmb_filt, aes(motif_activity_CTRL, motif_activity_LPS, col=bothsig)) + 
    geom_point()+
    facet_wrap(.~celltype,ncol=3)+
    geom_abline(slope=1)+
    scale_color_manual(values=c("purple", "red", "blue", "grey60") |> `names<-`(c("bothsig", "PSS_sig","ISEL_sig","not_var_sig"))) +
#png(width = 6, height = 6, file=paste0(outFolder,"figures/","pss_isel_comb-",contrast,"_logcpm.",run,"_scatter.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
#par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
png("01.1_ATAC_LPS_CTRL_motif_activity_PSS_ISEL_sig_colored.png", width = 1500, height = 1200, res = 240)
print(p)
dev.off()


p <- ggplot() + 
  # First, plot gray non-significant points
  geom_point(data = subset(cpmb_filt, bothsig == "not_var_sig"),
             aes(motif_activity_CTRL, motif_activity_LPS, color = "Not significant"),
             alpha = 0.5) +

  # Then, plot the ISEL-only significant points
  geom_point(data = subset(cpmb_filt, bothsig == "ISEL_sig"),
             aes(motif_activity_CTRL, motif_activity_LPS, color = "ISEL only")) +

  # Then, plot the PSS-only significant points
  geom_point(data = subset(cpmb_filt, bothsig == "PSS_sig"),
             aes(motif_activity_CTRL, motif_activity_LPS, color = "PSS only")) +

  # Finally, plot the double significant points
  geom_point(data = subset(cpmb_filt, bothsig == "bothsig"),
             aes(motif_activity_CTRL, motif_activity_LPS, color = "Both")) +

  facet_wrap(. ~ celltype, ncol = 3) +
  geom_abline(slope = 1) +
  scale_color_manual(
    name = "Significance",
    values = c("Not significant" = "grey60", "ISEL only" = "blue", "PSS only" = "red", "Both" = "purple")
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Save the PNG
png("01.1_ATAC_LPS_CTRL_motif_activity_PSS_ISEL_sig_colored.png", width = 1500, height = 1200, res = 240)
print(p)
dev.off()


########################################################
####################### MA plot ########################

#########################################################
########### response LPS DAMs

# treatment dams
#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_summary/2_response_DAMs.comb.txt.gz"
fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_treat/2_th20_plotData.comb.txt.gz"

treats <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

# assign atac metadata: 
treats$celltype <- "NA"
treats$celltype[treats$Cluster == "C0"] <- "A0 T CD4+"
treats$celltype[treats$Cluster == "C1"] <- "A1 T CD8+"
treats$celltype[treats$Cluster == "C2"] <- "A2 NK"
treats$celltype[treats$Cluster == "C3"] <- "A3 T CD4+"
treats$celltype[treats$Cluster == "C4"] <- "A4 Monocyte"
treats$celltype[treats$Cluster == "C5"] <- "A5 T CD4+"
treats$celltype[treats$Cluster == "C6"] <- "A6 B"
treats$celltype[treats$Cluster == "C7"] <- "A7 T CD4+"
treats$celltype[treats$Cluster == "C8"] <- "A8 T CD4+"
treats$celltype[treats$Cluster == "C9"] <- "A9 T CD4+"
treats$celltype[treats$Cluster == "C10"] <- "A10 DC"

treats <- treats %>% filter(padj < 0.1)
#subset(treats, grepl("IRF4", identifier, ignore.case = TRUE))
#subset(treats, grepl("CXCL10", identifier, ignore.case = TRUE))


PSSctrl_deseqres_sig <- vardams %>% filter(psycho_variable == "PSS_all_mean")
iselctrl_deseqres_sig <- vardams %>% filter(psycho_variable == "ISEL_Mean")

variable = "ISEL_Mean"
variable_name = "Social Support"

#variable = "PSS_all_mean"
#variable_name = "Psychological Stress"

#var_ctrl <- iselctrl_deseqres_sig
damslps <- treats %>% filter(!Cluster == "C9")
damslps$gene_cluster <- paste0(damslps$gene, "_", damslps$Cluster)

      # cpmb_filt$gene_cluster
      sub_cpmb_filt <- cpmb_filt %>% dplyr::select(gene_cluster, motif_activity_CTRL, motif_activity_LPS)


clusters <- unique(damslps$Cluster)
all_cluster_means <- lapply(clusters, function(c) {

      var_ctrl <- vardams %>% filter(psycho_variable == variable)

      lpsvctrl <- subset(damslps, Cluster==c)

      lpsvctrl$identifier <- lpsvctrl$gene
      ctrl_deseqres_sig <- subset(var_ctrl, padj_t<0.1 & Cluster==c)
      ctrl_deseqres_upsig <- subset(ctrl_deseqres_sig, estimate>0)
      ctrl_deseqres_downsig <- subset(ctrl_deseqres_sig, estimate<0)
      lpsvctrl <- transform(lpsvctrl, varDAM=ifelse(identifier %in% ctrl_deseqres_upsig$identifier,"increase",ifelse(identifier %in% ctrl_deseqres_downsig$identifier,"decrease","not_sig")))
      lpsvctrl$varDAM <- factor(lpsvctrl$varDAM, levels = c("not_sig","increase","decrease"))
      lpsvctrl$gene_cluster <- paste0(lpsvctrl$gene, "_", lpsvctrl$Cluster)
      merged <- merge(lpsvctrl, sub_cpmb_filt, by="gene_cluster")
      # cpmb_filt$gene_cluster


      return(merged)

      })
# Combine all clusters
motifb <- do.call(rbind, all_cluster_means)


# assign atac metadata: 
motifb$celltype <- "NA"
motifb$celltype[motifb$Cluster == "C0"] <- "A0 T CD4+"
motifb$celltype[motifb$Cluster == "C1"] <- "A1 T CD8+"
motifb$celltype[motifb$Cluster == "C2"] <- "A2 NK"
motifb$celltype[motifb$Cluster == "C3"] <- "A3 T CD4+"
motifb$celltype[motifb$Cluster == "C4"] <- "A4 Monocyte"
motifb$celltype[motifb$Cluster == "C5"] <- "A5 T CD4+"
motifb$celltype[motifb$Cluster == "C6"] <- "A6 B"
motifb$celltype[motifb$Cluster == "C7"] <- "A7 T CD4+"
motifb$celltype[motifb$Cluster == "C8"] <- "A8 T CD4+"
motifb$celltype[motifb$Cluster == "C9"] <- "A9 T CD4+"
motifb$celltype[motifb$Cluster == "C10"] <- "A10 DC"

motifb %>% filter(celltype == "A5 T CD4+") %>% filter(!varDAM == "not_sig")
motifb %>% filter(celltype == "A1 T CD8+") %>% filter(!varDAM == "not_sig")

#vardamsA5 <- vardams %>% filter(celltype == "A5 T CD4+") 
#treats %>% filter(celltype == "A5 T CD4+") %>% filter(gene %in% vardamsA5$gene)
#motifb %>% filter(celltype == "A5 T CD4+") %>% filter(gene %in% vardamsA5$gene)

table(motifb$Cluster)
motifb <- motifb %>% filter(!Cluster %in% c("C10", "C3", "C5", "C7"))

#motifb %>% filter(celltype == "A5 T CD4+") %>% filter(!varDAM == "not_sig")

p <- ggplot() + 
  # First: plot the not significant (gray) points
  geom_point(data = subset(motifb, varDAM == "not_sig"),
             #aes(motif_activity_LPS, beta, color = "not_sig")) +
             aes(motif_activity_CTRL, beta, color = "not_sig")) +

  # Then: plot the increased (red) points
  geom_point(data = subset(motifb, varDAM == "increase"),
             #aes(motif_activity_LPS, beta, color = "increase")) +
             aes(motif_activity_CTRL, beta, color = "increase")) +

  # Then: plot the decreased (blue) points
  geom_point(data = subset(motifb, varDAM == "decrease"),
             #aes(motif_activity_LPS, beta, color = "decrease")) +
             aes(motif_activity_CTRL, beta, color = "decrease")) +

  #scale_x_log10() + 
  #xlab("mean of motif activity LPS") + 
  xlab("mean of motif activity CTRL") + 
  ylab("beta") + 
  facet_wrap(.~celltype, ncol = 5) +
  ggtitle(paste0(variable_name)) +
  geom_hline(yintercept = 0, col = "grey40") + 
  geom_hline(yintercept = -2, col = "grey40", linetype = "dashed") + 
  geom_hline(yintercept = 2, col = "grey40", linetype = "dashed") + 
  scale_color_manual(
    name = "varDAM",
    values = c("increase" = "red", "decrease" = "blue", "not_sig" = "grey60")
  ) +
  theme_minimal() +
  theme(legend.position = "right", 
    axis.text.x = element_text(size = 6)) #+
  #scale_x_log10(breaks = scales::log_breaks(n = 5))


# Save the PNG
#png("02.1_ATAC_MAplot_LPS_ISEL_Mean_sig_colored.png", width = 1800, height = 600, res = 240)
#png("02.2_ATAC_MAplot_LPS_PSS_all_mean_sig_colored.png", width = 1800, height = 600, res = 240)
png("02.3_ATAC_MAplot_CTRL_ISEL_Mean_sig_colored.png", width = 1800, height = 600, res = 240)
#png("02.4_ATAC_MAplot_CTRL_PSS_all_mean_sig_colored.png", width = 1800, height = 600, res = 240)
print(p)
dev.off()








