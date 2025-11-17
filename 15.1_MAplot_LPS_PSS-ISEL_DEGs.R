library(EnhancedVolcano)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DESeq2)
library(data.table)
library(qvalue)
library(annotables)
library(tidyverse)
library(qqman)


#######################
#### scatter plot #####

rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/LPS_RNA/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 
variables <- c("PSS_all_mean", "ISEL_Mean")#, "cytocomp")


#2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_summary/2_th20_plotData.comb.txt.gz
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz"
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_summary/2_th20_plotData.comb.txt.gz"
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz"
#dams <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

# deseq degs
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resallc <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegsc <- resallc %>% filter(padj < 0.1)
resdegsc <- resdegsc %>% filter(var %in% variables)
resdegsc <- resdegsc %>% mutate(clust = gsub("^C", "", cluster))
  resdegsc$gene_cluster = paste0(resdegsc$identifier, "_",resdegsc$clust)


# deseq degs
#fname =  paste(dataDir, "allres_LPS_allvars_noCombat_DESeq.txt", sep = "")
#resallt <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
#resdegst <- resallt %>% filter(padj < 0.1)
#resdegst <- resdegst %>% filter(var %in% variables)
#resdegst <- resdegst %>% mutate(clust = gsub("^C", "", cluster))
#  resdegst$gene_cluster = paste0(resdegst$identifier, "_",resdegst$clust)


# norlamlized gene expression matrix
#normdir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data/"
normdir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/norm_counts/"

#motifmat <- readRDS(fname)


##############################################33
# prep for scatter plot
treat="LPS"
control="CTRL"
clusters <- as.character(0:4)
celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
  #all_cluster_means <- lapply(clusters, function(clust) {


allcounts <- list()

  for(c in 1:length(clusters)){

  clust <- clusters[c]
  celltype <- celltypes[c]

  # load CTRL  
  fname <- paste0(normdir, "C",clust, "_CTRL_rna.normal.rds")
  countsc <- readRDS(fname)
  cat("dim norm CTRL counts matrix ", clust,": ",dim(countsc),"\n")

  # make metadata
  metac <- tibble(sample = colnames(countsc)) %>%
    separate(sample, into = c("cluster", "batch", "dbgap.ID", "treat"), sep = "_", remove = FALSE) %>% as.data.frame()
  rownames(metac) <- metac$sample
  metac$celltype <- celltype

  # Load LPS
  fname <- paste0(normdir, "C",clust, "_LPS_rna.normal.rds")
  countst <- readRDS(fname)
  cat("dim norm LPS counts matrix ", clust,": ",dim(countst),"\n")
  
  # make metadata
  metat <- tibble(sample = colnames(countst)) %>%
    separate(sample, into = c("cluster", "batch", "dbgap.ID", "treat"), sep = "_", remove = FALSE) %>% as.data.frame()
  rownames(metat) <- metat$sample
  metat$celltype <- celltype

  #clusters <- unique(meta$cluster)
  #all_cluster_means <- lapply(clusters, function(clust) {
  
  # Control group
  #samples_c <- meta %>%filter(treat == control, cluster == clust) %>%ull(sample)
  samples_c <- metac$sample
  samples_t <- metat$sample

  mat_c <- countsc[, colnames(countsc) %in% samples_c, drop = FALSE]
  ave_c <- rowMeans(mat_c, na.rm = TRUE)

  mat_t <- countst[, colnames(countst) %in% samples_t, drop = FALSE]
  ave_t <- rowMeans(mat_t, na.rm = TRUE)

# Find common genes
  common_genes <- intersect(names(ave_t), names(ave_c))
  ave_t_sub <- ave_t[common_genes]
  ave_c_sub <- ave_c[common_genes]

  # Combine into one big data frame
  df <- data.frame(
    identifier = common_genes,
    cluster = clust,
    celltype = celltype,
    gene_expression_CTRL = ave_c_sub,
    gene_expression_LPS = ave_t_sub,
    stringsAsFactors = FALSE
  )
  df$gene_cluster = paste0(df$identifier, "_",df$cluster)


PSSctrl_deseqres_sig <- resdegsc %>% filter(var == "PSS_all_mean") %>% filter(clust == clust)
iselctrl_deseqres_sig <- resdegsc %>% filter(var == "ISEL_Mean") %>% filter(clust == clust)

# Annotate with significance
df$bothsig <- ifelse(
  df$gene_cluster %in% PSSctrl_deseqres_sig$gene_cluster & df$gene_cluster %in% iselctrl_deseqres_sig$gene_cluster, "bothsig",
  ifelse(df$gene_cluster %in% PSSctrl_deseqres_sig$gene_cluster, "PSS_sig",
         ifelse(df$gene_cluster %in% iselctrl_deseqres_sig$gene_cluster, "ISEL_sig", "not_var_sig"))
)
  df$bothsig <- factor(df$bothsig, levels = c("not_var_sig", "ISEL_sig", "PSS_sig", "bothsig"))

  #return(df)
  allcounts[[c]] <- df


}

dim(allcounts[[1]])
dim(allcounts[[2]])
dim(allcounts[[3]])
dim(allcounts[[4]])

table(allcounts[[1]]$bothsig)
table(allcounts[[2]]$bothsig)
table(allcounts[[3]]$bothsig)
table(allcounts[[4]]$bothsig)

#)

# Combine all clusters
cpmb <- do.call(rbind, allcounts)
#cpmb$gene_cluster <- paste0(cpmb$identifier, "_", cpmb$cluster)

cpmb_filt <- cpmb

#cpmb_filt$bothsig <- factor(cpmb_filt$bothsig, levels = c("not_var_sig", "ISEL_sig", "PSS_sig", "bothsig"))
cpmb_filt$bothsig <- factor(cpmb_filt$bothsig, levels = c("bothsig", "ISEL_sig", "PSS_sig", "not_var_sig"))


########33

p <- ggplot(cpmb_filt, aes(gene_expression_CTRL, gene_expression_LPS, col=bothsig)) + 
    geom_point()+
    facet_wrap(.~celltype,ncol=3)+
    geom_abline(slope=1)+
    scale_color_manual(values=c("purple", "red", "blue", "grey60") |> `names<-`(c("bothsig", "PSS_sig","ISEL_sig","not_var_sig"))) +
#png(width = 6, height = 6, file=paste0(outFolder,"figures/","pss_isel_comb-",contrast,"_logcpm.",run,"_scatter.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
#par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
png("01.1_RNA_LPS_CTRL_gene_expression_PSS_ISEL_sig_colored.png", width = 1500, height = 1200, res = 240)
print(p)
dev.off()


p <- ggplot() + 
  # First, plot gray non-significant points
  geom_point(data = subset(cpmb_filt, bothsig == "not_var_sig"),
             aes(gene_expression_CTRL, gene_expression_LPS, color = "Not significant"),
             alpha = 0.5) +

  # Then, plot the ISEL-only significant points
  geom_point(data = subset(cpmb_filt, bothsig == "ISEL_sig"),
             aes(gene_expression_CTRL, gene_expression_LPS, color = "ISEL only")) +

  # Then, plot the PSS-only significant points
  geom_point(data = subset(cpmb_filt, bothsig == "PSS_sig"),
             aes(gene_expression_CTRL, gene_expression_LPS, color = "PSS only")) +

  # Finally, plot the double significant points
  geom_point(data = subset(cpmb_filt, bothsig == "bothsig"),
             aes(gene_expression_CTRL, gene_expression_LPS, color = "Both")) +

  facet_wrap(. ~ celltype, ncol = 3) +
  geom_abline(slope = 1) +
  scale_color_manual(
    name = "Significance",
    values = c("Not significant" = "grey60", "ISEL only" = "blue", "PSS only" = "red", "Both" = "purple")
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Save the PNG
png("01.1_RNA_LPS_CTRL_gene_expression_PSS_ISEL_sig_colored.png", width = 1500, height = 1200, res = 240)
print(p)
dev.off()




########################################################
####################### load deseq results for vars ########################

#rm(list=ls())

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

#A
### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

resdegs <- resall %>% filter(padj < 0.1)

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("ISEL_Mean", "PSS_all_mean")

resdegs$celltype <- "NA"
resdegs$celltype[resdegs$cluster == "C0"] <- "R0 T CD4+"
resdegs$celltype[resdegs$cluster == "C1"] <- "R1 T CD8+"
resdegs$celltype[resdegs$cluster == "C2"] <- "R2 NK"
resdegs$celltype[resdegs$cluster == "C3"] <- "R3 Monocyte"
resdegs$celltype[resdegs$cluster == "C4"] <- "R4 B"
#resdegs$celltype[resdegs$cluster == "C5"] <- "R6 DC"

resdegs <- resdegs %>% filter(var %in% vars)



########################################################
####################### MA plot ########################

#PSSctrl_deseqres_sig <- vardams %>% filter(psycho_variable == "PSS_all_mean")
#iselctrl_deseqres_sig <- vardams %>% filter(psycho_variable == "ISEL_Mean")

### load LPS vs CTRL results
treatDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/treatment_noCOMBAT/"


platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 
#variables <- c("PSS_all_mean", "ISEL_Mean")#, "cytocomp")

#variable = "ISEL_Mean"
#variable = "PSS_all_mean"

##############################################33
# prep for scatter plot
treat="LPS"
control="CTRL"
clusters <- as.character(0:4)
celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
  #all_cluster_means <- lapply(clusters, function(clust) {


allres <- list()

  for(c in 1:length(clusters)){

  clust <- clusters[c]
  celltype <- celltypes[c]

  # load CTRLvsTreatment  
  fname <- paste0(treatDir, "deseqres/", platePrefix, "C",clust, ".deseqres.RNA-LPS_vs_RNA-CTRL.treatment_noCOMBAT.txt")
  res <- as.data.frame(fread(fname))
  #load(fname)
  #res <- results(dds)
  #cat("dim norm CTRLvsLPS counts matrix ", clust,": ",dim(res),"\n")
  ressig <- res %>% filter(padj < 0.1) 

  #return(df)
  allres[[c]] <- res

}

resall <- do.call(rbind, allres)

variable = "ISEL_Mean"
variable_name = "Social Support"

#variable = "PSS_all_mean"
#variable_name = "Psychological Stress"

      # cpmb_filt$gene_cluster
      sub_cpmb_filt <- cpmb_filt %>% dplyr::select(gene_cluster, gene_expression_CTRL, gene_expression_LPS)


clusters <- c("C0", "C1", "C2", "C3", "C4")
#all_cluster_means <- lapply(clusters, function(c) {
alltreats <- list()

  for(c in 1:length(clusters)){

  clust <- clusters[c]
  celltype <- celltypes[c]


      var_ctrl <- resdegs %>% filter(var == variable)

      lpsvctrl <- subset(resall, cluster==clust)

      #lpsvctrl$identifier <- lpsvctrl$gene
      ctrl_deseqres_sig <- subset(var_ctrl, padj<0.1 & cluster==clust)
      ctrl_deseqres_upsig <- subset(ctrl_deseqres_sig, logFC>0 )
      ctrl_deseqres_downsig <- subset(ctrl_deseqres_sig, logFC<0 )
      lpsvctrl <- transform(lpsvctrl, varDEG=ifelse(identifier %in% ctrl_deseqres_upsig$identifier,"increase",ifelse(identifier %in% ctrl_deseqres_downsig$identifier,"decrease","not_sig")))
      lpsvctrl$varDEG <- factor(lpsvctrl$varDEG, levels = c("not_sig","increase","decrease"))
      lpsvctrl <- lpsvctrl %>%mutate(Cluster = gsub("^C", "", cluster))
      lpsvctrl$gene_cluster <- paste0(lpsvctrl$identifier, "_", lpsvctrl$Cluster)
      merged <- merge(lpsvctrl, sub_cpmb_filt, by="gene_cluster")
      # cpmb_filt$gene_cluster


      #return(merged)

        alltreats[[c]] <- merged


      }
# Combine all clusters
exprb <- do.call(rbind, alltreats)


exprb$celltype <- "NA"
exprb$celltype[exprb$cluster == "C0"] <- "A0 T CD4+"
exprb$celltype[exprb$cluster == "C1"] <- "A1 T CD8+"
exprb$celltype[exprb$cluster == "C2"] <- "A2 NK"
exprb$celltype[exprb$cluster == "C3"] <- "A3 Monocyte"
exprb$celltype[exprb$cluster == "C4"] <- "A4 B"

exprb <- exprb %>% filter(!cluster %in% c("C4"))


p <- ggplot() + 
  # First: plot the not significant (gray) points
  geom_point(data = subset(exprb, varDEG == "not_sig"),
             aes(gene_expression_LPS, logFC, color = "not_sig")) +
             #aes(gene_expression_CTRL, logFC, color = "not_sig")) +

  # Then: plot the increased (red) points
  geom_point(data = subset(exprb, varDEG == "increase"),
             aes(gene_expression_LPS, logFC, color = "increase")) +
             #aes(gene_expression_CTRL, logFC, color = "increase")) +

  # Then: plot the decreased (blue) points
  geom_point(data = subset(exprb, varDEG == "decrease"),
             aes(gene_expression_LPS, logFC, color = "decrease")) +
             #aes(gene_expression_LPS, logFC, color = "decrease")) +

  #scale_x_log10() + 
  xlab("mean of gene expression LPS") + 
  #xlab("mean of gene expression CTRL") + 
  ylab("logFC") + 
  facet_wrap(.~celltype, ncol = 4) +
  ggtitle(paste0(variable_name)) +
  geom_hline(yintercept = 0, col = "grey40") + 
  geom_hline(yintercept = -2, col = "grey40", linetype = "dashed") + 
  geom_hline(yintercept = 2, col = "grey40", linetype = "dashed") + 
  scale_color_manual(
    name = "varDEG",
    values = c("increase" = "red", "decrease" = "blue", "not_sig" = "grey60")
  ) +
  theme_minimal() +
  theme(legend.position = "right", 
    axis.text.x = element_text(size = 6)) #+
  #scale_x_log10(breaks = scales::log_breaks(n = 5))


# Save the PNG
png("02.1_RNA_MAplot_LPS_ISEL_Mean_sig_colored.png", width = 1800, height = 600, res = 240)
#png("02.2_RNA_MAplot_LPS_PSS_all_mean_sig_colored.png", width = 1800, height = 600, res = 240)
#png("02.3_RNA_MAplot_CTRL_ISEL_Mean_sig_colored.png", width = 1800, height = 600, res = 240)
#png("02.4_RNA_MAplot_CTRL_PSS_all_mean_sig_colored.png", width = 1800, height = 600, res = 240)
print(p)
dev.off()


############################################################################################################################
############################### just for the 2 cell types ISEL #############################################################
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/LPS_RNA/")

variable = "ISEL_Mean"
variable_name = "Social Support"

#variable = "PSS_all_mean"
#variable_name = "Psychological Stress"

#)

#damslps <- resallc
#damslps$gene_cluster <- paste0(damslps$identifier, "_", damslps$Cluster)

#var_ctrl <- iselctrl_deseqres_sig
#damslps <- dams %>% filter(!Cluster == "C9")
#damslps$gene_cluster <- paste0(damslps$gene, "_", damslps$Cluster)

      # cpmb_filt$gene_cluster
      sub_cpmb_filt <- cpmb_filt %>% dplyr::select(gene_cluster, gene_expression_CTRL, gene_expression_LPS)


#clusters <- c("C0", "C1", "C2", "C3", "C4")
clusters <- c("C0", "C3")

#all_cluster_means <- lapply(clusters, function(c) {
alltreats <- list()

  for(c in 1:length(clusters)){

  clust <- clusters[c]
  celltype <- celltypes[c]


      var_ctrl <- resdegs %>% filter(var == variable)

      lpsvctrl <- subset(resall, cluster==clust)

      #lpsvctrl$identifier <- lpsvctrl$gene
      ctrl_deseqres_sig <- subset(var_ctrl, padj<0.1 & cluster==clust)
      ctrl_deseqres_upsig <- subset(ctrl_deseqres_sig, logFC>0 )
      ctrl_deseqres_downsig <- subset(ctrl_deseqres_sig, logFC<0 )
      lpsvctrl <- transform(lpsvctrl, varDEG=ifelse(identifier %in% ctrl_deseqres_upsig$identifier,"increase",ifelse(identifier %in% ctrl_deseqres_downsig$identifier,"decrease","not_sig")))
      lpsvctrl$varDEG <- factor(lpsvctrl$varDEG, levels = c("not_sig","increase","decrease"))
      lpsvctrl <- lpsvctrl %>%mutate(Cluster = gsub("^C", "", cluster))
      lpsvctrl$gene_cluster <- paste0(lpsvctrl$identifier, "_", lpsvctrl$Cluster)
      merged <- merge(lpsvctrl, sub_cpmb_filt, by="gene_cluster")
      # cpmb_filt$gene_cluster


      #return(merged)

        alltreats[[c]] <- merged


      }
# Combine all clusters
exprb <- do.call(rbind, alltreats)

exprb$celltype <- "NA"
exprb$celltype[exprb$cluster == "C0"] <- "A0 T CD4+"
exprb$celltype[exprb$cluster == "C3"] <- "A3 Monocyte"

exprb$dir <- "NA"
exprb$dir[exprb$varDEG == "not_sig"] <- "Non-significant"
exprb$dir[exprb$varDEG == "decrease"] <- "↑ Social Support & DEG ↓"
exprb$dir[exprb$varDEG == "increase"] <- "↑ Social Support & DEG ↑"

#exprb$dir <- "NA"
#exprb$dir[exprb$varDEG == "not_sig"] <- <- "Non-significant"
#exprb$dir[exprb$varDEG == "decrease"] <- "↑ Psychological Stress & DEG ↓"
#exprb$dir[exprb$varDEG == "increase"] <- "↑ Psychological Stress & DEG ↑"



p <- ggplot() + 
  # First: plot the not significant (gray) points
  geom_point(data = subset(exprb, dir == "Non-significant"),
             #aes(gene_expression_LPS, logFC, color = "not_sig")) +
             aes(gene_expression_CTRL, logFC, color = "Non-significant")) +

  # Then: plot the increased (red) points
  geom_point(data = subset(exprb, dir == "↑ Social Support & DEG ↑"),
             #aes(gene_expression_LPS, logFC, color = "↑ Social Support & DEG ↑")) +
             aes(gene_expression_CTRL, logFC, color = "↑ Social Support & DEG ↑")) +

  # Then: plot the decreased (blue) points
  geom_point(data = subset(exprb, dir == "↑ Social Support & DEG ↓"),
             #aes(gene_expression_LPS, logFC, color = "decrease")) +
             aes(gene_expression_LPS, logFC, color = "↑ Social Support & DEG ↓")) +

  #scale_x_log10() + 
  #xlab("mean of gene expression LPS") + 
  xlab("mean of gene expression CTRL") + 
  ylab("logFC") + 
  facet_wrap(.~celltype, ncol = 4) +
  ggtitle(paste0(variable_name)) +
  geom_hline(yintercept = 0, col = "grey40") + 
  geom_hline(yintercept = -2, col = "grey40", linetype = "dashed") + 
  geom_hline(yintercept = 2, col = "grey40", linetype = "dashed") + 
  scale_color_manual(
    name = NULL,
    values = c("↑ Social Support & DEG ↑" = "red", "↑ Social Support & DEG ↓" = "blue", "Non-significant" = "grey60")
  ) +
  theme_minimal() +
  theme(legend.position = "right", 
    axis.text.x = element_text(size = 6)) #+
  #scale_x_log10(breaks = scales::log_breaks(n = 5))


# Save the PNG
#png("03.1_RNA_MAplot_LPS_ISEL_Mean_sig_colored_2cells.png", width = 1800, height = 600, res = 240)
#png("03.2_RNA_MAplot_LPS_PSS_all_mean_sig_colored_2cells.png", width = 1800, height = 600, res = 240)
png("03.3_RNA_MAplot_CTRL_ISEL_Mean_sig_colored_2cells.png", width = 1300, height = 600, res = 240)
#png("03.4_RNA_MAplot_CTRL_PSS_all_mean_sig_colored_2cells.png", width = 1800, height = 600, res = 240)
print(p)
dev.off()



############################################################################################################################
############################### just for the 2 cell types PSS #############################################################
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/LPS_RNA/")

#variable = "ISEL_Mean"
#variable_name = "Social Support"

variable = "PSS_all_mean"
variable_name = "Psychological Stress"

#)

#damslps <- resallc
#damslps$gene_cluster <- paste0(damslps$identifier, "_", damslps$Cluster)

#var_ctrl <- iselctrl_deseqres_sig
#damslps <- dams %>% filter(!Cluster == "C9")
#damslps$gene_cluster <- paste0(damslps$gene, "_", damslps$Cluster)

      # cpmb_filt$gene_cluster
      sub_cpmb_filt <- cpmb_filt %>% dplyr::select(gene_cluster, gene_expression_CTRL, gene_expression_LPS)


#clusters <- c("C0", "C1", "C2", "C3", "C4")
clusters <- c("C0", "C3")

#all_cluster_means <- lapply(clusters, function(c) {
alltreats <- list()

  for(c in 1:length(clusters)){

  clust <- clusters[c]
  celltype <- celltypes[c]


      var_ctrl <- resdegs %>% filter(var == variable)

      lpsvctrl <- subset(resall, cluster==clust)

      #lpsvctrl$identifier <- lpsvctrl$gene
      ctrl_deseqres_sig <- subset(var_ctrl, padj<0.1 & cluster==clust)
      ctrl_deseqres_upsig <- subset(ctrl_deseqres_sig, logFC>0 )
      ctrl_deseqres_downsig <- subset(ctrl_deseqres_sig, logFC<0 )
      lpsvctrl <- transform(lpsvctrl, varDEG=ifelse(identifier %in% ctrl_deseqres_upsig$identifier,"increase",ifelse(identifier %in% ctrl_deseqres_downsig$identifier,"decrease","not_sig")))
      lpsvctrl$varDEG <- factor(lpsvctrl$varDEG, levels = c("not_sig","increase","decrease"))
      lpsvctrl <- lpsvctrl %>%mutate(Cluster = gsub("^C", "", cluster))
      lpsvctrl$gene_cluster <- paste0(lpsvctrl$identifier, "_", lpsvctrl$Cluster)
      merged <- merge(lpsvctrl, sub_cpmb_filt, by="gene_cluster")
      # cpmb_filt$gene_cluster


      #return(merged)

        alltreats[[c]] <- merged


      }
# Combine all clusters
exprb <- do.call(rbind, alltreats)

exprb$celltype <- "NA"
exprb$celltype[exprb$cluster == "C0"] <- "A0 T CD4+"
exprb$celltype[exprb$cluster == "C3"] <- "A3 Monocyte"

#exprb$dir <- "NA"
#exprb$dir[exprb$varDEG == "not_sig"] <- "Non-significant"
#exprb$dir[exprb$varDEG == "decrease"] <- "↑ Social Support & DEG ↓"
#exprb$dir[exprb$varDEG == "increase"] <- "↑ Social Support & DEG ↑"

exprb$dir <- "NA"
exprb$dir[exprb$varDEG == "not_sig"] <- "Non-significant"
exprb$dir[exprb$varDEG == "decrease"] <- "↑ Psychological Stress & DEG ↓"
exprb$dir[exprb$varDEG == "increase"] <- "↑ Psychological Stress & DEG ↑"



p <- ggplot() + 
  # First: plot the not significant (gray) points
  geom_point(data = subset(exprb, dir == "Non-significant"),
             #aes(gene_expression_LPS, logFC, color = "not_sig")) +
             aes(gene_expression_CTRL, logFC, color = "Non-significant")) +

  # Then: plot the increased (red) points
  geom_point(data = subset(exprb, dir == "↑ Psychological Stress & DEG ↑"),
             #aes(gene_expression_LPS, logFC, color = "↑ Psychological Stress & DEG ↑")) +
             aes(gene_expression_CTRL, logFC, color = "↑ Psychological Stress & DEG ↑")) +

  # Then: plot the decreased (blue) points
  geom_point(data = subset(exprb, dir == "↑ Psychological Stress & DEG ↓"),
             #aes(gene_expression_LPS, logFC, color = "decrease")) +
             aes(gene_expression_LPS, logFC, color = "↑ Psychological Stress & DEG ↓")) +

  #scale_x_log10() + 
  #xlab("mean of gene expression LPS") + 
  xlab("mean of gene expression CTRL") + 
  ylab("logFC") + 
  facet_wrap(.~celltype, ncol = 4) +
  ggtitle(paste0(variable_name)) +
  geom_hline(yintercept = 0, col = "grey40") + 
  geom_hline(yintercept = -2, col = "grey40", linetype = "dashed") + 
  geom_hline(yintercept = 2, col = "grey40", linetype = "dashed") + 
  scale_color_manual(
    name = NULL,
    values = c("↑ Psychological Stress & DEG ↑" = "red", "↑ Psychological Stress & DEG ↓" = "blue", "Non-significant" = "grey60")
  ) +
  theme_minimal() +
  theme(legend.position = "right", 
    axis.text.x = element_text(size = 6)) #+
  #scale_x_log10(breaks = scales::log_breaks(n = 5))



# Save the PNG
#png("03.1_RNA_MAplot_LPS_ISEL_Mean_sig_colored_2cells.png", width = 1800, height = 600, res = 240)
#png("03.2_RNA_MAplot_LPS_PSS_all_mean_sig_colored_2cells.png", width = 1300, height = 600, res = 240)
#png("03.3_RNA_MAplot_CTRL_ISEL_Mean_sig_colored_2cells.png", width = 1300, height = 600, res = 240)
png("03.4_RNA_MAplot_CTRL_PSS_all_mean_sig_colored_2cells.png", width = 1300, height = 600, res = 240)
print(p)
dev.off()










###################################33333


#separate clusters and for the colors maybe we can categorize up vs down regulated and do each variable separately 
lapply(c("PSS_all_mean","ISEL_Mean"), function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_deseqres <- fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",control,".",varrun,".txt"))
        lpsvctrl_a <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
  subvars <- ldply(lapply(names(counts_ls),function(c){
      lpsvctrl <- subset(lpsvctrl_a, cluster==c)
      ctrl_deseqres_sig <- subset(ctrl_deseqres, padj<0.1 & cluster==c)
      ctrl_deseqres_upsig <- subset(ctrl_deseqres_sig, logFC>0 )
      ctrl_deseqres_downsig <- subset(ctrl_deseqres_sig, logFC<0 )
      lpsvctrl <- transform(lpsvctrl, varDEG=ifelse(identifier %in% ctrl_deseqres_upsig$identifier,"upsig",ifelse(identifier %in% ctrl_deseqres_downsig$identifier,"downsig","not_var_sig")))
      lpsvctrl$varDEG <- factor(lpsvctrl$varDEG, levels = c("not_var_sig","upsig","downsig"))
      return(lpsvctrl)
      }), data.frame)

p <- ggplot(subvars[order(subvars$varDEG),], aes(baseMean, logFC, col=varDEG, label=identifier)) + 
        geom_point() + scale_x_log10() + 
        xlab("mean of normalized counts") + 
        ylab("log fold change") + 
        facet_wrap(.~cluster,ncol=3)+
        ggtitle(paste0(var))+
        geom_hline(yintercept=0, col="grey40") + 
        geom_hline(yintercept=-2, col="grey40",linetype="dashed") + 
        geom_hline(yintercept=2, col="grey40",linetype="dashed") + 
        scale_color_manual(values=c("red", "blue", "grey60") |> `names<-`(c("upsig","downsig","not_var_sig"))) 
        #geom_point(data=tab, shape=1, size=3, show.legend=FALSE) + 

png(width = 6, height = 6, file=paste0(outFolder,"figures/",var,"_comb-",contrast,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(p)
dev.off()
})











#####################################
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",11,0.2) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]


#ahh this is the folder locations etc
#plotting LPS vs CONTROL
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
run="treatment_withCOMBAT"
outFolder=paste0(baseoutFolder,run,"/")
threshold=0.10
treat="RNA-LPS"
control="RNA-CTRL"
contrastdf <- data.frame(control=c("RNA-CTRL","RNA-LPS"),treatment=c("RNA-LPS","RNA-LPS-DEX"))
contrastdf <- transform(contrastdf, contrast=paste0(treatment,"_vs_",control))
contrast=contrastdf[1,]$contrast

varrun="SES_PCs_sex_age_and_treats_generem"
varoutFolder=paste0(baseoutFolder,varrun,"/")




#separate clusters and for the colors maybe we can categorize up vs down regulated and do each variable separately 
lapply(c("PSS_all_mean","ISEL_Mean"), function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_deseqres <- fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",control,".",varrun,".txt"))
        lpsvctrl_a <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
  subvars <- ldply(lapply(names(counts_ls),function(c){
      lpsvctrl <- subset(lpsvctrl_a, cluster==c)
      ctrl_deseqres_sig <- subset(ctrl_deseqres, padj<0.1 & cluster==c)
      ctrl_deseqres_upsig <- subset(ctrl_deseqres_sig, logFC>0 )
      ctrl_deseqres_downsig <- subset(ctrl_deseqres_sig, logFC<0 )
      lpsvctrl <- transform(lpsvctrl, varDEG=ifelse(identifier %in% ctrl_deseqres_upsig$identifier,"upsig",ifelse(identifier %in% ctrl_deseqres_downsig$identifier,"downsig","not_var_sig")))
      lpsvctrl$varDEG <- factor(lpsvctrl$varDEG, levels = c("not_var_sig","upsig","downsig"))
      return(lpsvctrl)
      }), data.frame)

p <- ggplot(subvars[order(subvars$varDEG),], aes(baseMean, logFC, col=varDEG, label=identifier)) + 
        geom_point() + scale_x_log10() + 
        xlab("mean of normalized counts") + 
        ylab("log fold change") + 
        facet_wrap(.~cluster,ncol=3)+
        ggtitle(paste0(var))+
        geom_hline(yintercept=0, col="grey40") + 
        geom_hline(yintercept=-2, col="grey40",linetype="dashed") + 
        geom_hline(yintercept=2, col="grey40",linetype="dashed") + 
        scale_color_manual(values=c("red", "blue", "grey60") |> `names<-`(c("upsig","downsig","not_var_sig"))) 
        #geom_point(data=tab, shape=1, size=3, show.legend=FALSE) + 

png(width = 6, height = 6, file=paste0(outFolder,"figures/",var,"_comb-",contrast,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(p)
dev.off()
})



# for the left plot

# expression in treatment vs expression in control, and then colored by deg for psychosocial factors. average across samples? and have one dot per gene. log transformed values

opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
load(opfn)
counts_ls$C6 <- NULL
metadata_ls$C6 <- NULL
cluster="C0"
cpmb_ls <- ldply(lapply(names(counts_ls),function(c){
PSSctrl_deseqres_sig <- subset(fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_","PSS_all_mean","-",control,".",varrun,".txt")),padj<0.1 & cluster==c)
iselctrl_deseqres_sig <- subset(fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_","ISEL_Mean","-",control,".",varrun,".txt")),padj<0.1 & cluster==c)

opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",c,".",combatrun,".RData")
load(opfn)
cluster_metadata_sce <- metadata_ls[[c]]
cluster_metadata <- data.frame(cluster_metadata_sce)
cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
cluster_metadata_c <- subset(cluster_metadata, treats==control)
cluster_metadata_t <- subset(cluster_metadata, treats==treat)
cluster_counts_c <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_c))]
cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_t))]
dget <- DGEList(counts=cluster_counts_t)
dget <- calcNormFactors(dget)
cpmt <- log2(cpm(dget)+1)
rownames(cpmt) <- rownames(cluster_counts_t)
cpmt_avg <- as.data.frame(rowMeans(cpmt, na.rm = TRUE))
cpmt_avg$identifier <- rownames(cpmt_avg)
dgec <- DGEList(counts=cluster_counts_c)
dgec <- calcNormFactors(dgec)
cpmc <- log2(cpm(dgec)+1)
rownames(cpmc) <- rownames(cluster_counts_c)
cpmc_avg <- as.data.frame(rowMeans(cpmc, na.rm = TRUE))
cpmc_avg$identifier <- rownames(cpmc_avg)

cpmb <- merge(cpmt_avg,cpmc_avg,by="identifier")
names(cpmb)[c(2,3)] <- c(paste0("logcpm_",treat),paste0("logcpm_",control))
cpmb <- transform(cpmb, cluster=c, bothsig=ifelse(identifier %in% PSSctrl_deseqres_sig$identifier & identifier %in% iselctrl_deseqres_sig$identifier,"bothsig", ifelse(identifier %in% PSSctrl_deseqres_sig$identifier, "PSS_all_mean_sig", ifelse(identifier %in% iselctrl_deseqres_sig$identifier, "ISEL_Mean_sig","not_var_sig"))))
cpmb$bothsig <- factor(cpmb$bothsig, levels = c("not_var_sig","ISEL_Mean_sig","PSS_all_mean_sig","bothsig"))
return(cpmb)
}),data.frame)

p <- ggplot(cpmb_ls[order(cpmb_ls$bothsig),], aes(logcpm_RNA.CTRL, logcpm_RNA.LPS, col=bothsig)) + 
    geom_point()+
    facet_wrap(.~cluster,ncol=3)+
    geom_abline(slope=1)+
    scale_color_manual(values=c("purple", "red", "blue", "grey60") |> `names<-`(c("bothsig", "PSS_all_mean_sig","ISEL_Mean_sig","not_var_sig"))) +
png(width = 6, height = 6, file=paste0(outFolder,"figures/","pss_isel_comb-",contrast,"_logcpm.",run,"_scatter.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(p)
dev.off()



# deseq LPS vs control is in here:
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
run="treatment_withCOMBAT"
outFolder=paste0(baseoutFolder,run,"/")

###########################

cluster_metadata_c <- subset(meta, treat==control)
cluster_metadata_t <- subset(meta, treat==treat)
cluster_counts_c <- motifmat[,which(colnames(motifmat) %in% rownames(cluster_metadata_c))]
cluster_counts_t <- motifmat[,which(colnames(motifmat) %in% rownames(cluster_metadata_t))]

#dget <- DGEList(counts=cluster_counts_t)
#dget <- calcNormFactors(dget)
#cpmt <- log2(cpm(dget)+1)
#rownames(cpmt) <- rownames(cluster_counts_t)

motif_ave_c <- as.data.frame(rowMeans(cluster_counts_c, na.rm = TRUE))
motif_ave_c$identifier <- rownames(motif_ave_c)

#dgec <- DGEList(counts=cluster_counts_c)
#dgec <- calcNormFactors(dgec)
#cpmc <- log2(cpm(dgec)+1)
#rownames(cpmc) <- rownames(cluster_counts_c)

motif_ave_t <- as.data.frame(rowMeans(cluster_counts_t, na.rm = TRUE))
motif_ave_t$identifier <- rownames(motif_ave_t)

cpmb <- merge(motif_ave_t,motif_ave_c,by="identifier")
names(cpmb)[c(2,3)] <- c(paste0("motif_activity_",treat),paste0("motif_activity",control))
cpmb$gene_cluster <- paste(cpmb$identifier, "_", cpmb$cluster)


cpmb <- transform(cpmb, cluster=c, bothsig=ifelse(identifier %in% PSSctrl_deseqres_sig$identifier & identifier %in% iselctrl_deseqres_sig$identifier,"bothsig", ifelse(identifier %in% PSSctrl_deseqres_sig$identifier, "PSS_sig", ifelse(identifier %in% iselctrl_deseqres_sig$identifier, "ISEL_sig","not_var_sig"))))
cpmb$bothsig <- factor(cpmb$bothsig, levels = c("not_var_sig","ISEL_Mean_sig","PSS_all_mean_sig","bothsig"))
return(cpmb)
}),data.frame)

