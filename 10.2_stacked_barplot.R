library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)

rm(list=ls())


setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

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

resdegs <- resdegs %>% filter(var %in% vars) %>% filter(!cluster == "C5")
resdegs$gene_clust <-paste0(resdegs$identifier, "_", resdegs$cluster)


pss_genes <- resdegs %>% filter(var == "PSS_all_mean") %>%
  select(identifier, celltype) %>% distinct()

isel_genes <- resdegs %>% filter(var == "ISEL_Mean") %>%
  select(identifier, celltype) %>% distinct()

shared <- inner_join(pss_genes, isel_genes, by = c("identifier", "celltype"))
pss_only <- anti_join(pss_genes, shared, by = c("identifier", "celltype"))
isel_only <- anti_join(isel_genes, shared, by = c("identifier", "celltype"))

count_df <- bind_rows(
  shared %>% mutate(group = "Shared"),
  pss_only %>% mutate(group = "Psychological Stress only"),
  isel_only %>% mutate(group = "Social Support only")
) %>%
  group_by(celltype, group) %>%
  summarise(count = n(), .groups = "drop")

#count_df$group <- factor(count_df$group, levels= c("Social Support only", "Psychological Stress only", "Shared"))
count_df$group <- factor(count_df$group, levels= c("Shared", "Psychological Stress only", "Social Support only"))

celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")

count_df$group <- factor(count_df$group, levels= c("Shared", "Psychological Stress only", "Social Support only"))
count_df$celltype <- factor(count_df$celltype, levels= c("R3 Monocyte", "R0 T CD4+", "R1 T CD8+", "R2 NK", "R4 B"))##, "R6 DC"))

group_colors <- c("Psychological Stress only" = "#E34234", "Social Support only" = "#3A84D8", "Shared" = "#2B8C8C")  # you can adjust the Shared color as needed

# green 228B22
# pastel purple: B07BA3
# torq: 2B8C8C

fig0 <- ggplot(count_df, aes(x = count, y = celltype, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("DEGs Overlap Between Psychological Stress and Social Support", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
#png("nDEGs_2var_barplot_all_celltypesv3.png", width = 1600, height = 1000, res = 240)
png(paste0(outFolder, "/02.1-stacked_bar_plot_unique_shared_ISEL_PSS_DEGs_09-12.png"), width = 1200, height = 700, res = 240)
print(fig0)
dev.off()

pdf(paste0(outFolder, "/02.1-stacked_bar_plot_unique_shared_ISEL_PSS_DEGs_09-12.pdf"), width=5, height=3)
print(fig0)
dev.off()

fname <- paste0(outFolder, "/02.1-number_unique_shared_DEGs_ISEL_PSS.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)




##################################################
#################### for DAMs ####################



setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")

timestamp()

rm(list=ls())

#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz"
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/1.2_ArchR_process/2_Integrate_output/harmony_default/TableS_summ.tsv"

fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz"
dams <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
sigdams <- dams %>% filter(padj_t <0.1)

clusters <- c("C0", "C1", "C2", "C3", "C4", "C5")#, "C6")
vars <- c("ISEL_Mean", "PSS_all_mean")


dams$celltype <- "NA"
dams$celltype[dams$Cluster == "C0"] <- "A0 T CD4+"
dams$celltype[dams$Cluster == "C5"] <- "A5 T CD4+"
dams$celltype[dams$Cluster == "C6"] <- "A6 T CD4+"
dams$celltype[dams$Cluster == "C7"] <- "A7 T CD4+"
dams$celltype[dams$Cluster == "C1"] <- "A1 T CD8+"
dams$celltype[dams$Cluster == "C2"] <- "A2 NK"
dams$celltype[dams$Cluster == "C3"] <- "A3 Monocyte"
dams$celltype[dams$Cluster == "C4"] <- "A4 B"
dams$celltype[dams$Cluster == "C8"] <- "A8 DC"


vars <- c("ISEL_Mean", "PSS_all_mean")
dams <- dams %>% filter(psycho_variable %in% vars) %>% filter(padj_t < 0.1) #%>% filter(!Cluster == "A9")
#dams <- dams %>% filter(padj_t < 0.1)
dim(dams) #189
length(unique(dams$gene)) #158


#resdegs <- resdegs %>% filter(var %in% vars)


pss_genes <- dams %>% filter(psycho_variable == "PSS_all_mean") %>%
  select(gene, celltype) %>% distinct()

isel_genes <- dams %>% filter(psycho_variable == "ISEL_Mean") %>%
  select(gene, celltype) %>% distinct()

shared <- inner_join(pss_genes, isel_genes, by = c("gene", "celltype"))
pss_only <- anti_join(pss_genes, shared, by = c("gene", "celltype"))
isel_only <- anti_join(isel_genes, shared, by = c("gene", "celltype"))

count_df <- bind_rows(
  shared %>% mutate(group = "Shared"),
  pss_only %>% mutate(group = "Psychological Stress only"),
  isel_only %>% mutate(group = "Social Support only")
) %>%
  group_by(celltype, group) %>%
  summarise(count = n(), .groups = "drop")

#count_df$group <- factor(count_df$group, levels= c("Social Support only", "Psychological Stress only", "Shared"))
count_df$group <- factor(count_df$group, levels= c("Shared", "Psychological Stress only", "Social Support only"))

#celltype <- c("R0 CD4+", "R1 CD4+", "R2 CD8+", "R3 NK", "R4 Monocyte", "R5 B", "R6 DC")
celltype <- c("A0 T CD4+", "A5 T CD4+", "A6 T CD4+", "A7 T CD4+", "A1 T CD8+", "A2 NK", "A3 Monocyte", "A4 B")#, "A8 Dendritic cell")

count_df$group <- factor(count_df$group, levels= c("Shared", "Psychological Stress only", "Social Support only"))
count_df$celltype <- factor(count_df$celltype, levels= c("A0 T CD4+", "A3 Monocyte", "A2 NK", "A6 T CD4+","A5 T CD4+",  "A1 T CD8+", "A7 T CD4+", "A4 B"))

#group_colors <- c("Psychological Stress only" = "#E34234", "Social Support only" = "#3A84D8", "Shared" = "#228B22")  # you can adjust the Shared color as needed
group_colors <- c("Psychological Stress only" = "#E34234", "Social Support only" = "#3A84D8", "Shared" = "#2B8C8C")  # you can adjust the Shared color as needed

fig0 <- ggplot(count_df, aes(x = count, y = celltype, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  labs(
    x = "Number of DAMs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("DAMs Overlap Between Psychological Stress and Social Support", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
#png("nDEGs_2var_barplot_all_celltypesv3.png", width = 1600, height = 1000, res = 240)
png(paste0(getwd(), "/09.1.stacked_bar_plot_unique_shared_ISEL_PSS_DAMs_09-12.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()


fname <- paste0(getwd(), "/09.2-number_unique_shared_DAMs_ISEL_PSS.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


##############################################################
########### stacked bar plot SES covariate
##########################################################

rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_SEScov" 

#A
### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs1 <- resall1 %>% filter(padj < 0.1) %>% filter(!cluster == "C5")


fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_SEScov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(padj < 0.1) %>% filter(!cluster == "C5")



clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
#vars <- c("ISEL_Mean")#, "PSS_all_mean")
vars <- c("PSS_all_mean")

resdegs1$celltype <- "NA"
resdegs1$celltype[resdegs1$cluster == "C0"] <- "R0 T CD4+"
resdegs1$celltype[resdegs1$cluster == "C1"] <- "R1 T CD8+"
resdegs1$celltype[resdegs1$cluster == "C2"] <- "R2 NK"
resdegs1$celltype[resdegs1$cluster == "C3"] <- "R3 Monocyte"
resdegs1$celltype[resdegs1$cluster == "C4"] <- "R4 B"
#resdegs1$celltype[resdegs1$cluster == "C5"] <- "R6 DC"

resdegs2$celltype <- "NA"
resdegs2$celltype[resdegs2$cluster == "C0"] <- "R0 T CD4+"
resdegs2$celltype[resdegs2$cluster == "C1"] <- "R1 T CD8+"
resdegs2$celltype[resdegs2$cluster == "C2"] <- "R2 NK"
resdegs2$celltype[resdegs2$cluster == "C3"] <- "R3 Monocyte"
resdegs2$celltype[resdegs2$cluster == "C4"] <- "R4 B"
#resdegs2$celltype[resdegs2$cluster == "C5"] <- "R6 DC"


resdegs1 <- resdegs1 %>% filter(var %in% vars)
resdegs2 <- resdegs2 %>% filter(var %in% vars)

resdegs1$model <- "no ses"
resdegs2$model <- "ses cov"

resdegs1$gene_clust <- paste0(resdegs1$identifier, "_", resdegs1$cluster)
resdegs2$gene_clust <- paste0(resdegs2$identifier, "_", resdegs2$cluster)

resdegsall <- rbind(resdegs1, resdegs2)

noisel_genes <- resdegsall %>% filter(model == "no ses") %>%
  select(identifier, celltype, gene_clust) %>% distinct()

covisel_genes <- resdegsall %>% filter(model == "ses cov") %>%
  select(identifier, celltype, gene_clust) %>% distinct()

shared <- inner_join(noisel_genes, covisel_genes, by = c("identifier", "celltype"))
noisel_only <- anti_join(noisel_genes, covisel_genes, by = c("identifier", "celltype"))
covisel_only <- anti_join(covisel_genes, noisel_genes, by = c("identifier", "celltype"))



count_df <- bind_rows(
  shared %>% mutate(group = "Shared"),
  noisel_only %>% mutate(group = "Default model only"),
  covisel_only %>% mutate(group = "SES corrected only")
) %>%
  group_by(celltype, group) %>%
  summarise(count = n(), .groups = "drop")



#count_df$group <- factor(count_df$group, levels= c("Shared", "Psychological Stress only", "Social Support only"))
#count_df$celltype <- factor(count_df$celltype, levels= c("R3 Monocyte", "R0 T CD4+", "R1 T CD8+", "R2 NK", "R4 B"))##, "R6 DC"))
#group_colors <- c("Psychological Stress only" = "#E34234", "Social Support only" = "#3A84D8", "Shared" = "#2B8C8C")  # you can adjust the Shared color as needed

#count_df$group <- factor(count_df$group, levels= c("Social Support only", "Psychological Stress only", "Shared"))
count_df$group <- factor(count_df$group, levels= c("Shared", "Default model only", "SES corrected only"))

#celltype <- c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B")#, "R6 DC")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")

count_df$group <- factor(count_df$group, levels= c("Shared", "Default model only", "SES corrected only"))
count_df$celltype <- factor(count_df$celltype, levels= c("R3 Monocyte", "R0 T CD4+", "R1 T CD8+", "R2 NK", "R4 B"))##, "R6 DC"))

#group_colors <- c("Default model only" = "#E34234", "SES corrected only" = "#34E3D5", "Shared" = "#228B22")  # you can adjust the Shared color as needed
group_colors <- c("Default model only" = "#D53E4F", "SES corrected only" = "#008080", "Shared" = "#6A3D9A")  # you can adjust the Shared color as needed

fig0 <- ggplot(count_df, aes(x = count, y = celltype, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("Psychological Stress DEGs Overlap Between SES corrected and defualt models", width = 30)
    #title = str_wrap("Social Support DEGs Overlap Between SES corrected and defualt models", width = 30)

  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
#png("nDEGs_2var_barplot_all_celltypesv3.png", width = 1600, height = 1000, res = 240)
png(paste0(getwd(), "/04.1-stacked_bar_plot_PSS_SEScorrected_noSES_v2.png"), width = 1200, height = 700, res = 240)
#png(paste0(getwd(), "/03.1-stacked_bar_plot_ISEL_SEScorrected_noSES_v2.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

fname <- paste0(getwd(), "/04.2-number_unique_shared_PSS_DEGs_SEScorrected_noSES_v2.txt")
#fname <- paste0(getwd(), "/03.2-number_unique_shared_ISEL_DEGs_SEScorrected_noSES_v2.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


degs2 <- resdegs2 %>% dplyr::select(identifier, celltype, padj, logFC)
fname <- paste0(getwd(), "/04.3-nDEGs_PSS_SEScorrected_noSES_v2.txt")
#fname <- paste0(getwd(), "/03.3-nDEGs_ISEL_SEScorrected_noSES_v2.txt")
write.table(degs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


 shared_percent <- as.data.frame(resdegsall) %>%
  distinct(gene_clust, model) %>%
  count(gene_clust) %>%
  summarise(
    shared = sum(n == 2),
    Default_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "no ses") %>% pull(gene_clust))),
    SEScov_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "ses cov") %>% pull(gene_clust))),
    total_Default_DEGs = dim(resdegs1)[1],
    percent_shared = round(shared / total_Default_DEGs * 100, 1),
    percent_Default_only = round(Default_only / total_Default_DEGs * 100, 1),
    percent_SEScov_only = round(SEScov_only / total_Default_DEGs * 100, 1)
  )

#fname <- paste0(getwd(), "/03.4-percentage_shared_ISEL_nDEGs_SEScorrected_noSES_v2.txt")
fname <- paste0(getwd(), "/04.2-percentage_shared_PSS_nDEGs_SEScorrected_noSES_v2.txt")
write.table(shared_percent, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



##############################################################
########### stacked bar plot WHR covariate
##########################################################

rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_WHRcov" 

#A
### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs1 <- resall1 %>% filter(padj < 0.1) %>% filter(!cluster == "C5")


fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_WHRcov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(padj < 0.1) %>% filter(!cluster == "C5")



clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("PSS_all_mean")#"ISEL_Mean")#, "PSS_all_mean")
#vars <- c("ISEL_Mean")#, "PSS_all_mean")


resdegs1$celltype <- "NA"
resdegs1$celltype[resdegs1$cluster == "C0"] <- "R0 T CD4+"
resdegs1$celltype[resdegs1$cluster == "C1"] <- "R1 T CD8+"
resdegs1$celltype[resdegs1$cluster == "C2"] <- "R2 NK"
resdegs1$celltype[resdegs1$cluster == "C3"] <- "R3 Monocyte"
resdegs1$celltype[resdegs1$cluster == "C4"] <- "R4 B"
#resdegs1$celltype[resdegs1$cluster == "C5"] <- "R6 DC"

resdegs2$celltype <- "NA"
resdegs2$celltype[resdegs2$cluster == "C0"] <- "R0 T CD4+"
resdegs2$celltype[resdegs2$cluster == "C1"] <- "R1 T CD8+"
resdegs2$celltype[resdegs2$cluster == "C2"] <- "R2 NK"
resdegs2$celltype[resdegs2$cluster == "C3"] <- "R3 Monocyte"
resdegs2$celltype[resdegs2$cluster == "C4"] <- "R4 B"
#resdegs2$celltype[resdegs2$cluster == "C5"] <- "R6 DC"


resdegs1 <- resdegs1 %>% filter(var %in% vars)
resdegs2 <- resdegs2 %>% filter(var %in% vars)

resdegs1$model <- "no whr"
resdegs2$model <- "whr cov"


resdegs1$gene_clust <- paste0(resdegs1$identifier, "_", resdegs1$cluster)
resdegs2$gene_clust <- paste0(resdegs2$identifier, "_", resdegs2$cluster)

resdegsall <- rbind(resdegs1, resdegs2)


noisel_genes <- resdegsall %>% filter(model == "no whr") %>%
  select(identifier, celltype, gene_clust) %>% distinct()

covisel_genes <- resdegsall %>% filter(model == "whr cov") %>%
  select(identifier, celltype, gene_clust) %>% distinct()

shared <- inner_join(noisel_genes, covisel_genes, by = c("identifier", "celltype"))
noisel_only <- anti_join(noisel_genes, covisel_genes, by = c("identifier", "celltype"))
covisel_only <- anti_join(covisel_genes, noisel_genes, by = c("identifier", "celltype"))

count_df <- bind_rows(
  shared %>% mutate(group = "Shared"),
  noisel_only %>% mutate(group = "Default model only"),
  covisel_only %>% mutate(group = "WHR corrected only")
) %>%
  group_by(celltype, group) %>%
  summarise(count = n(), .groups = "drop")

#count_df$group <- factor(count_df$group, levels= c("Social Support only", "Psychological Stress only", "Shared"))
count_df$group <- factor(count_df$group, levels= c("Shared", "Default model only", "WHR corrected only"))

#celltype <- c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B")#, "R6 DC")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")


count_df$group <- factor(count_df$group, levels= c("Shared", "Default model only", "WHR corrected only"))
count_df$celltype <- factor(count_df$celltype, levels= c("R3 Monocyte", "R0 T CD4+", "R1 T CD8+", "R2 NK", "R4 B"))##, "R6 DC"))

group_colors <- c("Default model only" = "#D53E4F", "WHR corrected only" = "#008080", "Shared" = "#6A3D9A")  # you can adjust the Shared color as needed
#group_colors <- c("Default model only" = "#D53E4F", "SES corrected only" = "#008080", "Shared" = "#6A3D9A")  # you can adjust the Shared color as needed

fig0 <- ggplot(count_df, aes(x = count, y = celltype, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("Psychological Stress DEGs Overlap Between WHR corrected and defualt models", width = 30)
    #title = str_wrap("Social Support DEGs Overlap Between WHR corrected and defualt models", width = 30)

  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
#png("nDEGs_2var_barplot_all_celltypesv3.png", width = 1600, height = 1000, res = 240)
png(paste0(getwd(), "/06.1-stacked_bar_plot_PSS_WHRcorrected_v2.png"), width = 1200, height = 700, res = 240)
#png(paste0(getwd(), "/05.1-stacked_bar_plot_ISEL_WHRcorrected_v2.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

fname <- paste0(getwd(), "/06.2-number_unique_shared_PSS_DEGs_WHRcorrected_v2.txt")
#fname <- paste0(getwd(), "/05.2-number_unique_shared_ISEL_DEGs_WHRcorrected_v2.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


degs2 <- resdegs2 %>% dplyr::select(identifier, celltype, padj, logFC)
fname <- paste0(getwd(), "/06.3-nDEGs_PSS_WHRcorrected_v2.txt")
#fname <- paste0(getwd(), "/05.3-nDEGs_ISEL_WHRcorrected_v2.txt")
write.table(degs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


 shared_percent <- resdegsall %>%
  distinct(gene_clust, model) %>%
  count(gene_clust) %>%
  summarise(
    shared = sum(n == 2),
    Default_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "no whr") %>% pull(gene_clust))),
    WHRcov_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "whr cov") %>% pull(gene_clust))),
    total_Default_DEGs = dim(resdegs1)[1],
    percent_shared = round(shared / total_Default_DEGs * 100, 1),
    percent_Default_only = round(Default_only / total_Default_DEGs * 100, 1),
    percent_WHRcov_only = round(WHRcov_only / total_Default_DEGs * 100, 1)
  )

fname <- paste0(getwd(), "/06.4-percentage_shared_PSS_nDEGs_WHRcorrected_noWHR_v2.txt")
#fname <- paste0(getwd(), "/05.4-percentage_shared_ISEL_nDEGs_WHRcorrected_noWHR_v2.txt")
write.table(shared_percent, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)




##############################################################
########### stacked bar plot ISEL covariate # modify to run it for indepednent filtering off
##########################################################

rm(list=ls())


setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

timestamp()

platePrefix <- "ALL.0.1.13." 
plateSufix1 <- "-RNA-CTRL.noCombat_DESeq.noindfilt" 
plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_iselcov.noindfilt" 

#plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
#plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_iselcov" 

### all DEGs 
fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq.noindfilt.txt", sep = "")
#fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs1 <- resall1 %>% filter(padj < 0.1) %>% filter(!cluster == "C5")


fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_iselcov.noindfilt.txt", sep = "")
#fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_iselcov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(padj < 0.1) %>% filter(!cluster == "C5")

#fname=paste0(dataDir, "/", "nDEGs_PSS_ISEL_20cellfilter_PSScov.noindfilt.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_20cellfilter_PSScov.noindfilt.txt")
#fname=paste0(dataDir, "/", "nDEGs_PSS_ISEL_20cellfilter_iselcov.noindfilt.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_20cellfilter_iselcov.noindfilt.txt")
#nDEGs_PSS_ISEL_20cellfilter.noindfilt.txt
#allres_PSS_ISEL_20cellfilter.noindfilt.txt



clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("PSS_all_mean")#"ISEL_Mean")#, "PSS_all_mean")


resdegs1$celltype <- "NA"
resdegs1$celltype[resdegs1$cluster == "C0"] <- "R0 T CD4+"
resdegs1$celltype[resdegs1$cluster == "C1"] <- "R1 T CD8+"
resdegs1$celltype[resdegs1$cluster == "C2"] <- "R2 NK"
resdegs1$celltype[resdegs1$cluster == "C3"] <- "R3 Monocyte"
resdegs1$celltype[resdegs1$cluster == "C4"] <- "R4 B"
#resdegs1$celltype[resdegs1$cluster == "C5"] <- "R6 DC"

resdegs2$celltype <- "NA"
resdegs2$celltype[resdegs2$cluster == "C0"] <- "R0 T CD4+"
resdegs2$celltype[resdegs2$cluster == "C1"] <- "R1 T CD8+"
resdegs2$celltype[resdegs2$cluster == "C2"] <- "R2 NK"
resdegs2$celltype[resdegs2$cluster == "C3"] <- "R3 Monocyte"
resdegs2$celltype[resdegs2$cluster == "C4"] <- "R4 B"
#resdegs2$celltype[resdegs2$cluster == "C5"] <- "R6 DC"


resdegs1 <- resdegs1 %>% filter(var %in% vars)
resdegs2 <- resdegs2 %>% filter(var %in% vars)

resdegs1$model <- "no isel"
resdegs2$model <- "isel cov"


resdegs1$gene_clust <- paste0(resdegs1$identifier, "_", resdegs1$cluster)
resdegs2$gene_clust <- paste0(resdegs2$identifier, "_", resdegs2$cluster)

resdegsall <- rbind(resdegs1, resdegs2)


noisel_genes <- resdegsall %>% filter(model == "no isel") %>%
  select(identifier, celltype,gene_clust) %>% distinct()

covisel_genes <- resdegsall %>% filter(model == "isel cov") %>%
  select(identifier, celltype, gene_clust) %>% distinct()

shared <- inner_join(noisel_genes, covisel_genes, by = c("identifier", "celltype"))
noisel_only <- anti_join(noisel_genes, covisel_genes, by = c("identifier", "celltype"))
covisel_only <- anti_join(covisel_genes, noisel_genes, by = c("identifier", "celltype"))

count_df <- bind_rows(
  shared %>% mutate(group = "Shared"),
  noisel_only %>% mutate(group = "Default model only"),
  covisel_only %>% mutate(group = "ISEL corrected only")
) %>%
  group_by(celltype, group) %>%
  summarise(count = n(), .groups = "drop")

#count_df$group <- factor(count_df$group, levels= c("Social Support only", "Psychological Stress only", "Shared"))
count_df$group <- factor(count_df$group, levels= c("Shared", "Default model only", "ISEL corrected only"))

#celltype <- c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B")#, "R6 DC")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")


count_df$group <- factor(count_df$group, levels= c("Shared", "Default model only", "ISEL corrected only"))
#count_df$celltype <- factor(count_df$celltype, levels= c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B"))#, "R6 DC"))
count_df$celltype <- factor(count_df$celltype, levels= c("R3 Monocyte", "R0 T CD4+", "R1 T CD8+", "R2 NK", "R4 B"))##, "R6 DC"))

group_colors <- c("Default model only" = "#D53E4F", "ISEL corrected only" = "#008080", "Shared" = "#6A3D9A")  # you can adjust the Shared color as needed
#group_colors <- c("Default model only" = "#D53E4F", "SES corrected only" = "#008080", "Shared" = "#6A3D9A")  # you can adjust the Shared color as needed

fig0 <- ggplot(count_df, aes(x = count, y = celltype, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("Psychological Stress DEGs Overlap Between ISEL corrected and defualt models", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
#png("nDEGs_2var_barplot_all_celltypesv3.png", width = 1600, height = 1000, res = 240)
png(paste0(getwd(), "/09.1-stacked_bar_plot_unique_shared_PSS_noISEL_covISEL_DEGs_05-2025_noindfit.png"), width = 1200, height = 700, res = 240)
#png(paste0(getwd(), "/07.1-stacked_bar_plot_unique_shared_PSS_noISEL_covISEL_DEGs_05-2025_v2.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

fname <- paste0(getwd(), "/09.2-number_unique_shared_PSS_DEGs_noISEL_covISEL_noindfit.txt")
#fname <- paste0(getwd(), "/07.2-number_unique_shared_PSS_DEGs_noISEL_covISEL_v2.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


degs2 <- resdegs2 %>% dplyr::select(identifier, celltype, padj, logFC)
fname <- paste0(getwd(), "/09.3-nDEGs_PSS_iselcov_noindfit.txt")
#fname <- paste0(getwd(), "/07.3-nDEGs_PSS_iselcov_v2.txt")
write.table(degs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


 shared_percent <- resdegsall %>%
  distinct(gene_clust, model) %>%
  count(gene_clust) %>%
  summarise(
    shared = sum(n == 2),
    Default_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "no isel") %>% pull(gene_clust))),
    ISELcov_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "isel cov") %>% pull(gene_clust))),
    total_Default_DEGs = dim(resdegs1)[1],
    percent_shared = round(shared / total_Default_DEGs * 100, 1),
    percent_Default_only = round(Default_only / total_Default_DEGs * 100, 1),
    percent_ISELcov_only = round(ISELcov_only / total_Default_DEGs * 100, 1)
  )

fname <- paste0(getwd(), "/09.4-percentage_shared_PSS_nDEGs_ISELcorrected_noindfit.txt")
#fname <- paste0(getwd(), "/07.4-percentage_shared_PSS_nDEGs_ISELcorrected_v2.txt")
write.table(shared_percent, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



##############################################################
########### stacked bar plot PSS covariate # modify to run it for indepednent filtering off. reciprocal of ISEL 
##########################################################


rm(list=ls())


setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

timestamp()

platePrefix <- "ALL.0.1.13." 
plateSufix1 <- "-RNA-CTRL.noCombat_DESeq.noindfilt" 
plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_PSScov.noindfilt" 

#plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
#plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_PSScov" 

### all DEGs 
fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq.noindfilt.txt", sep = "")
#fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs1 <- resall1 %>% filter(padj < 0.1) %>% filter(!cluster == "C5")

fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_PSScov.noindfilt.txt", sep = "")
#fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_PSScov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(padj < 0.1) %>% filter(!cluster == "C5")

#fname=paste0(dataDir, "/", "nDEGs_PSS_ISEL_20cellfilter_PSScov.noindfilt.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_20cellfilter_PSScov.noindfilt.txt")
#fname=paste0(dataDir, "/", "nDEGs_PSS_ISEL_20cellfilter_iselcov.noindfilt.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_20cellfilter_iselcov.noindfilt.txt")
#nDEGs_PSS_ISEL_20cellfilter.noindfilt.txt
#allres_PSS_ISEL_20cellfilter.noindfilt.txt



clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("ISEL_Mean")



resdegs1$celltype <- "NA"
resdegs1$celltype[resdegs1$cluster == "C0"] <- "R0 T CD4+"
resdegs1$celltype[resdegs1$cluster == "C1"] <- "R1 T CD8+"
resdegs1$celltype[resdegs1$cluster == "C2"] <- "R2 NK"
resdegs1$celltype[resdegs1$cluster == "C3"] <- "R3 Monocyte"
resdegs1$celltype[resdegs1$cluster == "C4"] <- "R4 B"
#resdegs1$celltype[resdegs1$cluster == "C5"] <- "R6 DC"

resdegs2$celltype <- "NA"
resdegs2$celltype[resdegs2$cluster == "C0"] <- "R0 T CD4+"
resdegs2$celltype[resdegs2$cluster == "C1"] <- "R1 T CD8+"
resdegs2$celltype[resdegs2$cluster == "C2"] <- "R2 NK"
resdegs2$celltype[resdegs2$cluster == "C3"] <- "R3 Monocyte"
resdegs2$celltype[resdegs2$cluster == "C4"] <- "R4 B"
#resdegs2$celltype[resdegs2$cluster == "C5"] <- "R6 DC"


resdegs1 <- resdegs1 %>% filter(var %in% vars)
resdegs2 <- resdegs2 %>% filter(var %in% vars)

resdegs1$model <- "no PSS"
resdegs2$model <- "PSS cov"


resdegs1$gene_clust <- paste0(resdegs1$identifier, "_", resdegs1$cluster)
resdegs2$gene_clust <- paste0(resdegs2$identifier, "_", resdegs2$cluster)

resdegsall <- rbind(resdegs1, resdegs2)


noisel_genes <- resdegsall %>% filter(model == "no PSS") %>%
  select(identifier, celltype) %>% distinct()

covisel_genes <- resdegsall %>% filter(model == "PSS cov") %>%
  select(identifier, celltype) %>% distinct()

shared <- inner_join(noisel_genes, covisel_genes, by = c("identifier", "celltype"))
noisel_only <- anti_join(noisel_genes, covisel_genes, by = c("identifier", "celltype"))
covisel_only <- anti_join(covisel_genes, noisel_genes, by = c("identifier", "celltype"))

count_df <- bind_rows(
  shared %>% mutate(group = "Shared"),
  noisel_only %>% mutate(group = "Default model only"),
  covisel_only %>% mutate(group = "PSS corrected only")
) %>%
  group_by(celltype, group) %>%
  summarise(count = n(), .groups = "drop")

#count_df$group <- factor(count_df$group, levels= c("Social Support only", "Psychological Stress only", "Shared"))
count_df$group <- factor(count_df$group, levels= c("Shared", "Default model only", "PSS corrected only"))

celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")


count_df$group <- factor(count_df$group, levels= c("Shared", "Default model only", "PSS corrected only"))
#count_df$celltype <- factor(count_df$celltype, levels= c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B"))#, "R6 DC"))
count_df$celltype <- factor(count_df$celltype, levels= c("R3 Monocyte", "R0 T CD4+", "R1 T CD8+", "R2 NK", "R4 B"))##, "R6 DC"))

group_colors <- c("Default model only" = "#D53E4F", "PSS corrected only" = "#008080", "Shared" = "#6A3D9A")  # you can adjust the Shared color as needed
#group_colors <- c("Default model only" = "#D53E4F", "SES corrected only" = "#008080", "Shared" = "#6A3D9A")  # you can adjust the Shared color as needed

fig0 <- ggplot(count_df, aes(x = count, y = celltype, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("Social Support DEGs Overlap Between PSS corrected and defualt models", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
#png("nDEGs_2var_barplot_all_celltypesv3.png", width = 1600, height = 1000, res = 240)
#png(paste0(getwd(), "/13.1-stacked_bar_plot_unique_shared_ISEL_noPSS_covPSS_DEGs_05-2025.png"), width = 1200, height = 700, res = 240)
png(paste0(getwd(), "/10.1-stacked_bar_plot_unique_shared_ISEL_noPSS_covPSS_DEGs_noindfilt.png"), width = 1200, height = 700, res = 240)
#png(paste0(getwd(), "/08.1-stacked_bar_plot_unique_shared_ISEL_noPSS_covPSS_DEGs_v2.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

fname <- paste0(getwd(), "/10.2-number_unique_shared_DEGs_ISEL_noPSS_covPSS_noindfilt.txt")
#fname <- paste0(getwd(), "/08.2-number_unique_shared_DEGs_ISEL_noPSS_covPSS_v2.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


degs2 <- resdegs2 %>% dplyr::select(identifier, celltype, padj, logFC)
fname <- paste0(getwd(), "/10.3-nDEGs_ISEL_covPSS_noindfilt.txt")
#fname <- paste0(getwd(), "/08.3-nDEGs_ISEL_covPSS_v2.txt")
write.table(degs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



 shared_percent <- resdegsall %>%
  distinct(gene_clust, model) %>%
  count(gene_clust) %>%
  summarise(
    shared = sum(n == 2),
    Default_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "no PSS") %>% pull(gene_clust))),
    PSScov_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "PSS cov") %>% pull(gene_clust))),
    total_Default_DEGs = dim(resdegs1)[1],
    percent_shared = round(shared / total_Default_DEGs * 100, 1),
    percent_Default_only = round(Default_only / total_Default_DEGs * 100, 1),
    percent_PSScov_only = round(PSScov_only / total_Default_DEGs * 100, 1)
  )

fname <- paste0(getwd(), "/10.4-percentage_shared_ISEL_nDEGs_PSScorrected_noindfilt.txt")
#fname <- paste0(getwd(), "/08.4-percentage_shared_ISEL_nDEGs_PSScorrected_v2.txt")
write.table(shared_percent, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


