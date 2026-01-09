
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(reshape)
library(stringr)
library(data.table)
library(ggrepel)



#################################
### get all sig DEGs
###################################
rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"
outDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"

# load psych DEGS
fname="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/allres_allvars_noCombat_DESeq.txt"
res <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")

resdegs <- res %>% filter(padj < 0.1)

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("ISEL_Mean", "PSS_all_mean")

resdegs <- resdegs %>% filter(var %in% vars)

# load LPS degs
fname=paste0("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/treatment_noCOMBAT/deseqres/ALL.0.1.13.deseqres.RNA-LPS_vs_RNA-CTRL.treatment_noCOMBAT.txt")
#fname=paste0("/rs/rs_grp_schold/CZI/RNA/analysis_bk/fastdemux_pseudobulk_ctrl/nodex/treatment_noCOMBAT/deseqres/ALL.0.15.13.deseqres.RNA-LPS_vs_RNA-CTRL.treatment_noCOMBAT.txt")
treats <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")

treats <- treats %>% filter(padj < 0.1)
table(treats$cluster)
subset(treats, grepl("IRF4", identifier, ignore.case = TRUE))
subset(treats, grepl("CXCL10", identifier, ignore.case = TRUE))

treats$celltype <- "NA"
treats$celltype[treats$cluster == "C0"] <- "R0 T CD4+"
treats$celltype[treats$cluster == "C1"] <- "R1 T CD8+"
treats$celltype[treats$cluster == "C2"] <- "R2 NK"
treats$celltype[treats$cluster == "C3"] <- "R3 Monocyte"
treats$celltype[treats$cluster == "C4"] <- "R4 B"
treats$celltype[treats$cluster == "C5"] <- "R5 DC"

treats <- treats %>% filter(!celltype=="R5 DC")

# classify direction in resdegs (PSS / ISEL)
resdegs2 <- resdegs %>%
  mutate(Direction = ifelse(logFC > 0, "Up", "Down"),
         Variable = ifelse(var == "PSS_all_mean", "Psychological Stress", "Social Support"))

# classify direction in LPS DEGs
treats2 <- treats %>%
  mutate(LPS_dir = ifelse(logFC > 0, "Up-regulated", "Down-regulated")) %>%
  select(identifier, celltype, LPS_dir)

# merge, but keep all resdegs
merged <- resdegs2 %>%
  left_join(treats2, by = c("identifier", "celltype")) %>%
  mutate(LPS_dir = ifelse(is.na(LPS_dir), "Non-significant", LPS_dir))

# count for stacked barplot
plotdata <- merged %>%
  group_by(Variable, Direction, LPS_dir) %>%
  summarise(N = n(), .groups = "drop")

plotdata <- plotdata %>% mutate(Group = paste(Variable, Direction), Group = factor(Group, 
                        levels = c("Psychological Stress Up", 
                                   "Psychological Stress Down", 
                                   "Social Support Up", 
                                   "Social Support Down")))

#plotdata$LPS_dir <- factor(plotdata$LPS_dir, levels= c("Up-regulated", "Down-regulated", "Non-significant"))

# stacked bar plot
fig0 <- ggplot(plotdata, aes(x = Group, y = N, fill = LPS_dir)) +
#fig0 <- ggplot(plotdata, aes(x = arrow, y = N, fill = LPS_dir)) +

  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Up-regulated" = "#e31a1c", 
                               "Down-regulated" = "#1f78b4", 
                               "Non-significant" = "grey70")) +
  labs(x = "", y = "Number of DEGs",
       fill = "LPS response") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))
png(paste0(getwd(), "/11.1_stacked_LPSupdown_PSS-ISEL_updown.png"), width = 1200, height = 1000, res = 240)
print(fig0)
dev.off()


##################### facet by cell type
# count for stacked barplot by celltype
plotdata <- merged %>%
  group_by(celltype, Variable, Direction, LPS_dir) %>%
  summarise(N = n(), .groups = "drop")


#define the order of groups
plotdata <- plotdata %>% 
  mutate(Group = paste(Variable, Direction))#, 
#         Group = factor(Group, 
#                        levels = c("Psychological Stress Up", 
#                                   "Psychological Stress Down", 
#                                   "Social Support Up", 
#                                   "Social Support Down")))


plotdata$arrow <- "NA"
plotdata$arrow[plotdata$Group == "Psychological Stress Up"] <- "↑ Psychological Stress & DEG ↑"
plotdata$arrow[plotdata$Group == "Psychological Stress Down"] <- "↑ Psychological Stress & DEG ↓"
plotdata$arrow[plotdata$Group == "Social Support Up"] <- "↑ Social Support & DEG ↑"
plotdata$arrow[plotdata$Group == "Social Support Down"] <- "↑ Social Support & DEG ↓"


# define the order of groups
plotdata <- plotdata %>% 
  mutate(Group = paste(Variable, Direction), 
         Group = factor(arrow, 
                        levels = c("↑ Psychological Stress & DEG ↑", 
                                   "↑ Psychological Stress & DEG ↓", 
                                   "↑ Social Support & DEG ↑", 
                                   "↑ Social Support & DEG ↓")))

# define the order of groups
#plotdata <- plotdata %>% 
#  mutate(Group = paste(Variable, Direction), 
#         Group = factor(Group, 
#                        levels = c("Psychological Stress Up", 
#                                   "Psychological Stress Down", 
#                                   "Social Support Up", 
#                                   "Social Support Down")))


plotdata <- plotdata %>% filter(!celltype=="R5 DC")

# stacked bar plot faceted by cell type, horizontal orientation
#fig1 <- ggplot(plotdata, aes(y = Group, x = N, fill = LPS_dir)) +
fig1 <- ggplot(plotdata, aes(y = arrow, x = N, fill = LPS_dir)) +

  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Up-regulated" = "#E69F00",  # orange
                               "Down-regulated" = "#56B4E9",  # teal/light blue
                               "Non-significant" = "grey70")) +
  labs(y = "Variable", x = "Number of DEGs",
       fill = "LPS response") +
  facet_wrap(~celltype, ncol = 1, scales = "free_y") +  # one column per facet, free y
   theme_minimal(base_size = 13) +
   # theme_bw() +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        strip.text = element_text(size = 12, color = "black"), 
        legend.position = "bottom", 
         plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
         )+
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

png(paste0(getwd(), "/11.2_stacked_LPSupdown_PSS-ISEL_updown_horiz_facet_arrow_V2.png"), 
    width = 1300, height = 1800, res = 240)  # taller PNG for horizontal facets
print(fig1)
dev.off()


library(ggplot2)
library(cowplot)
library(grid)


## Main plot (no legend)
p_main <- ggplot(plotdata, aes(y = arrow, x = N, fill = LPS_dir)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "Up-regulated" = "#E69F00",
    "Down-regulated" = "#56B4E9",
    "Non-significant" = "grey70"
  )) +
  labs(y = "Variable", x = "Number of DEGs") +
  facet_wrap(~celltype, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    legend.position = "none"
  )

png(paste0(getwd(), "/11.2_stacked_LPSupdown_PSS-ISEL_updown_horiz_facet_arrow_V2.png"), width = 1200, height = 1800, res = 240)
print(p_main)
dev.off()


library(cowplot)
library(grid)
library(ggplot2)
library(grid)
library(gtable)

# Build a plot that DEFINITELY contains a legend

png(paste0(getwd(), "/11.3_legend_V2.png"),
    width = 1200, height = 150, res = 240)

p_leg <- ggplot(
  data.frame(
    LPS_dir = factor(
      c("Up-regulated", "Down-regulated", "Non-significant"),
      levels = c("Up-regulated", "Down-regulated", "Non-significant")
    ),
    x = 1:3, y = 1
  ),
  aes(x, y, fill = LPS_dir)
) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "Up-regulated" = "#E69F00",
      "Down-regulated" = "#56B4E9",
      "Non-significant" = "grey70"
    ),
    name = "LPS response"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(fill = guide_legend(nrow = 1))

# Extract legend safely
g <- ggplotGrob(p_leg)
leg <- gtable_filter(g, "guide-box")
grid.newpage()
grid.draw(leg)
dev.off()


##################################################
#################### for DAMs ####################
rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"
outDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"

#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz"
fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz
#res_motif <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz", header=T)
#motif2 <- sort(unique(res_motif$motif_name))
#motif2 <- sort(unique(res_motif$gene))


dams <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
sigdams <- dams %>% filter(padj_t <0.1)

clusters <- c("C0", "C1", "C2", "C3", "C4", "C5")#, "C6")
vars <- c("ISEL_Mean", "PSS_all_mean")


# assign atac metadata: 
dams$celltype <- "NA"
dams$celltype[dams$Cluster == "C0"] <- "A0 T CD4+"
dams$celltype[dams$Cluster == "C1"] <- "A1 T CD8+"
dams$celltype[dams$Cluster == "C2"] <- "A2 NK"
dams$celltype[dams$Cluster == "C3"] <- "A3 T CD4+"
dams$celltype[dams$Cluster == "C4"] <- "A4 Monocyte"
dams$celltype[dams$Cluster == "C5"] <- "A5 T CD4+"
dams$celltype[dams$Cluster == "C6"] <- "A6 B"
dams$celltype[dams$Cluster == "C7"] <- "A7 T CD4+"
dams$celltype[dams$Cluster == "C8"] <- "A8 T CD4+"
dams$celltype[dams$Cluster == "C9"] <- "A9 T CD4+"
dams$celltype[dams$Cluster == "C10"] <- "A10 DC"


vars <- c("ISEL_Mean", "PSS_all_mean")
dams <- dams %>% filter(psycho_variable %in% vars) %>% filter(padj_t < 0.1) #%>% filter(!Cluster == "A9")
#dams <- dams %>% filter(padj_t < 0.1)
dim(dams) #189
length(unique(dams$gene)) #158


# treatment dams
#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_summary/2_response_DAMs.comb.txt.gz"
#treats <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_modelInd_summary/2_th20_plotData.comb.txt.gz"
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

treats <- treats %>% filter(!celltype=="A8 DC")

treats <- treats %>% filter(padj < 0.1)
#subset(treats, grepl("IRF4", identifier, ignore.case = TRUE))
#subset(treats, grepl("CXCL10", identifier, ignore.case = TRUE))


resdegs2 <- dams %>%
  mutate(Direction = ifelse(estimate > 0, "Up", "Down"),
         Variable = ifelse(psycho_variable == "PSS_all_mean", "Psychological Stress", "Social Support"))

# classify direction in LPS DEGs
treats2 <- treats %>%
  mutate(LPS_dir = ifelse(beta > 0, "Increased Activity", "Decreased Activity")) %>%
  select(gene, celltype, LPS_dir)

# merge, but keep all resdegs
merged <- resdegs2 %>%
  left_join(treats2, by = c("gene", "celltype")) %>%
  mutate(LPS_dir = ifelse(is.na(LPS_dir), "Non-significant", LPS_dir))

# count for stacked barplot
plotdata <- merged %>%
  group_by(Variable, Direction, LPS_dir) %>%
  summarise(N = n(), .groups = "drop")

plotdata <- plotdata %>% mutate(Group = paste(Variable, Direction), Group = factor(Group, 
                        levels = c("Psychological Stress Up", 
                                   "Psychological Stress Down", 
                                   "Social Support Up", 
                                   "Social Support Down")))

#plotdata$LPS_dir <- factor(plotdata$LPS_dir, levels= c("Up-regulated", "Down-regulated", "Non-significant"))

# stacked bar plot
fig0 <- ggplot(plotdata, aes(x = Group, y = N, fill = LPS_dir)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Increased Activity" = "#E69F00", 
                               "Decreased Activity" = "#56B4E9", 
                               "Non-significant" = "grey70")) +
  labs(x = "", y = "Number of DAMs",
       fill = "LPS response") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))
png(paste0(getwd(), "/20.4_stacked_DAMs_LPSupdown_PSS-ISEL_updown.png"), width = 1200, height = 1000, res = 240)
print(fig0)
dev.off()


##################### facet by cell type
# count for stacked barplot by celltype
plotdata <- merged %>%
  group_by(celltype, Variable, Direction, LPS_dir) %>%
  summarise(N = n(), .groups = "drop")

# define the order of groups
plotdata <- plotdata %>% 
  mutate(Group = paste(Variable, Direction), 
         Group = factor(Group, 
                        levels = c("Psychological Stress Up", 
                                   "Psychological Stress Down", 
                                   "Social Support Up", 
                                   "Social Support Down")))

plotdata$arrow <- "NA"
plotdata$arrow[plotdata$Group == "Psychological Stress Up"] <- "↑ Psychological Stress & DAM ↑"
plotdata$arrow[plotdata$Group == "Psychological Stress Down"] <- "↑ Psychological Stress & DAM ↓"
plotdata$arrow[plotdata$Group == "Social Support Up"] <- "↑ Social Support & DAM ↑"
plotdata$arrow[plotdata$Group == "Social Support Down"] <- "↑ Social Support & DAM ↓"


# define the order of groups
plotdata <- plotdata %>% 
  mutate(Group = paste(Variable, Direction), 
         Group = factor(arrow, 
                        levels = c("↑ Psychological Stress & DAM ↑", 
                                   "↑ Psychological Stress & DAM ↓", 
                                   "↑ Social Support & DAM ↑", 
                                   "↑ Social Support & DAM ↓")))


#plotdata <- plotdata %>% filter(!celltype=="A8 DC")
plotdata <- plotdata %>% filter(!celltype %in% c("A3 T CD4+", "A6 B"))

# stacked bar plot faceted by cell type, horizontal orientation
#fig1 <- ggplot(plotdata, aes(y = Group, x = N, fill = LPS_dir)) +
fig1 <- ggplot(plotdata, aes(y = arrow, x = N, fill = LPS_dir)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Increased Activity" = "#E69F00", 
                               "Decreased Activity" = "#56B4E9", 
                               "Non-significant" = "grey70")) +
  labs(y = "Variable", x = "Number of DAMs",
       fill = "LPS response") +
  facet_wrap(~celltype, ncol = 1, scales = "free_y") +  # one column per facet, free y
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        strip.text = element_text(size = 12, color = "black"), 
        legend.position = "none")

png(paste0(getwd(), "/20.6_stacked_LPSupdown_PSS-ISEL_updown_horiz_facet_arrow_v2.png"), 
#png(paste0(getwd(), "/20.6_stacked_LPSupdown_PSS-ISEL_updown_horiz_facet_noA1CD8.png"), 
    width = 1200, height = 1800, res = 240)  # taller PNG for horizontal facets
print(fig1)
dev.off()


library(ggplot2)
library(grid)
library(gtable)

# -------------------------
# Legend-only plot
# -------------------------

png(paste0(getwd(), "/20.6_legend_LPS_response.png"),
    width = 1200, height = 150, res = 240)

p_leg <- ggplot(
  data.frame(
    LPS_dir = factor(
      c("Increased Activity", "Decreased Activity", "Non-significant"),
      levels = c("Increased Activity", "Decreased Activity", "Non-significant")
    ),
    x = 1:3, y = 1
  ),
  aes(x, y, fill = LPS_dir)
) +
  geom_col() +
  scale_fill_manual(
    values = c(
      "Increased Activity" = "#E69F00",
      "Decreased Activity" = "#56B4E9",
      "Non-significant" = "grey70"
    ),
    name = "LPS response"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(fill = guide_legend(nrow = 1))

# Extract and save legend
g <- ggplotGrob(p_leg)
leg <- gtable_filter(g, "guide-box")
grid.newpage()
grid.draw(leg)
dev.off()





    






























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
png(paste0(getwd(), "/04.1-stacked_bar_plot_unique_shared_ISEL_PSS_DEGs_09-12.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

fname <- paste0(getwd(), "/04.1-number_unique_shared_DEGs_ISEL_PSS.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F





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
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.15.13." 
plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_SEScov" 

#A
### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs1 <- resall1 %>% filter(padj < 0.1)


fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_SEScov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(padj < 0.1)



clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("ISEL_Mean")#, "PSS_all_mean")

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
    #title = str_wrap("Psychological Stress DEGs Overlap Between SES corrected and defualt models", width = 30)
    title = str_wrap("Social Support DEGs Overlap Between SES corrected and defualt models", width = 30)

  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
#png("nDEGs_2var_barplot_all_celltypesv3.png", width = 1600, height = 1000, res = 240)
#png(paste0(getwd(), "/10.1.1-stacked_bar_plot_PSS_SEScorrected_noSES.png"), width = 1200, height = 700, res = 240)
png(paste0(getwd(), "/10.2.1-stacked_bar_plot_ISEL_SEScorrected_noSES.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

#fname <- paste0(getwd(), "/10.1.2-number_unique_shared_PSS_DEGs_SEScorrected_noSES.txt")
fname <- paste0(getwd(), "/10.2.2-number_unique_shared_ISEL_DEGs_SEScorrected_noSES.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


degs2 <- resdegs2 %>% dplyr::select(identifier, celltype, padj, logFC)
#fname <- paste0(getwd(), "/10.1.3-nDEGs_PSS_SEScorrected_noSES.txt")
fname <- paste0(getwd(), "/10.2.3-nDEGs_ISEL_SEScorrected_noSES.txt")
write.table(degs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


 shared_percent <- resdegsall %>%
  distinct(gene_clust, model) %>%
  count(gene_clust) %>%
  summarise(
    shared = sum(n == 2),
    no_ses_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "no ses") %>% pull(gene_clust))),
    ses_cov_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "ses cov") %>% pull(gene_clust))),
    total_gene_clust = n(),
    percent_shared = round(shared / total_gene_clust * 100, 1),
    percent_no_ses_only = round(no_ses_only / total_gene_clust * 100, 1),
    percent_ses_cov_only = round(ses_cov_only / total_gene_clust * 100, 1)
  )

fname <- paste0(getwd(), "/10.1.4-percentage_shared_ISEL_nDEGs_SEScorrected_noSES.txt")
#fname <- paste0(getwd(), "/10.2.4-percentage_shared_PSS_nDEGs_SEScorrected_noSES.txt")
write.table(shared_percent, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



##############################################################
########### stacked bar plot WHR covariate
##########################################################

rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.15.13." 
plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_WHRcov" 

#A
### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs1 <- resall1 %>% filter(padj < 0.1)


fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_WHRcov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(padj < 0.1)



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
png(paste0(getwd(), "/11.1.1-stacked_bar_plot_PSS_WHRcorrected.png"), width = 1200, height = 700, res = 240)
#png(paste0(getwd(), "/11.2.1-stacked_bar_plot_ISEL_WHRcorrected.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

fname <- paste0(getwd(), "/11.1.2-number_unique_shared_PSS_DEGs_WHRcorrected.txt")
#fname <- paste0(getwd(), "/11.2.2-number_unique_shared_ISEL_DEGs_WHRcorrected.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


degs2 <- resdegs2 %>% dplyr::select(identifier, celltype, padj, logFC)
fname <- paste0(getwd(), "/11.1.3-nDEGs_PSS_WHRcorrected.txt")
#fname <- paste0(getwd(), "/11.2.3-nDEGs_ISEL_WHRcorrected.txt")
write.table(degs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


 shared_percent <- resdegsall %>%
  distinct(gene_clust, model) %>%
  count(gene_clust) %>%
  summarise(
    shared = sum(n == 2),
    no_whr_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "no whr") %>% pull(gene_clust))),
    whr_cov_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "whr cov") %>% pull(gene_clust))),
    total_gene_clust = n(),
    percent_shared = round(shared / total_gene_clust * 100, 1),
    percent_no_whr_only = round(no_whr_only / total_gene_clust * 100, 1),
    percent_whr_cov_only = round(whr_cov_only / total_gene_clust * 100, 1)
  )

fname <- paste0(getwd(), "/11.1.4-percentage_shared_PSS_nDEGs_WHRcorrected_noWHR.txt")
#fname <- paste0(getwd(), "/11.2.4-percentage_shared_ISEL_nDEGs_WHRcorrected_noWHR.txt")
write.table(shared_percent, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)






##############################################################
########### stacked bar plot ISEL covariate # modify to run it for indepednent filtering off
##########################################################

rm(list=ls())


setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

timestamp()

platePrefix <- "ALL.0.15.13." 
#plateSufix1 <- "-RNA-CTRL.noCombat_DESeq.noindfilt" 
#plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_iselcov.noindfilt" 

plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_iselcov" 

### all DEGs 
#fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq.noindfilt.txt", sep = "")
fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs1 <- resall1 %>% filter(padj < 0.1)

#fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_iselcov.noindfilt.txt", sep = "")
fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_iselcov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(padj < 0.1)

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
png(paste0(getwd(), "/14.1-stacked_bar_plot_unique_shared_PSS_noISEL_covISEL_DEGs_05-2025.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

fname <- paste0(getwd(), "/14.2-number_unique_shared_PSS_DEGs_noISEL_covISEL.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


degs2 <- resdegs2 %>% dplyr::select(identifier, celltype, padj, logFC)
fname <- paste0(getwd(), "/14.3-nDEGs_PSS_iselcov.txt")
write.table(degs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


 shared_percent <- resdegsall %>%
  distinct(gene_clust, model) %>%
  count(gene_clust) %>%
  summarise(
    shared = sum(n == 2),
    no_isel_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "no isel") %>% pull(gene_clust))),
    isel_cov_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "isel cov") %>% pull(gene_clust))),
    total_gene_clust = n(),
    percent_shared = round(shared / total_gene_clust * 100, 1),
    percent_no_isel_only = round(no_isel_only / total_gene_clust * 100, 1),
    percent_isel_cov_only = round(isel_cov_only / total_gene_clust * 100, 1)
  )

#fname <- paste0(getwd(), "/12.1.4-percentage_shared_ISEL_nDEGs_WHRcorrected_noWHR.txt")
fname <- paste0(getwd(), "/14.2.4-percentage_shared_PSS_nDEGs_ISELcorrected.txt")
write.table(shared_percent, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)




##############################################################
########### stacked bar plot PSS covariate # modify to run it for indepednent filtering off. reciprocal of ISEL 
##########################################################


rm(list=ls())


setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

timestamp()

platePrefix <- "ALL.0.15.13." 
#plateSufix1 <- "-RNA-CTRL.noCombat_DESeq.noindfilt" 
#plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_PSScov.noindfilt" 

plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_PSScov" 

### all DEGs 
#fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq.noindfilt.txt", sep = "")
fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs1 <- resall1 %>% filter(padj < 0.1)

#fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_PSScov.noindfilt.txt", sep = "")
fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_PSScov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(padj < 0.1)

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
png(paste0(getwd(), "/15.1-stacked_bar_plot_unique_shared_ISEL_noPSS_covPSS_DEGs_05-2025.png"), width = 1200, height = 700, res = 240)

print(fig0)
dev.off()

fname <- paste0(getwd(), "/15.2-number_unique_shared_DEGs_ISEL_noPSS_covPSS.txt")
write.table(count_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


degs2 <- resdegs2 %>% dplyr::select(identifier, celltype, padj, logFC)
fname <- paste0(getwd(), "/15.3-nDEGs_ISEL_covPSS.txt")
write.table(degs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



 shared_percent <- resdegsall %>%
  distinct(gene_clust, model) %>%
  count(gene_clust) %>%
  summarise(
    shared = sum(n == 2),
    no_pss_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "no PSS") %>% pull(gene_clust))),
    pss_cov_only = sum(n == 1 & gene_clust %in% (resdegsall %>% filter(model == "PSS cov") %>% pull(gene_clust))),
    total_gene_clust = n(),
    percent_shared = round(shared / total_gene_clust * 100, 1),
    percent_no_pss_only = round(no_pss_only / total_gene_clust * 100, 1),
    percent_pss_cov_only = round(pss_cov_only / total_gene_clust * 100, 1)
  )

fname <- paste0(getwd(), "/15.1.4-percentage_shared_ISEL_nDEGs_PSScorrected.txt")
#fname <- paste0(getwd(), "/12.2.4-percentage_shared_PSS_nDEGs_ISELcorrected.txt")
write.table(shared_percent, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



















setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/barplot/")

fname = paste0(getwd(), "/10.2.2-number_unique_shared_ISEL_DEGs_SEScorrected_noSES.txt")
tbl <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

tbls <- tbl %>%
  pivot_wider(names_from = group, values_from = count) %>%
  mutate(
    total_DEGs = `Default model only` + `SES corrected only` + Shared,
    percent_shared = round((Shared / total_DEGs) * 100, 1)
  ) %>% as.data.frame()

tbls

fname <- paste0(getwd(), "/10.2.4-percentage_unique_shared_ISEL_nDEGs_SEScorrected.txt")
write.table(tbls, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

  



overall_percent_shared <- tbl %>%
  group_by(group) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  summarise(
    shared = total[group == "Shared"],
    default_only = total[group == "Default model only"],
    ses_only = total[group == "SES corrected only"],
    total_unique = shared + default_only + ses_only,
    percent_shared = round((shared / total_unique) * 100, 1)
  ) %>% as.data.frame()

  overall_percent_shared




#ALL.0.2.11.C3.deseqres_ISEL_Mean-RNA-LPS.SES_PCs_sex_age_and_treats_generem_PSScov.noindfilt.txt
#ALL.0.2.11.C5.deseqres_PSS_all_mean-RNA-LPS.SES_PCs_sex_age_and_treats_generem_iselcov.noindfilt.txt



