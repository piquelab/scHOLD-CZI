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




#################################################################
################### facet wrapped qq plot  ######################
#################################################################

rm(list=ls())

outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/qqplot/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

#ALL.0.15.13.C0.deseqres_ISEL_Mean-RNA-CTRL.noCombat_DESeq.txt

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/barplot/")
#dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/deseqres/"

#platePrefix <- "ALL.0.2.11." 
#plateSufix <- "-RNA-CTRL.SES_PCs_sex_age_and_treats_generem" 
#plateSufix <- "-RNA-CTRL.SES_PCs_sex_age_and_treats_generem_iselcov" 


variables <- c("PSS_all_mean", "ISEL_Mean")#, "cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
variable_names <- c("Psychological Stress", "Social Support")#, "Cytokines") #"socioeconomic status", "diastolic blood pressure", 

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
#clusters <- c("0","4")

celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")


#col2 <- c("#E69F00","#009E73","#984EA3","#D55E00", "#CC79A7", "#56B4E9")#, "#F0E442", "#999999")
#col2 <- c("#E69F00","#984EA3","#D55E00", "#CC79A7", "#56B4E9")#, "#F0E442", "#999999")

col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#AA4B56")#, "#D4B9DA")

names(col2) <- celltypes


clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")

# Initialize list to store results
all_res <- list()

timestamp()

# Outer loop for variables
for (i in 1:length(variables)) {
  myvar <- variables[i]
  var_name <- variable_names[i]
  
  cat("######################################################################\n")
  cat("## Processing: ", var_name, "\n")

  res_list <- list()

  # Inner loop for clusters
  for (c in 1:length(clusters)) {
    cluster <- clusters[c]
    celltype <- celltypes[c]

    fname <- paste(dataDir, platePrefix, "C", cluster, ".deseqres_", myvar, plateSufix, ".txt", sep="")
    
    # Read data
    res_data <- read.table(fname, sep="\t", header=T, quote='"', comment="")
    res_data$celltype <- celltype
    res_data$Variable <- var_name  # Add variable name for faceting
    res_data <- res_data %>% dplyr::select(-any_of("baseMean"))

    res_list[[paste("cluster_", cluster)]] <- res_data
  }

  # Combine data for the current variable
  res_var <- do.call(rbind, res_list)
  all_res[[var_name]] <- res_var
}

# Combine all variables' data into one dataframe
resall <- do.call(rbind, all_res)

#celltype <- c("CD4+ T cells 1", "CD4+ T cells 2", "CD8+ T cells", "NK cells", "Monocytes", "B cells", "Dentritic Cells")
#resall$celltype <- factor(resall$celltype, levels= c("CD4+ T cells 1", "CD4+ T cells 2", "CD8+ T cells", "NK cells", "Monocytes", "B cells", "Dentritic Cells"))
#resall$Variable <- factor(resall$Variable, levels= c("Social Support", "Psychological Resources", "Perceived Stress", "Cytokines"))

# Define colors
#col2 <- c("#E69F00","#009E73","#984EA3","#D55E00", "#CC79A7", "#56B4E9")#, "#F0E442", "#999999")
#col2 <- c("#E69F00", "#984EA3","#D55E00", "#CC79A7", "#56B4E9")#, "#F0E442", "#999999")
col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#AA4B56")#, "#D4B9DA")

names(col2) <- celltypes

# Compute expected and observed p-values
res <- resall %>%
  filter(!is.na(pvalue)) %>%
  arrange(pvalue) %>%
  group_by(Variable, celltype) %>%
  mutate(r = rank(pvalue, ties.method = "random"), pexp = r / length(pvalue))

res$celltype <- factor(res$celltype, levels = celltypes)
#res$Variable <- factor(res$Variable, levels = c("Psychological Stress", "Social Support", "Cytokines"))
res$Variable <- factor(res$Variable, levels = c("Psychological Stress", "Social Support"))#, "Cytokines"))

# Create QQ plot
p1 <- ggplot(res, aes(x = -log10(pexp), y = -log10(pvalue), color = celltype)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlab(expression(Expected -log[10](p))) +
  ylab(expression(Observed -log[10](p))) +
  facet_wrap(~Variable, nrow = 1) +  # Facet by variable, in one row
  theme_bw() +
     ylim(0, 10) +
  theme(
    axis.text.x = element_text(size = 15),  
    axis.text.y = element_text(size = 15),                        
    text = element_text(size = 15),                              
    plot.title = element_text(size = 20)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = col2)  # Apply custom colors

# Save the combined plot
figfn <- paste0(outFolder, "01.1_qqplot_ISEL_PSS.png")
png(figfn, width = 1500, height = 1000, res = 240)  # Adjust width for side-by-side layout
print(p1)
dev.off()


figfn <- paste0(outFolder, "01.1_qqplot_ISEL_PSS.pdf")
pdf(figfn, width=7, height=7)
print(p1)
dev.off()




#################################################################
###################DARs ######################
#################################################################

rm(list=ls())

outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/qqplot/"
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/qqplot/")


fname="/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/1_DiffPeak.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
res_data <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")

ressig <- res_data %>% filter(p.adjusted < 0.1) %>% filter(psycho_variable %in% c("PSS_all_mean", "ISEL_Mean"))
ressig <- res_data %>% filter(psycho_variable %in% c("PSS_all_mean", "ISEL_Mean"))


# assign atac metadata: 
ressig$celltype <- "NA"
ressig$celltype[ressig$Cluster == "C0"] <- "A0 T CD4+"
ressig$celltype[ressig$Cluster == "C1"] <- "A1 T CD8+"
ressig$celltype[ressig$Cluster == "C2"] <- "A2 NK"
ressig$celltype[ressig$Cluster == "C3"] <- "A3 T CD4+"
ressig$celltype[ressig$Cluster == "C4"] <- "A4 Monocyte"
ressig$celltype[ressig$Cluster == "C5"] <- "A5 T CD4+"
ressig$celltype[ressig$Cluster == "C6"] <- "A6 B"
ressig$celltype[ressig$Cluster == "C7"] <- "A7 T CD4+"
ressig$celltype[ressig$Cluster == "C8"] <- "A8 T CD4+"
ressig$celltype[ressig$Cluster == "C9"] <- "A9 T CD4+"
ressig$celltype[ressig$Cluster == "C10"] <- "A10 DC"

ressig$Variable <- "NA"
ressig$Variable[ressig$psycho_variable == "PSS_all_mean"] <- "Psychological Stress"
ressig$Variable[ressig$psycho_variable == "ISEL_Mean"] <- "Social Support"

ressig <- ressig %>% filter(Cluster %in% c("C0", "C1", "C2", "C4", "C6"))


# Compute expected and observed p-values
res <- ressig %>%
  filter(!is.na(p.value)) %>%
  arrange(p.value) %>%
  group_by(Variable, celltype) %>%
  mutate(r = rank(p.value, ties.method = "random"), pexp = r / length(p.value))



#celltypes <- c("A0 T CD4+", "A1 T CD8+", "A2 NK", "A3 T CD4+", "A4 Monocyte", "A5 T CD4+", "A6 B", "A7 T CD4+")#, "A8 T CD4+", "A9 T CD4+", "A10 DC")#, "R5 DC")
celltypes <- c("A0 T CD4+", "A1 T CD8+", "A2 NK",  "A4 Monocyte", "A6 B")#, "A8 T CD4+", "A9 T CD4+", "A10 DC")#, "R5 DC")

col_cl <- c(
  "A0 T CD4+"   = "#FF7F00",  # orange
  "A1 T CD8+"   = "#E6E600",  # blue
  "A2 NK"      = "#4DAF4A",  # reddish orange
  "A4 Monocyte" = "#984EA3",  # purple/magenta
  "A6 B"       = "#D97986"#,  # light blue
 # "A3 T CD4+"  = darken("#FF7F00", 0.1),   #  teal green
 # "A5 T CD4+"   = lighten("#FF7F00", 0.2),   #  teal green
 # "A7 T CD4+"  = lighten("#FF7F00",0.5)   #  teal green
)
 
res$celltype <- factor(res$celltype, levels = celltypes)
#res$Variable <- factor(res$Variable, levels = c("Psychological Stress", "Social Support", "Cytokines"))
res$Variable <- factor(res$Variable, levels = c("Psychological Stress", "Social Support"))#, "Cytokines"))

# Create QQ plot
p1 <- ggplot(res, aes(x = -log10(pexp), y = -log10(p.value), color = celltype)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlab(expression(Expected -log[10](p))) +
  ylab(expression(Observed -log[10](p))) +
  facet_wrap(~Variable, nrow = 1) +  # Facet by variable, in one row
  theme_bw() +
     ylim(0, 10) +
  theme(
    axis.text.x = element_text(size = 15),  
    axis.text.y = element_text(size = 15),                        
    text = element_text(size = 15),                              
    plot.title = element_text(size = 20)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = col_cl)  # Apply custom colors

# Save the combined plot
figfn <- paste0(outFolder, "01.4_DARs_qqplot_ISEL_PSSv2.png")
png(figfn, width = 1500, height = 1000, res = 240)  # Adjust width for side-by-side layout
print(p1)
dev.off()


figfn <- paste0(outFolder, "01.4_DARs_qqplot_ISEL_PSSv2.pdf")
pdf(figfn, width=7, height=7)
print(p1)
dev.off()
