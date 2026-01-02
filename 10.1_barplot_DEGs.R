
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(data.table)
library(ggrepel)
library(qvalue)

################################################################################################################
################### facet wrapped barplots for 3 conditions, ISEL, PSS, cytocomp  ######################
################################################################################################################

rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 
#plateSufix <- "-RNA-LPS.noCombat_DESeq" 

#variables <- c("cytocomp","z_il6_0_log_w", "z_tnfa_0_log_w", "z_ifny_0_log_w","z_il12p70_0_log_w",    
#              "z_il1b_0_log_w", "z_il2_0_log_w", "z_il8_0_log_w", "z_il4_0_log_w", "z_il13_0_log_w","z_il10_0_log_w")
#variable_names <- c("Cytokines", "IL-6", "TNF-α", "IFN-γ", "IL-12p70", "IL-1β", "IL-2", "IL-8", "IL-4", "IL-13","IL-10")
#variables <- c("z_ifny_0_log_w")#, "cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
#variable_names <- c("IFN-γ")#, "Cytokines") #"socioeconomic status", "diastolic blood pressure", 

variables <- c("PSS_all_mean", "ISEL_Mean")#, "cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
variable_names <- c("Psychological Stress", "Social Support")#, "Cytokines") #"socioeconomic status", "diastolic blood pressure", 

#variables <- c("z_ifny_0_log_w")#, "cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
#variable_names <- c("IFN-γ")#, "Cytokines") #"socioeconomic status", "diastolic blood pressure", 

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R5 DC")

timestamp()
##################################################### bar plot by cell type

# Initialize an empty list to store data for all variables
all_data <- list()
sigsalldata <- list()
allressigs <- list()

# Loop over variables
for (i in seq_along(variables)) {
  myvar <- variables[i]
  var_name <- variable_names[i]
  
  cat("########################################################################################", "\n")
  cat("############### running ",myvar, " #######################", "\n")

  dfvar <- list()
  
  # Loop over clusters
  for (j in seq_along(clusters)) {
    clust <- clusters[j]
    clustname <- celltype[j]
      cat("## cluster ",clust, " var ", myvar, " ##", "\n")

    fname <- paste(dataDir, platePrefix, "C", clust, ".deseqres_", myvar, plateSufix, ".txt", sep = "")
    
    # Read data
    res_data <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")
    res_data$celltype <- clustname
          cat("#nDEGs: ", length(res_data %>% filter(padj<0.1) %>% pull(identifier)) , "\n")
    # Store each data frame
    dfvar[[clust]] <- res_data
  }
  
  # Combine results for all clusters
  res<- bind_rows(dfvar)
  ressig <- res %>% filter(padj < 0.1)

  deg_counts <- table(ressig$cluster)

  # Skip this variable if no cluster has at least 50 DEGs
  #if (all(deg_counts < 50)) {
  #  next
  #}

  # Summarize number of DEGs
  sigs <- ressig %>%
    mutate(direction = ifelse(logFC > 0, "1", "2")) %>%
    group_by(celltype, direction) %>%
    summarise(ngene = n(), .groups = "drop") %>%
    mutate(ngene2 = ifelse(direction == "2", -ngene, ngene),
           comb = paste(celltype, direction, sep = "_"),
           Variable = var_name) # Add variable name
  
  # Store data
  all_data[[i]] <- sigs
  allressigs[[i]] <- ressig


  # Summarize number of DEGs
  sigsall <- ressig %>%
    group_by(celltype) %>%
    summarise(ngene = n(), .groups = "drop") %>%
    mutate(
           Variable = var_name) # Add variable name
    # Store data
  sigsalldata[[i]] <- sigsall
}

all_data <- lapply(all_data, function(df) {
  df %>% mutate(direction = as.character(direction))
})
# Combine all data into one dataframe
plotdata <- bind_rows(all_data)
allplotdata <- bind_rows(sigsalldata)
ressigall <- bind_rows(allressigs)

fname <- paste0(outFolder, "01.1.number_DEGs_up_down_ISEL_PSS.txt")
write.table(plotdata, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

fname <- paste0("01.2.number_DEGs_all_per_ISEL_PSS.txt")
write.table(allplotdata, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

fname <- paste0("01.3.all_ressigs_ISEL_PSS.txt")
write.table(ressigall, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

plotdata$celltype <- factor(plotdata$celltype, levels= c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B"))#, "R5 DC"))

col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#AA4B56")#, "#D4B9DA")

names(col2) <- celltype
col2w <- colorspace::lighten(col2, 0.3)
col2comb <- c(col2, col2w)
names(col2comb) <- paste(celltype, rep(c(1, 2), each = 5), sep = "_")

# Set axis limits
maxy <- max(plotdata$ngene2) + 50
miny <- min(plotdata$ngene2) - 50
breaks_value <- pretty(c(miny, maxy), 10)

#breaks_value <- pretty(c(-600, 850), 10)
#breaks_value <- pretty(c(-100, 100), 10)

# Create plot with facet wrap
fig0 <- ggplot(plotdata, aes(x = celltype, y = ngene2, fill = comb)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col2comb) +
  geom_hline(yintercept = 0, color = "grey60") +
  ylim(1000, 1000) +
  geom_text(aes(label = abs(ngene2), vjust = ifelse(direction == "2", 1.2, -0.2)), size = 3) +
  #scale_y_continuous("Number of DEGs", breaks = breaks_value, labels = function(x) abs(x),  limits = c(-600, 850)) +
  scale_y_continuous("Number of DEGs", breaks = breaks_value, labels = function(x) abs(x)) +#,  limits = c(-100, 100)) +
  theme_bw() +
  xlab("Cell Type") +
  facet_wrap(~Variable, nrow = 1) +  # Facet by variable
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, hjust = 1, vjust = 1, angle = 45),
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size = 12),  
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14))

# Save the plot
png("01.1_nDEGs_PSS_ISEL.png", width = 1300, height = 1000, res = 240)
print(fig0)
dev.off()

pdf("01.1_nDEGs_PSS_ISEL.pdf", width=7, height=7)
print(fig0)
dev.off()

##########################







######################## DARs ##########################

rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/")
outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/"

fname="/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/1_DiffPeak.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
res_data <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")

ressig <- res_data %>% filter(p.adjusted < 0.1) %>% filter(psycho_variable %in% c("PSS_all_mean", "ISEL_Mean"))

# assign atac metadata: 
ressig$celltype <- "NA"
ressig$celltype[ressig$Cluster == "C0"] <- "A0 T CD4+"
ressig$celltype[ressig$Cluster == "C1"] <- "A1 T CD8+"
ressig$celltype[ressig$Cluster == "C2"] <- "A2 NK"
ressig$celltype[ressig$Cluster == "C4"] <- "A4 Monocyte"


ressig$Variable <- "NA"
ressig$Variable[ressig$psycho_variable == "PSS_all_mean"] <- "Psychological Stress"
ressig$Variable[ressig$psycho_variable == "ISEL_Mean"] <- "Social Support"

  # Summarize number of DEGs
  sigs <- ressig %>%
    mutate(direction = ifelse(estimate > 0, "1", "2")) %>%
    group_by(celltype, direction, Variable) %>%
    summarise(ngene = n(), .groups = "drop") %>%
    mutate(ngene2 = ifelse(direction == "2", -ngene, ngene),
           comb = paste(celltype, direction, sep = "_"))#,
           #Variable = Variable) # Add variable name
  




clusters <- c("C0", "C1", "C2", "C4")#, "5")#, "6")
celltype <- c("A0 T CD4+", "A1 T CD8+", "A2 NK", "A4 Monocyte")#, "R5 DC")

plotdata <- sigs

plotdata$celltype <- factor(plotdata$celltype, levels= c("A0 T CD4+", "A1 T CD8+", "A2 NK", "A4 Monocyte"))#, "R5 DC"))

col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3")#, "#AA4B56")#, "#D4B9DA")

names(col2) <- celltype
col2w <- colorspace::lighten(col2, 0.3)
col2comb <- c(col2, col2w)
names(col2comb) <- paste(celltype, rep(c(1, 2), each = 4), sep = "_")

# Set axis limits
maxy <- max(plotdata$ngene2) + 50
miny <- min(plotdata$ngene2) - 50
breaks_value <- pretty(c(miny, maxy), 10)

#breaks_value <- pretty(c(-600, 850), 10)
#breaks_value <- pretty(c(-100, 100), 10)

# Create plot with facet wrap
fig0 <- ggplot(plotdata, aes(x = celltype, y = ngene2, fill = comb)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col2comb) +
  geom_hline(yintercept = 0, color = "grey60") +
  ylim(1000, 1000) +
  geom_text(aes(label = abs(ngene2), vjust = ifelse(direction == "2", 1.2, -0.2)), size = 3) +
  #scale_y_continuous("Number of DEGs", breaks = breaks_value, labels = function(x) abs(x),  limits = c(-600, 850)) +
  scale_y_continuous("Number of DAMs", breaks = breaks_value, labels = function(x) abs(x)) +#,  limits = c(-100, 100)) +
  theme_bw() +
  xlab("Cell Type") +
  facet_wrap(~Variable, nrow = 1) +  # Facet by variable
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, hjust = 1, vjust = 1, angle = 45),
        axis.text.y = element_text(size=12), 
        axis.title.y = element_text(size = 12),  
        axis.title.x = element_blank(),
        strip.text = element_text(size = 14))

# Save the plot
png("01.4_nDAMs_PSS_ISEL.png", width = 1300, height = 1000, res = 240)
print(fig0)
dev.off()

pdf("01.4_nDAMs_PSS_ISEL.pdf", width=7, height=7)
print(fig0)
dev.off()
