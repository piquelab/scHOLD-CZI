library(ggplot2)
library(dplyr)
library(tidyr)
library(DESeq2)
library(data.table)
library(qvalue)
library(tidyverse)
library(pheatmap)

rm(list=ls())


outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/heatmap/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

variables <- c("PSS_all_mean", "ISEL_Mean")#, "cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
variable_names <- c("Psychological Stress", "Social Support")#, "Cytokines") #"socioeconomic status", "diastolic blood pressure", 

# Variables to loop over
#variables <- c("cytocomp", "PSS_all_mean", "isel")#, "cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
#variable_names <- c("Cytokines", "Psychological Stress", "Social Support")#, "Cytokines") #"socioeconomic status", "diastolic blood pressure", 

#fname =  paste(dataDir, platePrefix, "_all_genes_tested_celltype_varname.txt", sep = "")

### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
res$gene_cluster <- paste0(res$cluster, "_", res$identifier)


res$celltype <- "NA"
res$celltype[res$cluster == "C0"] <- "R0 T CD4+"
res$celltype[res$cluster == "C1"] <- "R1 T CD8+"
res$celltype[res$cluster == "C2"] <- "R2 NK"
res$celltype[res$cluster == "C3"] <- "R3 Monocyte"
res$celltype[res$cluster == "C4"] <- "R4 B"
#res$celltype[res$cluster == "C6"] <- "R6 DC"


res$variable[res$var == "ISEL_Mean"] <- "Social Support"
res$variable[res$var == "PSS_all_mean"] <- "Psychological Stress"
#res$variable[res$var == "cytocomp"] <- "Cytokines"
#res$variable[res$var == "BPd_avg"] <- "diastolic blood pressure"
#res$variable[res$var == ""] <- ""
#res$variable[res$var == ""] <- ""

ressub <- res %>% filter(var %in% variables) 

celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")

resdeg <- res %>% filter(padj < 0.1)
resdegsub <- resdeg %>% filter(var %in% variables)
length(unique(resdegsub$identifier)) 
# 2746 unique DEGs when only looking at PSS and isel

resdegsub$gene_cluster <- paste0(resdegsub$cluster, "_", resdegsub$identifier)

uniondegs <- ressub %>% filter(identifier %in% resdegsub$identifier) 

table(resdegsub$var, resdegsub$cluster)

uniondegs$variable_celltype  <- paste0(uniondegs$variable, " - ", uniondegs$celltype)
uniondegs$zscore  <- uniondegs$logFC / uniondegs$SE

#uniondegs <- uniondegs %>% filter(!cluster == "C6") #27207

uniondegs$variable <- factor(uniondegs$variable, levels= variable_names)


############# generate heatmap

# Convert to wide format: Rows = variable_celltype, Columns = Genes
heatmap_data <- uniondegs %>%
  dplyr::select(variable_celltype, identifier, zscore) %>% # Use logFC or another numeric value
  pivot_wider(names_from = identifier, values_from = zscore, values_fill = 0) %>%
  column_to_rownames("variable_celltype")

# Convert to matrix for heatmap plotting
heatmap_matrix <- as.matrix(heatmap_data)
rownames(heatmap_matrix)

 
row_order <- c( #"Cytokines - CD4+ T cells 1", "Cytokines - CD4+ T cells 2", "Cytokines - CD8+ T cells", "Cytokines - NK cells", "Cytokines - Monocytes", "Cytokines - B cells",
              "Psychological Stress - R0 T CD4+", "Psychological Stress - R1 T CD8+", "Psychological Stress - R2 NK", "Psychological Stress - R3 Monocyte", "Psychological Stress - R4 B", #"Psychological Stress - R6 Dendritic cell",
              "Social Support - R0 T CD4+", "Social Support - R1 T CD8+", "Social Support - R2 NK", "Social Support - R3 Monocyte", "Social Support - R4 B")#, "Social Support - R6 Dendritic cell")#, 
                           

heatmap_matrix_ordered <- heatmap_matrix[row_order, colnames(heatmap_matrix)]


###### need to add cell type and variable annotations on the heatmap
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")

 #col2 <- c("#E69F00","#984EA3","#D55E00", "#CC79A7", "#56B4E9")
 #col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#AA4B56")
col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#D97986")

names(col2) <- celltype

celltype <- unique(uniondegs$celltype)
variables <- unique(uniondegs$variable)
expanded_annotations <- rep(variables, each = length(celltype))
expanded_var_name <- paste0(uniondegs$variable, " - ", uniondegs$celltype)

expanded_var_name <- unique(expanded_var_name)

# Create row and column annotations based on the expanded variable names
row_annotation <- data.frame(variable = expanded_annotations)

rownames(row_annotation) <- expanded_var_name
row_annotation$CellType <- paste0(gsub(".* - ", "", rownames(row_annotation)))

var_color <- c("#E34234", "#3A84D8")#, "#FFB6C1")#, "#E34234", "#ADD8E6",) 
#var_color <- c("#FFA500", "#E34234", "#29465B")#, "#FFB6C1")#, "#E34234", "#ADD8E6",) 

#cell_colors <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#AA4B56")
cell_colors <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#D97986")


names(var_color) <- variable_names
names(cell_colors) <- celltype

annotation_colors <- list(variable  = var_color, CellType = cell_colors)

#figfn <- paste0("02.1.heatmap_zscore_unionDEGs_2Var_6celltypes_annotated_colorupd.png")
figfn <- paste0("01.3.heatmap_zscore_unionDEGs_2Var_5celltypes_annotated_celltype_ordered.png")

png(figfn, width = 3000, height = 1300, res = 240)

# Define a custom color palette with a more gradual transition near zero
breaks <- seq(-4, 4, length.out = 101) # Adjust range based on your data
colors <- colorRampPalette(c("blue", "white", "red"))(100)

heatmap_matrix_ordered <- heatmap_matrix[row_order, colnames(heatmap_matrix)]

# Plot heatmap
p1 <- pheatmap(heatmap_matrix_ordered,
         color = colors, 
         breaks = breaks, # Ensure better transition near zero
         cluster_rows = FALSE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE,
         show_colnames = FALSE, # Can be set to TRUE if gene names are not too many
         fontsize = 14, fontsize_number = 5.5,
         #scale = "row",  # Normalize within rows for better visualization
         annotation_row = row_annotation,  # Add row annotations
         #annotation_col = col_annotation,  # Add column annotations
         annotation_colors = annotation_colors # Add colors for annotations
         )  

  print(p1)
  dev.off()




########################################### for logFC/zscore

# Convert to wide format: Rows = variable_celltype, Columns = Genes
heatmap_data <- uniondegs %>%
  #dplyr::select(variable_celltype, identifier, logFC) %>% # Use logFC or another numeric value
  dplyr::select(variable_celltype, identifier, zscore) %>% # Use logFC or another numeric value
  #pivot_wider(names_from = identifier, values_from = logFC, values_fill = 0) %>%
  pivot_wider(names_from = identifier, values_from = zscore, values_fill = 0) %>%
  column_to_rownames("variable_celltype")

# Convert to matrix for heatmap plotting
heatmap_matrix <- as.matrix(heatmap_data)
rownames(heatmap_matrix)


row_order <- c( #"Cytokines - CD4+ T cells 1", "Cytokines - CD4+ T cells 2", "Cytokines - CD8+ T cells", "Cytokines - NK cells", "Cytokines - Monocytes", "Cytokines - B cells",
              "Psychological Stress - R0 T CD4+", "Psychological Stress - R1 T CD8+", "Psychological Stress - R2 NK", "Psychological Stress - R3 Monocyte", "Psychological Stress - R4 B", #"Psychological Stress - R6 Dendritic cell",
              "Social Support - R0 T CD4+", "Social Support - R1 T CD8+", "Social Support - R2 NK", "Social Support - R3 Monocyte", "Social Support - R4 B")#, "Social Support - R6 Dendritic cell")#, 

heatmap_matrix_ordered <- heatmap_matrix[row_order, colnames(heatmap_matrix)]

# Plot the clustered heatmap
#figfn <- paste0("heatmap_logFC_unionDEGs_4Var_7celltypes.png")
#figfn <- paste0("02.1.heatmap_logFC_unionDEGs_4Var_6celltypes3_annotations_all_genes.png")
#figfn <- paste0("02.2.heatmap_logFC_unionDEGs_2Var_6celltypes_annotations_all_genes_scalerowoff.png")
#figfn <- paste0(outFolder, "01.4.heatmap_logFC_unionDEGs_2Var_5celltypes_annotations_all_genes_scalerowoff_celltypeorder.png")
figfn <- paste0(outFolder, "01.3.heatmap_zscore_unionDEGs_2Var_5celltypes_annotations_all_genes_scalerowoff_celltypeorder.png")
png(figfn, width = 3000, height = 1300, res = 240)

# Define a custom color palette with a more gradual transition near zero
breaks <- seq(-1, 1, length.out = 101) # Adjust range based on your data
colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Plot heatmap
p1 <- pheatmap(heatmap_matrix_ordered,
         color = colors, 
         breaks = breaks, # Ensure better transition near zero
         cluster_rows = FALSE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE,
         show_colnames = FALSE, # Can be set to TRUE if gene names are not too many
         fontsize = 14, fontsize_number = 5.5,
         #scale = "row",  # Normalize within rows for better visualization
         annotation_row = row_annotation,  # Add row annotations
         #annotation_col = col_annotation,  # Add column annotations
         annotation_colors = annotation_colors # Add colors for annotations
         )  

  print(p1)
  dev.off()

#figfn <- paste0(outFolder, "01.4.heatmap_logFC_unionDEGs_2Var_5celltypes_annotations_all_genes_scalerowoff_celltypeorder.pdf")
figfn <- paste0(outFolder, "01.3.heatmap_zscore_unionDEGs_2Var_5celltypes_annotations_all_genes_scalerowoff_celltypeorder.pdf")
pdf(figfn, width=14, height=7)
print(p1)
dev.off()
