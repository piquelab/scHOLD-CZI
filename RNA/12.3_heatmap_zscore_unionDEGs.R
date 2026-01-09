library(scales)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sva)

rm(list=ls())

outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/heatmap/"
#if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13."
plateSufix <- "-RNA-CTRL.noCombat_DESeq"

variables <- c("cytocomp", "PSS_all_mean", "ISEL_Mean")
variable_names <- c("Cytokines", "Psychological Stress", "Social Support")

cluster_number <- c("0","1","2","3","4")
celltypes <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")

col2 <- c("#FF7F00", "#E6E600", "#4DAF4A", "#984EA3", "#D97986")
names(col2) <- celltypes

zscore_list <- list()
gene_list <- list()
siggene_list <- list()

expanded_var <- c()
expanded_var_name <- c()

for (i in seq_along(variables)) {

  myvar <- variables[i]
  var_name <- variable_names[i]

  for (j in seq_along(cluster_number)) {

    cell <- cluster_number[j]
    celltype <- celltypes[j]

    fname <- paste0(
      dataDir, platePrefix, "C", cell,
      ".deseqres_", myvar, plateSufix, ".txt"
    )

    resC <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    resC <- resC %>% dplyr::select(-any_of("baseMean"))

    ressigC <- resC %>% filter(padj < 0.1)

    resC$zscore <- resC$logFC / resC$SE
    resC$symb <- resC$identifier
    ressigC$symb <- ressigC$identifier

    key <- paste0(cell, "_", myvar)

    zscore_list[[key]] <- resC$zscore
    gene_list[[key]] <- resC$symb
    siggene_list[[key]] <- ressigC$symb

    expanded_var <- c(expanded_var, key)
    expanded_var_name <- c(expanded_var_name,
                           paste0(var_name, " - ", celltype))
  }
}

n <- length(expanded_var)

zscore_matrix <- matrix(NA, nrow=n, ncol=n)
p_matrix <- matrix(NA, nrow=n, ncol=n)

rownames(zscore_matrix) <- expanded_var_name
colnames(zscore_matrix) <- expanded_var_name

for (v1 in seq_len(n)) {
  for (v2 in seq_len(n)) {

    union_genes <- union(
      siggene_list[[expanded_var[v1]]],
      siggene_list[[expanded_var[v2]]]
    )

    union_genes_tested <- union_genes[
      union_genes %in% gene_list[[expanded_var[v1]]] &
      union_genes %in% gene_list[[expanded_var[v2]]]
    ]

    z1 <- zscore_list[[expanded_var[v1]]][
      gene_list[[expanded_var[v1]]] %in% union_genes_tested
    ]

    z2 <- zscore_list[[expanded_var[v2]]][
      gene_list[[expanded_var[v2]]] %in% union_genes_tested
    ]

    if (length(z1) < 3 || length(z2) < 3) {
      zscore_matrix[v1, v2] <- 0
      p_matrix[v1, v2] <- 0
      next
    }

    ct <- cor.test(z1, z2, method="pearson")

    zscore_matrix[v1, v2] <- ifelse(ct$p.value <= 0.05, ct$estimate, 0)
    p_matrix[v1, v2] <- ct$p.value
  }
}

write.table(
  zscore_matrix,
  file=paste0(outFolder, "02.6.zscore-correlations_unionDEGs.txt"),
  sep="\t", quote=FALSE, col.names=NA
)

expanded_annotations <- rep(variable_names, each=length(celltypes))

row_annotation <- data.frame(
  Variable = expanded_annotations,
  CellType = gsub(".* - ", "", expanded_var_name)
)

rownames(row_annotation) <- expanded_var_name

var_color <- c("#FFA500", "#E34234", "#3A84D8")
names(var_color) <- variable_names

cell_colors <- col2

annotation_colors <- list(
  Variable = var_color,
  CellType = cell_colors
)

png(
  paste0(outFolder,
         "02.1.zscore-correlation_unionDEGs_heatmap.png"),
  width=2300, height=1500, res=240
)

p2 <- pheatmap(
  zscore_matrix,
  main="z-score correlation (union DEGs)",
  cluster_rows=FALSE,
  cluster_cols=FALSE,
  color=colorRampPalette(c("blue","white","red"))(100),
  breaks=seq(-1,1,length.out=101),
  display_numbers=TRUE,
  fontsize=10,
  fontsize_number=7,
  na_col="white",
  annotation_row=row_annotation,
  annotation_colors=annotation_colors
)

print(p2)
dev.off()

pdf(
  paste0(outFolder,
         "02.1.zscore-correlation_unionDEGs_heatmap.pdf"),
  width=14, height=7
)

print(p2)
dev.off()



png(
  paste0(outFolder,
         "02.3.zscore-correlation_unionDEGs_heatmap_simplified.png"),
  width=2300, height=1500, res=240
)

p2 <- pheatmap(
  zscore_matrix,
  main="z-score correlation (union DEGs)",
  cluster_rows=FALSE,
  cluster_cols=FALSE,
  color=colorRampPalette(c("blue","white","red"))(100),
  breaks=seq(-1,1,length.out=101),
  display_numbers=FALSE,
  fontsize=10,
  fontsize_number=7,
  na_col="white",
  annotation_row=row_annotation,
  annotation_colors=annotation_colors
)

print(p2)
dev.off()
