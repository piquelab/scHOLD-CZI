####
####
library(reactome.db)
library(org.Hs.eg.db)
library(dplyr)

library(clusterProfiler)

###
###
rm(list=ls())


#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/up/")
#dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

#platePrefix <- "ALL.0.15.13." 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/GO/enrichPath/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/pathway_scores/")
outdir="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/pathway_scores/"


### save all pathway file
fname = paste0(getwd(), "/reactomePA_all_geneIDs.txt")
dfpath <- read.table(fname, sep="\t", header=T, quote='"', comment="")

#colnames(dfpath)
#head(dfpath)[1:3,1:2]

# load normalized gene expression counts: 
normdir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/norm_counts/"
#normdir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data/"
#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data
#C0_CTRL_rna.normal.rds

#fname <- paste0(normdir, "C0_CTRL_rna.normal.rds")
#counts_c0 <- readRDS(fname)
#head(counts_c0)[1:5,1:5]

PATHID <- dfpath$PATHID
PATHNAME <- dfpath$PATHNAME

clusters <- c("C0", "C1", "C2", "C3", "C4")#, "C5", "C6")
#clusters <- c("C2", "C3", "C4", "C5", "C6")


#load each cluster count 

for(c in 1:length(clusters)){

	clust <- clusters[c]
	fname <- paste0(normdir, clust, "_CTRL_rna.normal.rds")
	counts <- readRDS(fname)


# Initialize output dataframe
pathway_scores <- data.frame(
  PATHID = dfpath$PATHID,
  row.names = dfpath$PATHID
)

pathway_genes <- data.frame(PATHID = character(), n_genes = numeric(), genes_used = character(), stringsAsFactors = FALSE)

# Loop through each pathway
for (i in 1:nrow(dfpath)) {
  pathway_id <- dfpath$PATHID[i]
  genes_entrez <- unlist(strsplit(dfpath$ENTREZID[i], ","))

  # Convert ENTREZID to ENSEMBL
  gene_map <- try(bitr(genes_entrez, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Hs.eg.db), silent=TRUE)
  
  if (inherits(gene_map, "try-error") || nrow(gene_map) == 0) {
    next  # Skip if conversion failed
  }

  gene_ensembl <- unique(gene_map$SYMBOL)

  # Subset counts matrix for those genes
  matched_genes <- intersect(gene_ensembl, rownames(counts))
  if (length(matched_genes) == 0) {
    next
  }

  # Calculate column-wise sum (one score per individual)
  pathway_score <- colSums(counts[matched_genes, , drop = FALSE])

  # Add scores as a new row
  pathway_scores[pathway_id, names(pathway_score)] <- pathway_score

  # Add pathway gene summary
  pathway_genes <- rbind(
    pathway_genes,
    data.frame(
      PATHID = pathway_id,
      n_genes = length(matched_genes),
      genes_used = paste(matched_genes, collapse = ","),
      stringsAsFactors = FALSE
    )
  )
}

# Optional: convert rownames back to column if needed
pathway_scores$PATHID <- rownames(pathway_scores)

### save all pathway file
fname = paste0(outdir, clust, "_reactome_pathway_scores_CTRL_HOLD_normalized_counts.txt")
write.table(pathway_scores, fname, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

### save patwhay genes
fname = paste0(outdir, "gene_counts_", clust, "_reactome_pathway_scores_CTRL_HOLD.txt")
write.table(pathway_genes, fname, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}




