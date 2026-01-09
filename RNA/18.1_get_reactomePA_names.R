###################################################
library(reactome.db)
library(org.Hs.eg.db)
library(dplyr)


#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/GO/enrichPath/")
#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/GO/enrichPath/")


setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/pathway_scores/")
outdir="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/pathway_scores/"

# Step 1: Get all gene-to-pathway mappings from reactome.db
gene2path <- AnnotationDbi::select(
  reactome.db,
  keys = keys(reactome.db, keytype = "ENTREZID"),
  columns = c("PATHID", "PATHNAME"),
  keytype = "ENTREZID"
)

# Step 2: Add gene symbols from org.Hs.eg.db
gene2path <- gene2path %>%
  left_join(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(gene2path$ENTREZID),
      columns = "SYMBOL",
      keytype = "ENTREZID"
    ),
    by = "ENTREZID"
)

# Step 3 (optional): Remove any duplicates or NA values
gene2path <- gene2path %>%
  distinct() %>%
  filter(!is.na(SYMBOL))

# Step 4: View or save
head(gene2path)
# View(gene2path)  # Optional if using RStudio
# write.csv(gene2path, "all_reactome_pathways_with_genes.csv", row.names = FALSE)


### save gene to pathways file
fname = paste0(getwd(), "/gene_symbol_to_reactomePA_IDs.txt")
write.table(gene2path, fname, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

############################
library(reactome.db)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(stringr)

# Step 1: Get all human Entrez Gene IDs
human_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Step 2: Map only human Entrez IDs to Reactome pathways
path2gene_df <- AnnotationDbi::select(
  reactome.db,
  keys = human_entrez,
  keytype = "ENTREZID",
  columns = c("PATHID", "PATHNAME")
)

# Step 3: Remove any rows with NA (some genes may not have pathway annotation)
path2gene_df <- na.omit(path2gene_df)

# Step 4: Collapse gene IDs per pathway
pathway_gene_df <- path2gene_df %>%
  group_by(PATHID, PATHNAME) %>%
  summarise(
    ENTREZIDS = paste(unique(ENTREZID), collapse = ","),
    n_genes = n_distinct(ENTREZID),
    .groups = "drop"
  )

# Step 5 (Optional): Add gene symbols
library(org.Hs.eg.db)
entrez2symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(path2gene_df$ENTREZID),
  keytype = "ENTREZID",
  columns = "SYMBOL"
)

# Merge to get gene symbols
path2gene_df <- left_join(path2gene_df, entrez2symbol, by = "ENTREZID")

# Collapse symbols per pathway
pathway_gene_df <- path2gene_df %>%
  group_by(PATHID, PATHNAME) %>%
  summarise(
    ENTREZIDS = paste(unique(ENTREZID), collapse = ","),
    SYMBOLS   = paste(unique(SYMBOL), collapse = ","),
    n_genes   = n_distinct(ENTREZID),
    .groups = "drop"
  ) 

# View result
head(pathway_gene_df)

### save all pathway file
fname = paste0(getwd(), "/reactomePA_all_geneIDs.txt")
write.table(pathway_gene_df, fname, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

