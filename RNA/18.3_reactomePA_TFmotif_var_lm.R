####
####
library(reactome.db)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

library(clusterProfiler)

###
###
rm(list=ls())



######################################################################## get table of results

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/GO/enrichPath/for_motifs/")

#atacdir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/Example_show/"
atacdir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/"

#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz
#outdir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/GO/enrichPath/for_motifs/"

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/ATAC/")
outdir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/ATAC/pathway/"

pathwaydir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/pathway_scores/"

# INPUTS
cluster_rna <- c("C0")
cluster_atac <- c("C0")
celltypes <- c("CD4+ T cell")
variable <-  c("PSS_all_mean", "ISEL_Mean")
variable_name <-c("Psychological Stress", "Social Support")

# Read list of motifs with negative zval_diff
fname = "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/ATAC/00_C0_DAMs_opposite_sign.txt"
motifs_df <- read.table(fname, sep="\t", header=TRUE, quote='"', comment="")
motifs <- motifs_df %>% filter(sign(estimate_ISEL_Mean) == -1) %>% pull(gene)


# Output list to collect results
all_results <- list()

for(i in seq_along(cluster_rna)) {
  cl_rna <- cluster_rna[i]
  cl_atac <- cluster_atac[i]
  celltype <- celltypes[i]

#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/psycho_var_dir/0.1_covariates_varSel.txt
  # Load motif list
  #res_motif <- read.table(paste0(atacdir, "2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz"), header=TRUE)
  res_motif <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz", header=T)
  #motif2 <- sort(unique(res_motif$motif_name))
  motif2 <- sort(unique(res_motif$gene))

  for(j in seq_along(variable)) {
    var0 <- variable[j]
    var00 <- variable_name[j]

    # Load covariates
    cv0 <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/psycho_var_dir/0.1_covariates_varSel.txt", header=TRUE, sep="\t")
    colSel <- c("sampleID", "Batch", "sex_alph", "PC1", "PC2", "age", var0)
    cv <- cv0 %>% dplyr::select(all_of(colSel)) %>% drop_na(all_of(colSel))

    # Load pathway scores
    #pathwaydir <- "PATH/TO/YOUR/PATHWAYDIR/"  # Replace with actual path
    fn_score <- paste0(pathwaydir, cl_rna, "_reactome_pathway_scores_CTRL_HOLD_normalized_counts.txt")
    pathwyscore <- read.table(fn_score, sep="\t", header=TRUE)

    # Load gene counts
    fn_count <- paste0(pathwaydir, "gene_counts_", cl_rna, "_reactome_pathway_scores_CTRL_HOLD.txt")
    pathway_genes <- read.table(fn_count, sep="\t", header=TRUE)

    # Load gene-to-pathway annotation
    fname_annot <- paste0(pathwaydir, "reactomePA_all_geneIDs.txt")
    pathwayname <- read.table(fname_annot, sep="\t", header=TRUE)

    for (motif0 in motifs) {
      motif00 <- unlist(str_split(gsub("_.*", "", motif0), "\\.\\."))

      # Get pathway IDs linked to the motif
      PATHIDs <- unique(pathwayname$PATHID[grepl(paste(motif00, collapse="|"), pathwayname$SYMBOL)])
      PATHnames <- unique(pathwayname$PATHNAME[grepl(paste(motif00, collapse="|"), pathwayname$SYMBOL)])
      PATHnames <- gsub("^Homo sapiens: ", "", PATHnames)

      for (p in seq_along(PATHIDs)) {
        pID <- PATHIDs[p]
        pname <- PATHnames[p]

        # Reshape pathway scores
        yscore_df <- pathwyscore %>%
          filter(PATHID == pID) %>%
          pivot_longer(
            cols = -PATHID,
            names_to = "bti",
            values_to = "pathway_score"
          ) %>%
          dplyr::select(bti, pathway_score) %>%
          mutate(
            sampleID = str_extract(bti, "HO\\.\\d+") %>%
              str_replace("\\.", "-")
          )

        # Merge with covariates
        cv2 <- merge(yscore_df, cv, by = "sampleID")

        #if (nrow(cv2) < 10) next  # Skip if not enough data

        # Regress out covariates
        formu <- as.formula(paste("pathway_score ~ factor(Batch) + factor(sex_alph) + PC1 + PC2 + age"))
        lm_cov <- lm(formu, data = cv2)

        # Use residuals to model variable of interest
        plotDF <- cv2 %>% dplyr::select(all_of(c("sampleID", var0))) %>%
          mutate(y = residuals(lm_cov))
        names(plotDF)[2] <- "x"

        lm_final <- lm(y ~ x, data = plotDF)
        coef_x <- summary(lm_final)$coefficients["x", ]
        b0 <- round(as.numeric(coef_x[1]), 4)
        p0 <- as.numeric(coef_x[4])

        # Save result
        all_results[[length(all_results) + 1]] <- data.frame(
          motif_name = motif0,
          variable = var0,
          variable_name = var00,
          PATHID = pID,
          PATHNAME = pname,
          beta = b0,
          pvalue = p0
        )
      }
    }
  }
}

# Combine all rows into one data frame
final_df <- do.call(rbind, all_results)

# View or write output
head(final_df)
fname = paste0(outdir, "TFmotif_reactome_pathway_score_lm_results.txt")
write.table(final_df, fname, sep="\t", row.names=FALSE, quote=FALSE)



# Ensure 'final_df' is available from previous step

# Filter for only the two variables of interest and p < 0.05
df_filtered <- final_df %>%
  filter(variable %in% c("ISEL_Mean", "PSS_all_mean") & pvalue < 0.05)

# Spread to wide format to compare betas side-by-side
df_wide <- df_filtered %>%
  dplyr::select(motif_name, PATHID, PATHNAME, variable, beta, pvalue) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(beta, pvalue),
    names_sep = "_"
  )

# Keep only rows where beta signs are opposite
df_opposite <- df_wide %>%
  filter(
    sign(beta_ISEL_Mean) != sign(beta_PSS_all_mean)
  ) %>%
  dplyr::select(
    motif_name, PATHID, PATHNAME,
    beta_ISEL_Mean, pvalue_ISEL_Mean,
    beta_PSS_all_mean, pvalue_PSS_all_mean
  ) %>% as.data.frame()

# View the result
print(df_opposite)
length(unique(df_opposite$PATHNAME))

fname = paste0(outdir, "TFmotif_reactome_pathway_score_lm_opposite_direction_pathways.txt")
write.table(df_opposite, fname, sep="\t", row.names=FALSE, quote=FALSE)

df_opposite <- df_opposite %>% arrange(pvalue_ISEL_Mean, pvalue_PSS_all_mean)

df_opposite %>% dplyr::select(motif_name, PATHNAME) 















df_opposite <- df_opposite %>% arrange(pvalue)



# order by pvalue
final_df <- final_df %>% arrange(pvalue)

final_df %>% dplyr::select(motif_name, PATHNAME) %>% head(20)

# top 3: 
1        IRF4_308           Signaling by Interleukins
2  STAT1..STAT2_36           Signaling by Interleukins
3  STAT1..STAT2_36     Interleukin-12 family signaling





outdir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/ATAC/pathway/"

#outdir="/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/GO/enrichPath/for_motifs/"
fname=paste0(outdir, "/TFmotif_reactome_pathway_score_lm_results.txt")
paths <- read.table(fname, sep="\t", header=TRUE, quote='"', comment="")

length(unique(paths$motif_name)) #9
length(unique(paths$PATHNAME)) #47 (122 total)

# rank by pvalue: 

#df_opposite <- df_opposite %>% arrange(pvalue_ISEL_Mean, pvalue_PSS_all_mean)
paths <- paths %>% arrange(pvalue)
paths %>% dplyr::select(motif_name, variable, PATHNAME, pvalue, beta) %>% head(20)


# Filter for only the two variables of interest and p < 0.05
df_filtered <- paths %>%
  filter(variable %in% c("ISEL_Mean", "PSS_all_mean") & pvalue < 0.05)

# Spread to wide format to compare betas side-by-side
df_wide <- df_filtered %>%
  dplyr::select(motif_name, PATHID, PATHNAME, variable, beta, pvalue) %>%
  pivot_wider(
    names_from = variable,
    values_from = c(beta, pvalue),
    names_sep = "_"
  )

df_wide <- df_wide %>% arrange(pvalue_ISEL_Mean, pvalue_PSS_all_mean) %>% as.data.frame()

df_wide %>% dplyr::select(PATHNAME, motif_name, pvalue_PSS_all_mean, pvalue_ISEL_Mean, beta_PSS_all_mean, beta_ISEL_Mean) %>% head(20)

nona_df_wide <- df_wide %>% drop_na()


length(unique(nona_df_wide$motif_name)) #9
length(nona_df_wide$motif_name) #31

length(unique(df_wide$motif_name)) #9
length(df_wide$motif_name) #62

fname = paste0(outdir, "02_TFmotif_reactome_pathway_score_lm_pvalue_ranked_na_removed_n31_uniq9.txt")
write.table(nona_df_wide, fname, sep="\t", row.names=FALSE, quote=FALSE)


fname = paste0(outdir, "01_TFmotif_reactome_pathway_score_lm_pvalue_ranked_withNAs_n62_uniq9.txt")
write.table(df_wide, fname, sep="\t", row.names=FALSE, quote=FALSE)


