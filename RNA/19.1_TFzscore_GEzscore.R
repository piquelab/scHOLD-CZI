library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

rm(list=ls())

setwd("/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/ATAC/TFmotif/")

# load data
sigmotifs <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/888_pub.outs/plots_main/TableS_175motif_sigs.infor.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="\"")
archrmotifs <- read.table("/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/1.2_ArchR_process/4_motif_output/Motif_jaspar2022.motifinfor.ArchR.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="\"")
resall <- read.table("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/allres_allvars_noCombat_DESeq.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="\"")

# prep motifs
sigmotifs <- sigmotifs %>%
  mutate(motif_clean = sub("_.*", "", motif_name),
         var_clust = paste0(psycho_variable, "_", cluster_rna, "_", cluster_atac))

# add symbol info
sigmotifs <- left_join(sigmotifs, archrmotifs, by="motif_name")

# split complex motifs
motifs_expanded <- sigmotifs %>%
  mutate(TF_raw = name) %>%
  separate_rows(TF_raw, sep = "\\.{1,2}|::|-") %>%
  mutate(
    TF_raw = str_trim(TF_raw),
    tf_var_rna = paste0(TF_raw, "_", psycho_variable, "_", cluster_rna)
  ) %>%
  filter(!is.na(TF_raw) & TF_raw != "")

# prepare DEseq
resall <- resall %>%
  mutate(tf_var_rna = paste0(identifier, "_", var, "_", cluster))

# merge TF motif activity + TF expression
merged <- left_join(motifs_expanded, resall, by="tf_var_rna", suffix=c("_ATAC","_RNA"))
merged <- as.data.frame(merged)

# correct RNA→ATAC cluster mapping
clmap <- data.frame(
  cluster_rna = paste0("C",0:4),
  cluster_atac_corrected = c("C0","C1","C2","C4","C6")
)

merged <- merged %>%
  select(-cluster_atac) %>%
  left_join(clmap, by="cluster_rna") %>%
  mutate(var_clust_corrected = paste0(psycho_variable, "_", cluster_rna, "_", cluster_atac_corrected)) %>%as.data.frame()

# save
write.table(merged, "final_TF_motif_expression_merged.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


sigge <- merged %>% filter(padj<0.1)


df <- merged

# do it for enrichments at 5% FDR DAMs
df <- merged %>% filter(padj_diff<0.05)

var_color <- c("PSS_all_mean" = "#E34234",
               "ISEL_Mean"    = "#3A84D8")

# clean cell type: remove any "C0 ", "C1 ", or "CX" prefix
df <- df %>%
  mutate(CellType = str_replace(celltype, "^C[0-9X]+\\s*", ""))

df$CellType <- "NA"
df$CellType[df$cluster_rna=="C0"] <- "R0 A0 CD4+ T"
df$CellType[df$cluster_rna=="C1"] <- "R1 A1 CD8+ T"
df$CellType[df$cluster_rna=="C2"] <- "R2 A2 NK"
df$CellType[df$cluster_rna=="C3"] <- "R3 A4 Monocyte"

# filter RNA-significant TFs
df2 <- df %>%
  filter(!is.na(padj) & padj < 0.1)

# compute z-scores
df2 <- df2 %>%
  mutate(
    z_TF = beta_diff / std_diff,
    z_DEG = logFC / SE
  )

# remove duplicate TFs from complex motifs
df2 <- df2 %>%
  distinct(TF_raw, psycho_variable, cluster_rna, .keep_all = TRUE)

# plot
p <- ggplot(df2, aes(x = z_DEG, y = z_TF,
                     color = psycho_variable,
                     shape = CellType)) +
  geom_point(size = 3) +
  geom_abline(linetype = 2, color = "grey50") +
  scale_color_manual(values = var_color) +
  theme_minimal() +
  labs(
    x = "RNA TF gene expression z-score",
    y = "ATAC TF motif activity z-score",
    color = "Variable",
    shape = "Cell Type",
    title = "TF Motif Activity vs TF Gene Expression"
  )

#png("TF_zscore_DEG_vs_motif_activity.png", width = 1000, height = 800, res = 150)
png("TF_zscore_DEG_vs_motif_activity_FDR5.png", width = 1000, height = 800, res = 150)
print(p)
dev.off()



oppsign <- df2 %>% filter(z_DEG>2 & z_TF<2) %>% select(motif_name, identifier, var_clust, CellType, padj_diff, beta_diff, z_TF, padj, logFC, z_DEG)
# save
#write.table(oppsign, "TF_DAMzscore_vs_GEzscore_opposite_signs.txt",
write.table(oppsign, "TF_DAMzscore_vs_GEzscore_opposite_signs_FDR5.txt",
            sep="\t", quote=FALSE, row.names=FALSE)


df2select <- df2 %>% select(motif_name, identifier, var_clust, CellType, padj_diff, beta_diff, z_TF, padj, logFC, z_DEG)
# save
#write.table(df2select, "TF_DAMzscore_vs_GEzscore_all_significants.txt",
 write.table(df2select, "TF_DAMzscore_vs_GEzscore_all_significants_FDR5.txt",
           sep="\t", quote=FALSE, row.names=FALSE)



##### just for CD4+ cells: 



# remove duplicate TFs from complex motifs
df3 <- df2 %>%
   filter(cluster_rna == "C0")

# plot
p <- ggplot(df3, aes(x = z_DEG, y = z_TF,
                     color = psycho_variable)) +
  geom_point(size = 3) +
  geom_abline(linetype = 2, color = "grey50") +
  scale_color_manual(values = var_color) +
  theme_minimal() +
  labs(
    x = "RNA TF gene expression z-score",
    y = "ATAC TF motif activity z-score",
    color = "Variable",
    shape = "Cell Type",
    title = "TF Motif Activity vs TF Gene Expression"
  )

#png("TF_zscore_DEG_vs_motif_activity_C0_CD4.png", width = 1000, height = 800, res = 150)
png("TF_zscore_DEG_vs_motif_activity_C0_CD4_FDR5.png", width = 1000, height = 800, res = 150)
print(p)
dev.off()


oppsigndf3 <- df3 %>%
  group_by(TF_raw, CellType, cluster_rna) %>%
  filter(n_distinct(psycho_variable) == 2) %>%
  mutate(sign_TF = sign(z_TF)) %>%
  # check if signs are opposite
  filter(abs(sum(sign_TF)) == 0) %>%
  ungroup() %>% as.data.frame()

# plot
p <- ggplot(oppsigndf3, aes(x = z_DEG, y = z_TF,
                     color = psycho_variable)) +
  geom_point(size = 3) +
  geom_abline(linetype = 2, color = "grey50") +
  scale_color_manual(values = var_color) +
  theme_minimal() +
  labs(
    x = "RNA TF gene expression z-score",
    y = "ATAC TF motif activity z-score",
    color = "Variable",
    shape = "Cell Type",
    title = "TF Motif Activity vs TF Gene Expression"
  )

#png("TF_zscore_DEG_vs_motif_activity_C0_CD4_oppos_ISEL_PSS.png", width = 800, height = 600, res = 150)
png("TF_zscore_DEG_vs_motif_activity_C0_CD4_oppos_ISEL_PSS_FDR5.png", width = 800, height = 600, res = 150)
print(p)
dev.off()

