library(dplyr)
library(tidyr)

rm(list=ls())


############## ISEL and PSS DEGs
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/tables/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

#platePrefix <- "ALL.0.15.13." 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

resdegs <- resall #%>% filter(padj < 0.1)

psdegs <- resdegs %>% filter(var %in% c("ISEL_Mean", "PSS_all_mean"))
psdegs <- psdegs %>% select(!treats)

#fname <- paste0("PSS_ISEL_DEGs_FDR10.txt")
fname <- paste0("PSS_ISEL_DESeq_results.txt")
write.table(psdegs, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

cytall <- resall %>% filter(var == "cytocomp")
fname <- paste0("cytocomp_all_DESeq_results.txt")
write.table(cytall, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


############## ISEL and PSS DAMs
rm(list=ls())

#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
#fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"#

fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
dams <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
sigdams <- dams #%>% filter(padj_t <0.1)

clusters <- c("C0", "C1", "C2", "C3", "C4", "C5")#, "C6")
vars <- c("ISEL_Mean", "PSS_all_mean")

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
dams <- dams %>% filter(psycho_variable %in% vars) #%>% filter(padj_t < 0.1) #%>% filter(!Cluster == "A9")
#dams <- dams %>% filter(padj_t < 0.1)
dim(dams) #189
length(unique(dams$gene)) #158

damsig <- dams %>% select(!c(treat, nind))

dams <- dams %>% filter(!Cluster %in% c("C10","C9", "C8", "C3","C5", "C7")) #%>% filter(padj_t < 0.1) #%>% filter(!Cluster == "A9")
#dams <- dams %>% filter(!Cluster %in% c("C3","C5", "C7")) #%>% filter(padj_t < 0.1) #%>% filter(!Cluster == "A9")

table(dams$psycho_variable, dams$Cluster)


#fname <- paste0("PSS_ISEL_DAMs_FDR10.txt")
fname <- paste0("PSS_ISEL_Motif_activity_results.txt")
write.table(damsig, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

############## LPS DEGs
rm(list=ls())

# load LPS degs
fname=paste0("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/treatment_noCOMBAT/deseqres/ALL.0.1.13.deseqres.RNA-LPS_vs_RNA-CTRL.treatment_noCOMBAT.txt")
treats <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")

treats <- treats #%>% filter(padj < 0.1)

treats$celltype <- "NA"
treats$celltype[treats$cluster == "C0"] <- "R0 T CD4+"
treats$celltype[treats$cluster == "C1"] <- "R1 T CD8+"
treats$celltype[treats$cluster == "C2"] <- "R2 NK"
treats$celltype[treats$cluster == "C3"] <- "R3 Monocyte"
treats$celltype[treats$cluster == "C4"] <- "R4 B"
treats$celltype[treats$cluster == "C5"] <- "R5 DC"

treats <- treats %>% select(!contrast)
treats <- treats %>% filter(!cluster %in% c("C5")) #%>% filter(padj_t < 0.1) #%>% filter(!Cluster == "A9")

#fname <- paste0("LPS_DEGs_FDR10.txt")
fname <- paste0("LPS_DESeq_results.txt")
write.table(treats, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



############## LPS DAMs
rm(list=ls())

# load LPS degs
# treatment dams
#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_summary/2_response_DAMs.comb.txt.gz"
fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_treat/2_th20_plotData.comb.txt.gz"

treats <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

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

#treats <- treats %>% filter(!celltype=="A8 DC")

treats <- treats %>% filter(padj < 0.1)
treats <- treats %>% select(!nind)


treats <- treats %>% filter(!Cluster %in% c("C10","C9", "C8", "C3","C5", "C7")) #%>% filter(padj_t < 0.1) #%>% filter(!Cluster == "A9")

table(treats$Cluster)

#fname <- paste0("LPS_DAMs_FDR10.txt")
fname <- paste0("LPS_Motif_activity_results.txt")
write.table(treats, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)




########################################## GO results
rm(list=ls())

#/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/up/
#fname = "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/up/enriched_pathways_FDR10_allcelltypes_2var_ISEL_updated.txt"
fname = "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/updown/plots/enriched_pathways_FDR10_allcelltypes_2var_ISEL_updated.txt"

goup <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

goup <- goup %>% select(!c(treat, Variable))

fname <- paste0("Reactome_PSSup_ISELdown_FDR10.txt")
write.table(goup, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



######################################## WHR corrected 
rm(list=ls())

dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"
#platePrefix <- "ALL.0.15.13." 
#plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
#plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_WHRcov" 

fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_WHRcov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(!cluster == "C5") #%>% filter(padj < 0.1)

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
#vars <- c("PSS_all_mean")#"ISEL_Mean")#, "PSS_all_mean")
vars <- c("ISEL_Mean", "PSS_all_mean")

resdegs2$celltype <- "NA"
resdegs2$celltype[resdegs2$cluster == "C0"] <- "R0 T CD4+"
resdegs2$celltype[resdegs2$cluster == "C1"] <- "R1 T CD8+"
resdegs2$celltype[resdegs2$cluster == "C2"] <- "R2 NK"
resdegs2$celltype[resdegs2$cluster == "C3"] <- "R3 Monocyte"
resdegs2$celltype[resdegs2$cluster == "C4"] <- "R4 B"
#resdegs2$celltype[resdegs2$cluster == "C5"] <- "R6 DC"

resdegs2 <- resdegs2 %>% filter(var %in% vars)

#resdegs2 <- resdegs2 %>% select(!treats)

fname <- paste0("PSS_ISEL_WHR-corrected_all_DESeq_results.txt")
write.table(resdegs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

######################################## SES corrected
rm(list=ls())
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"
#platePrefix <- "ALL.0.15.13." 
#plateSufix1 <- "-RNA-CTRL.noCombat_DESeq" 
#plateSufix2 <- "-RNA-CTRL.noCombat_DESeq_SEScov" 

fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_SEScov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(!cluster == "C5") #%>% filter(padj < 0.1)


clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("ISEL_Mean", "PSS_all_mean")
#vars <- c("PSS_all_mean")

resdegs2$celltype <- "NA"
resdegs2$celltype[resdegs2$cluster == "C0"] <- "R0 T CD4+"
resdegs2$celltype[resdegs2$cluster == "C1"] <- "R1 T CD8+"
resdegs2$celltype[resdegs2$cluster == "C2"] <- "R2 NK"
resdegs2$celltype[resdegs2$cluster == "C3"] <- "R3 Monocyte"
resdegs2$celltype[resdegs2$cluster == "C4"] <- "R4 B"
#resdegs2$celltype[resdegs2$cluster == "C5"] <- "R6 DC"

resdegs2 <- resdegs2 %>% filter(var %in% vars)

#resdegs2 <- resdegs2 %>% select(!treats)

#fname <- paste0("PSS_ISEL_SES-corrected_DEGs_FDR10.txt")
fname <- paste0("PSS_ISEL_SES-corrected_all_DESeq_resluts.txt")
write.table(resdegs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

######################################## ISEL corrected

rm(list=ls())
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_iselcov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(!cluster == "C5") #%>% filter(padj < 0.1)

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("PSS_all_mean")#"ISEL_Mean")#, "PSS_all_mean")

resdegs2$celltype <- "NA"
resdegs2$celltype[resdegs2$cluster == "C0"] <- "R0 T CD4+"
resdegs2$celltype[resdegs2$cluster == "C1"] <- "R1 T CD8+"
resdegs2$celltype[resdegs2$cluster == "C2"] <- "R2 NK"
resdegs2$celltype[resdegs2$cluster == "C3"] <- "R3 Monocyte"
resdegs2$celltype[resdegs2$cluster == "C4"] <- "R4 B"
#resdegs2$celltype[resdegs2$cluster == "C5"] <- "R6 DC"

resdegs2 <- resdegs2 %>% filter(var %in% vars)
#resdegs2 <- resdegs2 %>% select(!treats)


#fname <- paste0("PSS_ISEL-corrected_DEGs_FDR10.txt")
fname <- paste0("PSS_ISEL-corrected_all_DESeq_results.txt")
write.table(resdegs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

######################################## PSS corrected 

rm(list=ls())
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

fname =  paste(dataDir, "allres_PSS_ISEL_noCombat_DESeq_PSScov.txt", sep = "")
resall2 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
resdegs2 <- resall2 %>% filter(!cluster == "C5") # %>% filter(padj < 0.1)

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("ISEL_Mean")

resdegs2$celltype <- "NA"
resdegs2$celltype[resdegs2$cluster == "C0"] <- "R0 T CD4+"
resdegs2$celltype[resdegs2$cluster == "C1"] <- "R1 T CD8+"
resdegs2$celltype[resdegs2$cluster == "C2"] <- "R2 NK"
resdegs2$celltype[resdegs2$cluster == "C3"] <- "R3 Monocyte"
resdegs2$celltype[resdegs2$cluster == "C4"] <- "R4 B"
#resdegs2$celltype[resdegs2$cluster == "C5"] <- "R6 DC"

resdegs2 <- resdegs2 %>% filter(var %in% vars)
#resdegs2 <- resdegs2 %>% select(!treats)



#fname <- paste0("ISEL_PSS-corrected_DEGs_FDR10.txt")
fname <- paste0("ISEL_PSS-corrected_all_DESeq_results.txt")
write.table(resdegs2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



################################# DARs: 


############## LPS DAMs
rm(list=ls())

# load LPS degs
# treatment dams
#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/treat_Cluster_res0.07_th20_summary/2_response_DAMs.comb.txt.gz"
fname="/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/1_DiffPeak.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
dars <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')


fname="/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/1_DiffPeak.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
res_data <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")

dars <- res_data

dars$celltype <- "NA"
dars$celltype[dars$Cluster == "C0"] <- "A0 T CD4+"
dars$celltype[dars$Cluster == "C1"] <- "A1 T CD8+"
dars$celltype[dars$Cluster == "C2"] <- "A2 NK"
dars$celltype[dars$Cluster == "C3"] <- "A3 T CD4+"
dars$celltype[dars$Cluster == "C4"] <- "A4 Monocyte"
dars$celltype[dars$Cluster == "C5"] <- "A5 T CD4+"
dars$celltype[dars$Cluster == "C6"] <- "A6 B"
dars$celltype[dars$Cluster == "C7"] <- "A7 T CD4+"
dars$celltype[dars$Cluster == "C8"] <- "A8 T CD4+"
dars$celltype[dars$Cluster == "C9"] <- "A9 T CD4+"
dars$celltype[dars$Cluster == "C10"] <- "A10 DC"

vars <- c("ISEL_Mean", "PSS_all_mean")
dars <- dars %>% filter(psycho_variable %in% vars)


dars <- dars %>% filter(!Cluster %in% c("C10","C9", "C8", "C3","C5", "C7")) #%>% filter(padj_t < 0.1) #%>% filter(!Cluster == "A9")


table(dars$psycho_variable, dars$Cluster)

#treats <- treats %>% filter(!celltype=="A8 DC")

darss <- dars %>% filter(p.adjusted < 0.1)
table(darss$psycho_variable, darss$Cluster)
#treats <- treats %>% select(!nind)

#fname <- paste0("LPS_DAMs_FDR10.txt")
fname <- paste0("DARs_DESeq_results.txt")
write.table(dars, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)




###################### CellRanger RNA library stats

fname = "/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/counts_cellranger_2024-04-19/all_summary.csv"
res_data <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")

#columns to select: 
# colnames(res_data)
# "Library.ID"                     "Estimated.Number.of.Cells"
# "Mean.Reads.per.Cell"            "Median.Genes.per.Cell"
#"Number.of.Reads"                "Valid.Barcodes"
# "Sequencing.Saturation"          "Reads.Mapped.to.Genome"
# "Reads.Mapped.Antisense.to.Gene" "Fraction.Reads.in.Cells"
# "Total.Genes.Detected"           "Median.UMI.Counts.per.Cell"




###################### CellRanger ATAC library stats






