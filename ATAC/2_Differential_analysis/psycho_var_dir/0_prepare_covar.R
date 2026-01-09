####
####
library(tidyverse)
library(Matrix)
library(data.table)


##
## library(cowplot)
## library(RColorBrewer)
## library(scales)
## library(viridis)
## library(circlize)
## library(ComplexHeatmap)
## library(openxlsx)
## library(ggrastr)

###

rm(list=ls())


###
### combine covariates, batch and PCs into a single file


node_grid <- "/rs/rs_grp_scatac/"
path_wk <- paste(node_grid, "schold/ATAC/sc-atac-cziHOLD/", sep="")


###
### covariates files name
  
infn_cv <- paste(path_wk, "HOLD-CZI_covariates_Ali_updated_20250528.txt", sep="")
infn_batch <- paste(path_wk, "HOLD-CZI_dbgapIDs_batch_20240917.txt", sep="")
infn_pc <- paste(path_wk, "HOLD-CZI_geno_pc_nokin.txt", sep="")

###
### psychosocial factors

cv <- read.table(infn_cv, header=T)
##
## col2 <- names(cv)
## var2 <- sort(c(col2[grepl("^z_", col2)], "cytocomp"))
## write.table(var2, file="cytokines_var11.txt", sep="", quote=F, row.names=F, col.names=F)
## var2 <- col2[grepl("PF", col2)]
## var2 <- c("PFOA", "PFOA_LOD", "Log_PFOA_1", "PFHxS", "PFHxS_LOD", "Log_PFHxS_1", "PFNA", "PFNA_LOD", "Log_PFNA_1",
##      "LPFHpS", "LPFHpS_LOD", "Log_LPFHpS_1", "PFDA", "PFDA_LOD",
##      "Log_PFDA_1", "PFOS", "PFOS_LOD", "Log_PFOS_1", "PFUdA", "PFUdA_LOD", "Log_PFUdA_1")
PF_comb <- c("PFOA", "PFHxS", "PFNA", "LPFHpS", "PFDA", "PFOS", "PFUdA")
cv <- cv%>%mutate(sumPFAS=cv%>%dplyr::select(all_of(PF_comb))%>%rowSums(na.rm=T),
                  log_sumPFAS=log2(sumPFAS+1))
## write.table(sort(c(var2, "sumPFAS", "log_sumPFAS")), file="PFAS_vars.txt", sep="", quote=F, row.names=F, col.names=F)


###
### variables
varSel0 <- read.table("psychosocial_top10.txt", header=F)$V1
varSel2 <- read.table("cytokines_var11.txt", header=F)$V1
varSel3 <- read.table("PFAS_vars.txt")$V1

varSel <- c(varSel0, varSel2, varSel3)

cv2 <- cv%>%dplyr::select(sampleID=dbgap.ID, all_of(c("age", "sex_alph", varSel)))


###
### get only batch 
batch_df <- read.table(infn_batch, header=T)%>%dplyr::select(sampleID=dbgap.ID, Batch)

### get only PCs infor
pcs <- read.table(infn_pc, header=T)%>%dplyr::select(sampleID=Sample_ID, PC1, PC2, PC3)

## combine
cv_full <- cv2%>%inner_join(batch_df, by="sampleID")%>%inner_join(pcs, by="sampleID")
cv_full <- cv_full%>%dplyr::select(all_of(c("sampleID", "Batch", "sex_alph", "PC1", "PC2", "PC3", "age", varSel)))
  
write.table(cv_full, file="0.1_covariates_varSel.txt", quote=F, sep="\t", row.names=F)


###
### 


###
### END

## ###
## ### matrix for covariates
## x1 <- model.matrix(~0+cv$sex_alph)
## colnames(x1) <- gsub("cv\\$", "", colnames(x1))

## ##
## x2 <- model.matrix(~0+cv$Batch)
## colnames(x2) <- gsub("cv\\$", "", colnames(x2))

## cv_new <- cbind(x1, x2, cv[,c("age", "PC1", "PC2", "PC3", "sampleID")]) 

## ###
## #### save txt file
## write.table(cv_new, file="0.2_covariates_design.txt", quote=F, sep="\t", row.names=F)





