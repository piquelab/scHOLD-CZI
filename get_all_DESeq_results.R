library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(ggplot2)
library(data.table)
library(reshape)
library(stringr)
library(data.table)
library(ggrepel)
library(qvalue)


#################################
### get all sig DEGs
###################################
rm(list=ls())

outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/tables/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq" 
#plateSufix <- "-RNA-LPS.noCombat_DESeq" 

#plateSufix <- "-RNA-CTRL.SES_PCs_sex_age_and_treats_generem" 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq_SEScov" 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq_iselcov" 
plateSufix <- "-RNA-CTRL.noCombat_DESeq_WHRcov" 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq_PSScov" 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq_iselcov.noindfilt" 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq_PSScov.noindfilt" 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq.noindfilt" 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq_iselcov" 
#plateSufix <- "-RNA-CTRL.noCombat_DESeq_PSScov" 

#################################################
################# noCombat DESeq
#################################################

# Variables to loop over
variables <- c("PSS_all_mean", "ISEL_Mean")#, "cytocomp")
variable_names <- c("Perceived Stress", "Social Support")#, "Cytokines") 

# Variables to loop over
#variables <- c("PSS_all_mean")
#variable_names <- c("Perceived Stress") 

# Variables to loop over
#variables <- c("ISEL_Mean")#, "cytocomp")
#variable_names <- c("Social Support")#, "Cytokines") 

clusters <- c("0", "1", "2", "3", "4")#, "5")#  , "6")
#clusters <- c("0","4")

#celltype <- c("CD4+ T cells 1", "CD4+ T cells 2", "CD8+ T cells", "NK cells", "Monocytes", "B cells", "Dendritic  Cells")
#celltype <- c("R0 T CD4+", "R1 T CD4+", "R2 T CD8+", "R3 NK", "R4 Monocyte", "R5 B")#, "R6 Dendritic cell")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R5 DC")


###########################
# List of variables and their names
#variables <- c("isel", "pr_comp", "PSS_all_mean")#, "cytocomp")
#variable_names <- c("Social Support", "Psychological Resources", "Perceived Stress")#, "Cytokines")

# Initialize an empty list to store data for all variables
all_data <- list()
allres <- list()

# Loop over variables
for (i in seq_along(variables)) {
  myvar <- variables[i]
  var_name <- variable_names[i]
  
  dfvar <- list()
  
  # Loop over clusters
  for (j in seq_along(clusters)) {
    clust <- clusters[j]
    clustname <- celltype[j]
    
    fname <- paste(dataDir, platePrefix, "C", clust, ".deseqres_", myvar, plateSufix, ".txt", sep = "")
    
    # Read data
    res_data <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")
    res_data$celltype <- clustname
    
    # Store each data frame
    dfvar[[clust]] <- res_data
  }
  
  # Combine results for all clusters
  res<- bind_rows(dfvar)
  ressig <- res %>% filter(padj < 0.1)

  # Store data
  all_data[[i]] <- ressig
    allres[[i]] <- res

}


allsig <- bind_rows(all_data)
resall <- bind_rows(allres)


#fname=paste0(dataDir, "/", "allres_allvars_noCombat_DESeq.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq_SEScov.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq_iselcov.txt")
fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq_WHRcov.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq_PSScov.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq_iselcov.noindfilt.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq_PSScov.noindfilt.txt")
#fname=paste0(dataDir, "/", "allres_allvars_noCombat_DESeq.txt")
#fname=paste0(dataDir, "/", "allres_LPS_allvars_noCombat_DESeq.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq.noindfilt.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq_iselcov.txt")
#fname=paste0(dataDir, "/", "allres_PSS_ISEL_noCombat_DESeq_PSScov.txt")

    write.table(resall, fname, sep = "\t", row.names = FALSE, quote = FALSE)

#allres_PSS_ISEL_20cellfilter_iselcov.txt  nDEGs_PSS_ISEL_20cellfilter_iselcov.txt
#allres_PSS_ISEL_20cellfilter_SEScov.txt   nDEGs_PSS_ISEL_20cellfilter_SEScov.txt
#allres_PSS_ISEL_20cellfilter.txt          nDEGs_PSS_ISEL_20cellfilter_.txt
#allres_PSS_ISEL_20cellfilter_WHRcov.txt   nDEGs_PSS_ISEL_20cellfilter_WHRcov.txt
