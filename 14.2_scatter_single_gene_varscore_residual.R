library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)

timestamp()

rm(list=ls())

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/gene/residual/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 


#countSufix <- ".SES_PCs_sex_age_and_treats_adjusted"
deseqSuffix <- ".noCombat_DESeq"
deseqDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/"

#normdata <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data/"
normdata <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/norm_counts/"

#/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/ALL.0.2.11.DESeq_output-RNA-CTRL-PSS_all_meanC0.SES_PCs_sex_age_and_treats_generem.RDS
#load("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/ALL.0.2.11.DESeq_output-RNA-CTRL-PSS_all_meanC0.SES_PCs_sex_age_and_treats_generem.RData")
#variables <- c("pr_comp")
#variable_names <- c("psychological resources composite")


variables <- c("PSS_all_mean", "ISEL_Mean", "cytocomp" )#,"cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
variable_names <- c("Psychological Stress", "Social Support","Cytokines")# "Cytokines") #"socioeconomic status", "diastolic blood pressure", 

#mygene <- "CD52"
#myensemble <- "ENSG00000169442"
#variables <- c("pnsi")
#variable_names <- c("Neighborhood Stress")


########################
### defined function ###
########################

adjGene <- function(cvt, center=T){
   cvt <- cvt%>%mutate(comb=paste(MCls, Batch, sep="_"))
   if(center){
      cvt <- cvt%>%group_by(comb)%>%mutate(yscale=y-mean(y,na.rm=T))%>%ungroup()
   }else{
      cvt <- cvt%>%mutate(yscale=y)
   }
}
###

#c0r <- load("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.ComBat_seq.C0.SES_PCs_sex_age_and_treats_adjusted.RData")
#writing R output to this textfile
#sink(paste0(platePrefix,"new_norm_cd52_pnsi.txt"))

# _pulm_functions
timestamp()


for(i in 1:length(variables)){
  myvar <- variables[i]
  var_name <- variable_names[i]
  cat("######################################################################\n")
  cat("##Processing: ", var_name, "\n")
  
  ## assign variables, load data, and load experiment information
  topDirectory <- getwd()
  outDir <- paste(topDirectory, "/",myvar, sep='')
  system(paste("mkdir -p",outDir))
  ##
  #plotsDir <- paste(outDir, '/plots', sep='')
  #system(paste("mkdir -p", plotsDir))
  ##
  #statsDir <- paste(outDir, '/stats', sep='')
  #system(paste("mkdir -p", statsDir))
  ##
  #dataDir <- paste(outDir, '/data_objects', sep='')
  #system(paste("mkdir -p", dataDir))
  
  ##################################################################
  

  # Initialize list to hold the data frames for the current variable
  dfvar <- list()
  dfsig <- list()

  ################ Inner Loop: Perform erichGO for each cluster ###################
  for (cluster in 0:4){
    #tryCatch({
    #fname <- paste(dataDir, platePrefix, ".C", cluster, ".deseqres_", myvar, plateSufix, ".txt", sep="")
    fname <- paste(dataDir, platePrefix, "C", cluster, ".deseqres_", myvar, plateSufix, ".txt", sep = "")

    # Read in the data
    res_data <- read.table(fname, sep="\t", header=T, quote='"', comment="")
    resSig <- res_data %>% filter(padj < 0.1)
    rownames(resSig) <- resSig$identifier
    resSig <- resSig[order(resSig$padj), ]

    # Store each data frame in the list with the name "resX"
    dfvar[[paste0(cluster)]] <- res_data
    dfsig[[paste0(cluster)]] <- res_data %>% filter(padj < 0.1)

    # read count data
    #ALL.ComBat_seq.C1.SES_PCs_sex_age_and_treats_adjusted.RData
    #cname <- paste(dataDir, platePrefix, ".ComBat_seq", ".C", cluster, countSufix, ".RData", sep="")
    #load(cname)
    #count_data <- adjusted_counts

    # read DESeq object for each cluster
    #deseqvar <- readRDS("ALL.DESeq_output-RNA-CTRL-iselC0.SES_PCs_sex_age_and_treats_adjusted.RDS")
    #dname <- paste(deseqDir, platePrefix, ".DESeq_output-RNA-CTRL-",myvar, "C", cluster, deseqSuffix, ".RDS", sep="")
    #dname <- paste(deseqDir, platePrefix, "DESeq_output-RNA-CTRL-",myvar, "C", cluster, deseqSuffix, ".RDS", sep="")
    dname <- paste(deseqDir, platePrefix, "DESeq_output-RNA-CTRL-",myvar, "C", cluster, deseqSuffix, ".RData", sep="")
    load(dname)    
   
    dname <- paste(normdata, "C", cluster, "_CTRL_rna.normal", ".rds", sep="")
    # C0_CTRL_rna.normal.rds
    normcnt <- readRDS(dname)

    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS
    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS    
    cat("## total genes present in cluster", cluster, "for", var_name,":", nrow(dfvar[[paste0(cluster)]]),  "\n")
    cat("## numb DEGs present in cluster", cluster, "for", var_name,":", nrow(dfsig[[paste0(cluster)]]),  "\n")
    #cat("## number of genes for ", cluster, " count data:", nrow(counts),  "\n")
    cat("## dimention of DEseq count ", cluster, "for", myvar, dim(counts(dds)),  "\n")
    cat("## dimention of DEseq count ", cluster, "for", myvar, dim(normcnt),  "\n")

    cat("#################################################\n")


    head(rownames(normcnt))
    head(rownames(dds))

    head(colnames(normcnt))
    head(colnames(dds))

  #}
  
#/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/deseqres/ALL.0.2.11..DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_generem.RDS
#/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/
  cat("Now making plots\n")
  ###Make gene plots colored by wave
  ##normalized gene counts plots for top 10 significant genes
  ##NK cells
  cat("Generating significant gene scatter plots for cluster ", cluster, " - variable ", myvar, "\n")
  if(nrow(resSig) > 0 ) {
    DEGs <- rownames(resSig)
    stopifnot(all(DEGs %in% names(dds)))
    if(length(DEGs) > 20) {
      DEGs <- DEGs[1:20]
    }
    # normalized counts for the genes
    #tcounts <- t(log2((counts(dds[DEGs, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    tcounts <- t(normcnt[DEGs, ]) %>%
      merge(colData(dds), ., by="row.names") %>%
      gather(gene, expression, (ncol(.)-length(DEGs)+1):ncol(.))
    fname=paste(outDir, '/', "cluster_", cluster, "_", myvar, "_top10_significant_genes", ".pdf", sep="")
    pdf(fname)
    for(j in DEGs) {
      #geneID <- substr(j, 0, 15)
      geneID <- j
      #geneSymbol <- grch38$symbol[grch38$ensgene == geneID]
        df_gene <- filter(tcounts, gene == j)

        regressed <- lm(df_gene$expression ~ df_gene$age + df_gene$sex_alph + df_gene$BATCH + df_gene$PC1 + df_gene$PC2)
        df_gene$residual <- residuals(regressed)
        df_gene$expression <-  df_gene$residual

        # Calculate Pearson correlation
        cor_test <- cor.test(df_gene[[myvar]], df_gene$expression)
        cor_value <- round(cor_test$estimate, 2)  # Pearson correlation
        p_value <- format(cor_test$p.value, format = "f", digits = 4)  # P-value

      p <- ggplot(df_gene, aes(get(myvar), expression)) +
        geom_point(position=position_jitter(w=0.1, h=0)) +
        geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
        theme_minimal() +
        #ggtitle(geneID) +
        ggtitle(geneID, subtitle = paste0("r=", cor_value, ", p=", p_value)) +
        labs(x=var_name,
             y="Residual Expression") +
        theme(plot.title = element_text(hjust=0.5, size = rel(1.3)))
      print(p)
    }     
    dev.off() 
  }
  }
}




######################################################################
#### scatter plot for a single gene of interest
######################################################################

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/gene/residual/singlegene/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/singlegenes/residual/cytocomp_nom_fdr_2var/")

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/gene/nominal/")


#variables <- c("pr_comp","isel","PSS_all_mean")
#variable_names <- c("Psychological Resources", "Social Support", "Perceived Stress")

#variables <- c("chronic_sum", "cytocomp", "BPd_avg", "SES")
#variable_names <- c("Chronic conditions",  "Cytokines composite", "SES composite")

variables <- c("PSS_all_mean", "ISEL_Mean", "cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
variable_names <- c("Psychological Stress", "Social Support", "Cytokines") #"socioeconomic status", "diastolic blood pressure", 


#variables <- c("PSS_all_mean", "ISEL_Mean", "z_ifny_0_log_w" )#,"cytocomp") # "pr_comp", "SNI_NoP")#, "SES", "BPd_avg", 
#variable_names <- c("Psychological Stress", "Social Support","IFN-γ")# "Cytokines") #"socioeconomic status", "diastolic blood pressure", 

mygenes <- c("STAT1", "STAT2")#, "STAT3", "JAK1", "JAK2", "IRF2", "IRF9", "INFA1", "IFNB1", "IFNGR1", "IFNGR2", "IFIT1", "IFITM2", "IFNA21", "IFI44", "ISG15")
mygenes <- c("ISG15", "STAT2")#, "STAT3", "JAK1", "JAK2", "IRF2", "IRF9", "INFA1", "IFNB1", "IFNGR1", "IFNGR2", "IFIT1", "IFITM2", "IFNA21", "IFI44", "ISG15")
mygenes <- c("GBP4")#, "STAT3", "JAK1", "JAK2", "IRF2", "IRF9", "INFA1", "IFNB1", "IFNGR1", "IFNGR2", "IFIT1", "IFITM2", "IFNA21", "IFI44", "ISG15")


mygenes <- c("AC092821.3", "USP18",  "INFA1", "IFNB1", "LIPA")
# C0 genes: code bellow to get most signficant shared genes 
#mygenes <- topC0


# C0 genes: code bellow to get most signficant shared genes 
mygenes <- topC0
#mygenes <- topC4
#mygenes <- topC3
#mygenes <- topC2
#mygenes <- topC1

#mygenes <- ifnygenes 

#myensemble <- "ENSG00000169442"
#variables <- c("pnsi")
#variable_names <- c("Neighborhood Stress")

# IFN-α, IFN-β, IFIT1, IFIT3, IFITM1, IFITM2, IFNA21, IFNA5, IRF1, and IRF9, INFA1, IFNB1
# IFNGR1, IFNGR2, JAK1, JAK2, STAT1, STAT2, and STAT3

# _pulm_functions
timestamp()


for(i in 1:length(variables)){
     # tryCatch({

  myvar <- variables[i]
  var_name <- variable_names[i]
  cat("######################################################################\n")
  cat("##Processing: ", var_name, "\n")
  
  ## assign variables, load data, and load experiment information
  topDirectory <- getwd()
  outDir <- paste(topDirectory, "/",myvar, sep='')
  system(paste("mkdir -p",outDir))
  ##
  #plotsDir <- paste(outDir, '/plots', sep='')
  #system(paste("mkdir -p", plotsDir))
  ##
  #statsDir <- paste(outDir, '/stats', sep='')
  #system(paste("mkdir -p", statsDir))
  ##
  #dataDir <- paste(outDir, '/data_objects', sep='')
  #system(paste("mkdir -p", dataDir))
  
  ##################################################################
  

  # Initialize list to hold the data frames for the current variable
  dfvar <- list()
  dfsig <- list()

  ################ Inner Loop: Perform erichGO for each cluster ###################
  for (cluster in 0){
    #tryCatch({
    
    #outDir <- paste(outDir, "/",myvar,"/C", cluster,  sep='')
    #system(paste("mkdir -p",outDir))
 
    fname <- paste(dataDir, platePrefix, "C", cluster, ".deseqres_", myvar, plateSufix, ".txt", sep = "")
    
    # Read in the data
    res_data <- read.table(fname, sep="\t", header=T, quote='"', comment="")
    #resSig <- res_data %>% filter(padj < 0.1)
    resSig <- res_data 
    rownames(resSig) <- resSig$identifier
    resSig <- resSig[order(resSig$padj), ]

    # Store each data frame in the list with the name "resX"
    dfvar[[paste0(cluster)]] <- res_data
    dfsig[[paste0(cluster)]] <- res_data %>% filter(padj < 0.1)

    # read count data
    #ALL.ComBat_seq.C1.SES_PCs_sex_age_and_treats_adjusted.RData
    #cname <- paste(dataDir, platePrefix, ".ComBat_seq", ".C", cluster, countSufix, ".RData", sep="")
    #load(cname)
    #count_data <- adjusted_counts

    # read DESeq object for each cluster
    #deseqvar <- readRDS("ALL.DESeq_output-RNA-CTRL-iselC0.SES_PCs_sex_age_and_treats_adjusted.RDS")
    #dname <- paste(deseqDir, platePrefix, "DESeq_output-RNA-CTRL-",myvar, "C", cluster, deseqSuffix, ".RDS", sep="")
    #dds <- readRDS(dname)

    dname <- paste(deseqDir, platePrefix, "DESeq_output-RNA-CTRL-",myvar, "C", cluster, deseqSuffix, ".RData", sep="")
    load(dname)    
   
    dname <- paste(normdata, "C", cluster, "_CTRL_rna.normal", ".rds", sep="")
    # C0_CTRL_rna.normal.rds
    normcnt <- readRDS(dname)

    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS
    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS    
    cat("## total genes present in cluster", cluster, "for", var_name,":", nrow(dfvar[[paste0(cluster)]]),  "\n")
    cat("## numb DEGs present in cluster", cluster, "for", var_name,":", nrow(dfsig[[paste0(cluster)]]),  "\n")
    cat("## number of genes for ", cluster, " count data:", nrow(counts),  "\n")
    cat("## dimention of DEseq count ", cluster, "for", myvar, dim(counts(dds)),  "\n")
    cat("## dimention of DEseq count ", cluster, "for", myvar, dim(normcnt),  "\n")

    cat("#################################################\n")

  #}
  

  cat("Now making plots\n")
  ###Make gene plots colored by wave
  ##normalized gene counts plots for top 10 significant genes
  ##NK cells
  cat("Generating significant gene scatter plots for cluster ", cluster, " - variable ", myvar, "\n")

    for (g in 1:length(mygenes)) {
          #tryCatch({
        mygene <- mygenes[g]
         cat("gene ", mygene, "\n")

  if(nrow(resSig) > 0 ) {
    DEGs <- rownames(resSig)
    stopifnot(all(DEGs %in% names(dds)))
    if(length(DEGs) > 10) {
      DEGs <- DEGs[DEGs == mygene]
    }
    # normalized counts for the genes
    #tcounts <- t(log2((counts(dds[DEGs, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    tcounts <- t(normcnt[DEGs, ]) #%>%
      tcounts <- setNames(as.data.frame(normcnt[DEGs, ]), DEGs) %>%
      merge(colData(dds), ., by="row.names") %>% as.data.frame() %>%
      gather(mygene, expression, (ncol(.)-length(DEGs)+1):ncol(.))

      #tcounts <- t(normcnt[DEGs, , drop = FALSE]) %>%
      #merge(colData(dds), ., by = "row.names") %>%
      #rename(expression = all_of(DEGs))
      tcounts$gene <- DEGs

    fname=paste(outDir, '/', "C", cluster, "_",mygene, "_", myvar, "_significant_gene_", mygene, ".png", sep="")
    #pdf(fname)
    png(fname, width = 2200, height = 2200, res = 300)
    #png(fname, width=600, height=600, res=120)        

    for(j in DEGs) {
      #geneID <- substr(j, 0, 15)
      geneID <- j
        df_gene <- filter(tcounts, gene == j)

        regressed <- lm(df_gene$expression ~ df_gene$age + df_gene$sex_alph + df_gene$BATCH + df_gene$PC1 + df_gene$PC2)
        df_gene$residual <- residuals(regressed)
        df_gene$expression <-  df_gene$residual

              # Calculate Pearson correlation
        cor_test <- cor.test(df_gene[[myvar]], df_gene$expression)
        cor_value <- round(cor_test$estimate, 2)  # Pearson correlation
        p_value <- format(cor_test$p.value, format = "f", digits = 4)  # P-value

      #geneSymbol <- grch38$symbol[grch38$ensgene == geneID]
      p <- ggplot(df_gene, aes(get(myvar), expression)) +
        #geom_point(position=position_jitter(w=0.1, h=0)) +
        #geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
        geom_point(size = 3, alpha = 0.8, color = "steelblue", 
             position = position_jitter(width = 0.15, height = 0.15)) +
        geom_smooth(method = 'lm', color = "black", se = TRUE, size = 1.2) +
        #theme_minimal() +
        #ggtitle(geneID) +
        #ggtitle(paste0(geneID, " - C", cluster), subtitle = paste0("r=", cor_value, ", p=", p_value)) +
        #labs(x=var_name,
        #     y="Normalized Expression") +
        #theme(plot.title = element_text(hjust=0.5, size = rel(1.3)))

        labs(
        title = paste0(geneID, " - C", cluster),
        subtitle = paste0("r=", cor_value, ", p=", p_value),
        x = var_name,
        y = paste0(geneID, " Normalized Expression")
          ) +
          theme_minimal(base_size = 16) +
          theme(
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none"
      )

      print(p)
    }     
    dev.off() 
  }
  
  
  }
 }
}




#            }, 
#            error = function(e) {
#             #Handle the error, e.g., log it or print a message
#            print(paste("first mess:Error at iteration for variable ", myvar ," in cell type ", cluster, ":", e$message))
#            },
#            finally = {
#            # Optional: code to run after each iteration, regardless of success or error
#            print(paste("first mess: Completed iteration for variable ", myvar ," in cell type ", cluster))
#            })
#  }
# }
#}






############### shared top genes nominal for cytocomp and padj sig for PSS and ISEL

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/gene/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 


#countSufix <- ".SES_PCs_sex_age_and_treats_adjusted"
deseqSuffix <- ".noCombat_DESeq"
deseqDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/"

### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

resall <-  resall %>% filter(var %in% target_vars)

library(dplyr)

# Step 1. Separate the three variable types
isel <- resall %>%
  filter(var == "ISEL_Mean", padj < 0.1) %>%
  dplyr::select(identifier, cluster)

pss <- resall %>%
  filter(var == "PSS_all_mean", padj < 0.1) %>%
  dplyr::select(identifier, cluster)

cyto <- resall %>%
  filter(var == "cytocomp", pvalue < 0.05) %>%
   dplyr::select(identifier, cluster)

# Step 2. Find overlap per cluster
sig_overlap <- isel %>%
  inner_join(pss, by = c("identifier", "cluster")) %>%
  inner_join(cyto, by = c("identifier", "cluster"))

# Step 3. Extract full rows from resall for those overlapping genes
sig_resall <- resall %>%
  semi_join(sig_overlap, by = c("identifier", "cluster"))


topC0 <- sig_resall %>% filter(cluster=="C0") %>% pull(identifier) %>% unique() %>% head(10)
topC4 <- sig_resall %>% filter(cluster=="C4") %>% pull(identifier) %>% unique() %>% head(10)
topC3 <- sig_resall %>% filter(cluster=="C3") %>% pull(identifier) %>% unique() %>% head(15)
topC2 <- sig_resall %>% filter(cluster=="C2") %>% pull(identifier) %>% unique() %>% head(10)
topC1 <- sig_resall %>% filter(cluster=="C1") %>% pull(identifier) %>% unique() %>% head(10)


# Step 4. Split into cluster-specific dataframes
cluster_list <- split(sig_resall, sig_resall$cluster)
list2env(cluster_list, envir = .GlobalEnv)


# Save the combined dataframe as a single text file
fname <- paste0(getwd(), "/top_shared_DEGs_3var_pvalue_nomin.txt")
write.table(topgenes, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)





############### shared top genes:

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/gene/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.1.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 

target_vars <- c("ISEL_Mean", "PSS_all_mean", "cytocomp")

#countSufix <- ".SES_PCs_sex_age_and_treats_adjusted"
deseqSuffix <- ".noCombat_DESeq"
deseqDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/"

### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

resall <-  resall %>% filter(var %in% target_vars)

resdegs <- resall %>% filter(padj < 0.1) %>% filter(var %in% target_vars)
library(dplyr)

resdegs %>% filter(identifier=="ATG3")
resdegs %>% filter(identifier=="RBX1")
resdegs %>% filter(identifier=="SAP18")

# Step 1: Define the target variables
#target_vars <- c("ISEL_Mean", "PSS_all_mean", "z_ifny_0_log_w")

# Step 2: Filter resdegs to only these variables
filtered <- resdegs %>%
  filter(var %in% target_vars)

# Step 3: Create a unique key: gene + celltype
filtered <- filtered %>%
  mutate(gene_celltype = paste(identifier, celltype, sep = "::"))

filtered <- as.data.frame(filtered)

# Step 4: For each gene+celltype, count how many of the 3 variables it's found in
#shared_gene_ct <- filtered %>%
#  distinct(gene_celltype, var) %>%
#  count(gene_celltype, name = "n_vars") %>%
#  #filter(n_vars == 3)
#  filter(n_vars == 3)

  shared_gene_ct <- filtered %>%
  distinct(gene_celltype, var) %>%
  group_by(gene_celltype) %>%
  summarise(n_vars = n_distinct(var), .groups = "drop") %>%
  filter(n_vars == 2)


# Step 5: Filter original data to keep only rows for shared gene+celltype
shared_data <- filtered %>%
  filter(gene_celltype %in% shared_gene_ct$gene_celltype)

# Step 6: (Optional) Separate identifier and celltype back out
topgenes <- shared_data %>%
  tidyr::separate(gene_celltype, into = c("identifier", "celltype"), sep = "::", remove = FALSE) %>%
  as.data.frame() %>%
  dplyr::select(identifier, padj, var, celltype, cluster, logFC, SE) 
   
topgenes <- topgenes[order(topgenes$padj), ]
#topgenes <- topgenes[order(topgenes$padj), ]


# Save the combined dataframe as a single text file
fname <- paste0(getwd(), "/top_shared_DEGs_ISEL_PSS_ISEL_ifny_padj.txt")
write.table(topgenes, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


topC0 <- topgenes %>% filter(cluster=="C0") %>% pull(identifier) %>% unique() %>% head(10)
topC4 <- topgenes %>% filter(cluster=="C4") %>% pull(identifier) %>% unique() %>% head(10)
topC3 <- topgenes %>% filter(cluster=="C3") %>% pull(identifier) %>% unique() %>% head(10)
topC2 <- topgenes %>% filter(cluster=="C2") %>% pull(identifier) %>% unique() %>% head(10)
topC1 <- topgenes %>% filter(cluster=="C1") %>% pull(identifier) %>% unique() %>% head(10)

#ifnygenes <- topgenes %>% filter(var == "z_ifny_0_log_w") %>% pull(identifier)

table(topgenes$var)
table(topgenes$celltype)

    topgenes <- topgenes[order(topgenes$padj), ]
unique(head(topgenes, n=20)$identifier)

#AC004551.1

topgenes %>% filter(identifier=="AC004551.1")




###### fin the top nominal genes to include with cytokines: 


############### shared top genes:
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/gene/nominal/")

#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/gene/")
dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

platePrefix <- "ALL.0.15.13." 
plateSufix <- "-RNA-CTRL.noCombat_DESeq" 


#countSufix <- ".SES_PCs_sex_age_and_treats_adjusted"
deseqSuffix <- ".noCombat_DESeq"
deseqDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/"

### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

resdegs <- resall %>% filter(pvalue < 0.05)
library(dplyr)

# Step 1: Define the target variables
target_vars <- c("ISEL_Mean", "PSS_all_mean", "cytocomp")

# Step 2: Filter resdegs to only these variables
filtered <- resdegs %>%
  filter(var %in% target_vars)

# Step 3: Create a unique key: gene + celltype
filtered <- filtered %>%
  mutate(gene_celltype = paste(identifier, celltype, sep = "::"))

filtered <- as.data.frame(filtered)

  shared_gene_ct <- filtered %>%
  distinct(gene_celltype, var) %>%
  group_by(gene_celltype) %>%
  summarise(n_vars = n_distinct(var), .groups = "drop") %>%
  filter(n_vars == 3)


# Step 5: Filter original data to keep only rows for shared gene+celltype
shared_data <- filtered %>%
  filter(gene_celltype %in% shared_gene_ct$gene_celltype)

# Step 6: (Optional) Separate identifier and celltype back out
topgenes <- shared_data %>%
  tidyr::separate(gene_celltype, into = c("identifier", "celltype"), sep = "::", remove = FALSE) %>%
  as.data.frame() %>%
  dplyr::select(identifier, padj, pvalue, var, celltype, cluster, logFC, SE) 
   
topgenes <- topgenes[order(topgenes$pvalue), ]
#topgenes <- topgenes[order(topgenes$padj), ]

topgenes <- topgenes %>% arrange(pvalue)



# Save the combined dataframe as a single text file
fname <- paste0(getwd(), "/top_shared_DEGs_3var_pvalue_nomin.txt")
write.table(topgenes, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


topC0 <- topgenes %>% filter(cluster=="C0") %>% pull(identifier) %>% unique() %>% head(10)
topC4 <- topgenes %>% filter(cluster=="C4") %>% pull(identifier) %>% unique()
topC3 <- topgenes %>% filter(cluster=="C3") %>% pull(identifier) %>% unique() %>% head(10)
topC2 <- topgenes %>% filter(cluster=="C2") %>% pull(identifier) %>% unique() %>% head(10)
topC1 <- topgenes %>% filter(cluster=="C1") %>% pull(identifier) %>% unique() %>% head(10)







iselgenes <- ressig %>% filter(var == "ISEL_Mean") %>% pull(identifier) %>% unique()
pssgenes <- ressig %>% filter(var == "PSS_all_mean") %>% pull(identifier) %>% unique()

length(intersect(iselgenes, pssgenes))





