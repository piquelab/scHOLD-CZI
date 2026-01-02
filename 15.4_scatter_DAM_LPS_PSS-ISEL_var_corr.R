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

#datasc <- readRDS("/rs/rs_grp_schold/CZI/RNA/analysis/5b_IdenCelltype_fastdemux/ALL.seuratObj-.harmony-sctype-2024-05-15.rds")
#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/scatter/CTRLvsLPS/interaction_DEGs/interaction_and_CTRL")

setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/")
#setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/list1/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/topsigall/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/LPSgenes/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/topISEL/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/topPSS/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/no_voom/opptop/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/all_treats/")
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS/DAMs/")

#dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"
#normDir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data"
#normDir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data/no_voom"
#normDir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data/alltreat"

fname <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/option_nFeature15K_cluster_res0.12/1.3_YtX_ave.th20.clean.rds"
motifmat <- readRDS(fname)

# make metadata
meta <- tibble(sample = colnames(motifmat)) %>%
  separate(sample, into = c("cluster", "treat", "dbgap.ID"), sep = "_", remove = FALSE) %>% as.data.frame()

rownames(meta) <- meta$sample

# assign atac metadata: 
meta$celltype <- "NA"
meta$celltype[meta$Cluster == "C0"] <- "A0 T CD4+"
meta$celltype[meta$Cluster == "C1"] <- "A1 T CD8+"
meta$celltype[meta$Cluster == "C2"] <- "A2 NK"
meta$celltype[meta$Cluster == "C3"] <- "A3 T CD4+"
meta$celltype[meta$Cluster == "C4"] <- "A4 Monocyte"
meta$celltype[meta$Cluster == "C5"] <- "A5 T CD4+"
meta$celltype[meta$Cluster == "C6"] <- "A6 B"
meta$celltype[meta$Cluster == "C7"] <- "A7 T CD4+"
meta$celltype[meta$Cluster == "C8"] <- "A8 T CD4+"
meta$celltype[meta$Cluster == "C9"] <- "A9 T CD4+"
meta$celltype[meta$Cluster == "C10"] <- "A10 DC"

   normCTRL <- motifmat[, grep("CTRL", colnames(motifmat))]
   normLPS <- motifmat[, grep("LPS", colnames(motifmat))]


# load DAMs for variables
fname <-  "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
vardams <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
vardams <- vardams %>% filter(psycho_variable %in% c("PSS_all_mean", "ISEL_Mean")) %>% filter(padj_t < 0.1)
vardams$identifier <- vardams$gene 
vardams$gene_cluster <- paste0(vardams$identifier, "_", vardams$Cluster)

# assign atac metadata: 
vardams$celltype <- "NA"
vardams$celltype[vardams$Cluster == "C0"] <- "A0 T CD4+"
vardams$celltype[vardams$Cluster == "C1"] <- "A1 T CD8+"
vardams$celltype[vardams$Cluster == "C2"] <- "A2 NK"
vardams$celltype[vardams$Cluster == "C3"] <- "A3 T CD4+"
vardams$celltype[vardams$Cluster == "C4"] <- "A4 Monocyte"
vardams$celltype[vardams$Cluster == "C5"] <- "A5 T CD4+"
vardams$celltype[vardams$Cluster == "C6"] <- "A6 B"
vardams$celltype[vardams$Cluster == "C7"] <- "A7 T CD4+"
vardams$celltype[vardams$Cluster == "C8"] <- "A8 T CD4+"
vardams$celltype[vardams$Cluster == "C9"] <- "A9 T CD4+"
vardams$celltype[vardams$Cluster == "C10"] <- "A10 DC"
#trl_deseqres_sig <- vardams %>% filter(psycho_variable == "PSS_all_mean")
#iselctrl_deseqres_sig <- vardams %>% filter(psycho_variable == "ISEL_Mean")

subset(vardams, grepl("IRF4", gene, ignore.case = TRUE))
subset(vardams, grepl("IRF2", gene, ignore.case = TRUE))
subset(vardams, grepl("REL", gene, ignore.case = TRUE))
subset(vardams, grepl("IRF7", gene, ignore.case = TRUE))
subset(vardams, grepl("IRF7", gene, ignore.case = TRUE))
subset(vardams, grepl("STAT1", gene, ignore.case = TRUE))

subset(treats, grepl("STAT1", gene, ignore.case = TRUE))
subset(treats, grepl("IRF4", gene, ignore.case = TRUE))
subset(treats, grepl("ARNT..HIF1A_26", gene, ignore.case = TRUE))
subset(treats, grepl("ERF..NHLH1_584", gene, ignore.case = TRUE))


varmon <- vardams %>% filter(Cluster=="C3")
subset(varmon, grepl("STAT", gene, ignore.case = TRUE))
subset(varmon, grepl("ARNT..HIF1A_26", gene, ignore.case = TRUE))
subset(varmon, grepl("RELA_16", gene, ignore.case = TRUE))
subset(varmon, grepl("REL_15", gene, ignore.case = TRUE))
subset(varmon, grepl("ERF..NHLH1_584", gene, ignore.case = TRUE))


#dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/"
timestamp()

#ALL.ComBat_seq.C1.SES_PCs_sex_age_and_treats_adjusted.RData


#platePrefix <- "ALL.0.15.13." 
#plateSufixC <- "-RNA-CTRL.noCombat_DESeq" 
#plateSufixL <- "-RNA-LPS.noCombat_DESeq" 
#plateSufixI <- "RNA.CTRL.SES_PCs_sex_age_and_treats_generem_treatint" 

#dataDirI <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem_treatint/deseqres/"

#countSufix <- ".SES_PCs_sex_age_and_treats_adjusted"
#deseqSuffix <- ".noCombat_DESeq"
#deseqSuffixI <- ".SES_PCs_sex_age_and_treats_generem_treatint"
#deseqDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/"
#deseqDirI <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem_treatint/"


variables <- c("PSS_all_mean", "ISEL_Mean")#, "cytocomp")
variable_names <- c("Perceived Stress",  "Social Support")#, "Cytokine Composite")


mygenes <- c("IRF4_308", "IRF2_3", "IRF7_157", "STAT1..STAT2_36", "STAT1_41", "ARNT..HIF1A_26", "ERF..NHLH1_584")#, "RSAD2", "CMPK2", "IRF4", "CXCL11", "MT1F", "P2RY14")

#ARNT..HIF1A_26
#ERF..NHLH1_584
#mygenes <- c("CXCL10")#, "RSAD2", "CMPK2", "IRF4", "CXCL11", "MT1F", "P2RY14")

#variables <- c("PSS_all_mean")#, "ISEL_Mean")#, "cytocomp")
#variable_names <- c("Perceived Stress")#,  "Social Support")#, "Cytokine Composite")

#myensemble <- "ENSG00000169442"
#variables <- c("pnsi")
#variable_names <- c("Neighborhood Stress")

# IFN-α, IFN-β, IFIT1, IFIT3, IFITM1, IFITM2, IFNA21, IFNA5, IRF1, and IRF9, INFA1, IFNB1
# IFNGR1, IFNGR2, JAK1, JAK2, STAT1, STAT2, and STAT3

# _pulm_functions
timestamp()


#atacdir = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/"
#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/psycho_var_dir
atacdir = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/"


setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/scatter/CTRLvsLPS_DAM/")


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
  statsDir <- paste(outDir, '/stats', sep='')
  system(paste("mkdir -p", statsDir))
  ##
  #dataDir <- paste(outDir, '/data_objects', sep='')
  #system(paste("mkdir -p", dataDir))
  
  ##################################################################
  

  # Initialize list to hold the data frames for the current variable
  dfvarC <- list()
  dfsigC <- list()

  dfvarL <- list()
  dfsigL <- list()

  dfvarI <- list()
  dfsigI <- list()

  ################ Inner Loop: Perform erichGO for each cluster ###################
  for (cluster in c(0,4)){
    #tryCatch({

    cat("######## loading CTRL for cluster ", cluster, "\n")
    # Read in the data CTRL
    #fname <- paste(dataDir, platePrefix, "C", cluster, ".deseqres_", myvar, plateSufixC, ".txt", sep = "")
    #fname <- paste(dataDir, platePrefix, "C", clust, ".deseqres_", myvar, plateSufix, ".txt", sep = "")

    #res_dataC <- read.table(fname, sep="\t", header=T, quote='"', comment="")
    #resSigC <- res_dataC %>% filter(padj < 0.1)
    #rownames(resSigC) <- resSigC$identifier
    #resSigC <- resSigC[order(resSigC$padj), ]

    # Store each data frame in the list with the name "resX"
    #dfvarC[[paste0(cluster)]] <- res_dataC
    #dfsigC[[paste0(cluster)]] <- res_dataC %>% filter(padj < 0.1)

    #cat("## numb DEGs present in cluster", cluster, "for", var_name,":", nrow(dfsigC[[paste0(cluster)]]),  "\n")

    # read count data
    #ALL.ComBat_seq.C1.SES_PCs_sex_age_and_treats_adjusted.RData
    #cname <- paste(dataDir, platePrefix, ".ComBat_seq", ".C", cluster, countSufix, ".RData", sep="")
    #load(cname)
    #count_data <- adjusted_counts

    # read DESeq object for each cluster
    #deseqvar <- readRDS("ALL.DESeq_output-RNA-CTRL-iselC0.SES_PCs_sex_age_and_treats_adjusted.RDS")
    #ALL.0.15.13.DESeq_output-RNA-CTRL-PSS_all_meanC3.noCombat_DESeq.RData
    #dname <- paste(deseqDir, platePrefix, "DESeq_output-RNA-CTRL-",myvar, "C", cluster, deseqSuffix, ".RData", sep="")
    #load(dname)
    #ddsC <- dds
    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS
    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS    
    #cat("## total genes present in cluster", cluster, "for", var_name,":", nrow(dfvarC[[paste0(cluster)]]),  "\n")
    #cat("## numb DEGs present in cluster", cluster, "for", var_name,":", nrow(dfsigC[[paste0(cluster)]]),  "\n")
    ##cat("## number of genes for ", cluster, " count data:", nrow(counts),  "\n")
    #cat("## dimention of DEseq count ", cluster, "for", myvar, dim(counts(ddsC)),  "\n")

    # load normalized count data
    #fname = paste0(normDir, "/C", cluster, "_CTRL_rna.normal.rds")
    #fname = paste0(normDir, "/C", cluster, "_CTRL_rna.normal.no_voom.rds")
    #fname = paste0(normDir, "/C", cluster, "_all_rna.normal.rds")
    #norm <- readRDS(fname)

    # filter DAMs for the cluster

    #filter motif activity for the cluster
    normC <- normCTRL[, grep(paste0("C",cluster), colnames(normCTRL))]
    normL <- normLPS[, grep(paste0("C",cluster), colnames(normLPS))]

    # ctrl 
    x0 <- as.data.frame(str_split(colnames(normC), "_", simplify = T))
    names(x0) <- c("cluster", "treat", "sampleID")
    x0$bti <- colnames(normC)
    #x <- x0%>%dplyr::filter(cluster==cluster, treat=="CTRL")
    motif_dat2C <- normC[, x0$bti]
    colnames(motif_dat2C) <- x0$sampleID

    #LPS
    x0 <- as.data.frame(str_split(colnames(normL), "_", simplify = T))
    names(x0) <- c("cluster", "treat", "sampleID")
    x0$bti <- colnames(normL)
    #x <- x0%>%dplyr::filter(cluster==cluster, treat=="CTRL")
    motif_dat2L <- normL[, x0$bti]
    colnames(motif_dat2L) <- x0$sampleID


    ### covariates
    fn0 <-paste0(atacdir, "psycho_var_dir/0.1_covariates_varSel.txt")
    cv0 <- read.table(fn0, header=T, sep="\t")
    colSel <- c("sampleID", "Batch", "sex_alph", "PC1", "PC2", "age", myvar)
    cv <- cv0%>%dplyr::select(all_of(colSel))%>%drop_na(all_of(colSel))



    ###
    ### combine data together

    #shared <- Reduce(intersect, list(colnames(motif_dat2), colnames(rna_norm), colnames(TF_score), cv$sampleID))
    shared <- Reduce(intersect, list(colnames(motif_dat2C),  colnames(motif_dat2L), cv$sampleID)) #yscore_df2$sampleID, 


               #motif_regulated_score=TF_score[motif0, sampleID], 
               #motif_regulated_score_all=TF_score_all[motif0, sampleID])
    #motif0="IRF4_308"

 
               #pathway_score_2=pathscore2[pID2, sampleID])
    #cv2 <- cbind(cv2, t(rna_norm[motif00, cv2$sampleID, drop=F]))


    #cat("## merged four data completed: ", motif00, "\n")


    #normC <- readRDS(fname)

    #C0_CTRL_rna.normal.rds
    cat("################# CTRL LOADED ##################\n")

    cat("######## loading LPS for cluster ", cluster, "\n")
    # Read in the data LPS
    #fname <- paste(dataDir, platePrefix, "C", cluster, ".deseqres_", myvar, plateSufixL, ".txt", sep = "")
    
    #res_dataL <- read.table(fname, sep="\t", header=T, quote='"', comment="")
    #resSigL <- res_dataL %>% filter(padj < 0.1)
    #rownames(resSigL) <- resSigL$identifier
    #resSigL <- resSigL[order(resSigL$padj), ]

    # Store each data frame in the list with the name "resX"
    #dfvarL[[paste0(cluster)]] <- res_dataL
    #dfsigL[[paste0(cluster)]] <- res_dataL %>% filter(padj < 0.1)

    # read count data
    #ALL.ComBat_seq.C1.SES_PCs_sex_age_and_treats_adjusted.RData
    #cname <- paste(dataDir, platePrefix, ".ComBat_seq", ".C", cluster, countSufix, ".RData", sep="")
    #load(cname)
    #count_data <- adjusted_counts

    # read DESeq object for each cluster
    #deseqvar <- readRDS("ALL.DESeq_output-RNA-CTRL-iselC0.SES_PCs_sex_age_and_treats_adjusted.RDS")
    #dname <- paste(deseqDir, platePrefix, "DESeq_output-RNA-LPS-",myvar, "C", cluster, deseqSuffix, ".RData", sep="")
    #load(dname)
    #ddsL <- dds

    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS
    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS    
    #cat("## total genes present in cluster", cluster, "for", var_name,":", nrow(dfvarL[[paste0(cluster)]]),  "\n")
    #cat("## numb DEGs present in cluster", cluster, "for", var_name,":", nrow(dfsigL[[paste0(cluster)]]),  "\n")
    #cat("## number of genes for ", cluster, " count data:", nrow(counts),  "\n")
    #cat("## dimention of DEseq count ", cluster, "for", myvar, dim(counts(ddsL)),  "\n")

    # load normalized count data
    #fname = paste0(normDir, "/C", cluster, "_LPS_rna.normal.rds")
    #fname = paste0(normDir, "/C", cluster, "_LPS_rna.normal.no_voom.rds")
    #normL <- readRDS(fname)

    cat("##################### LPS LOADED ############################\n")


  cat("Now making plots\n")
  ###Make gene plots colored by wave
  ##normalized gene counts plots for top 10 significant genes
  ##NK cells
  cat("Generating significant gene scatter plots for cluster ", cluster, " - variable ", myvar, "\n")

    for (g in 1:length(mygenes)) {
          #tryCatch({
        mygene <- mygenes[g]
         cat("gene ", mygene, "\n")
         motif0 <- mygene


    #cv2$gene <- motif0

    #cv2 <- cv%>%filter(sampleID%in%shared)%>%
    #            mutate(motif_activity_CTRL=motif_dat2C[motif0, sampleID],
    #                  motif_activity_LPS=motif_dat2L[motif0, sampleID])#, 
               #pathway_score_2=pathscore2[pID2, sampleID])
    #cv2 <- cbind(cv2, t(rna_norm[motif00, cv2$sampleID, drop=F]))

   cv2C <- cv%>%filter(sampleID%in%shared)%>%
                mutate(motif_activity=motif_dat2C[motif0, sampleID],
                      treats="CTRL",
                      gene=motif0)#, 

    cv2L <- cv%>%filter(sampleID%in%shared)%>%
                mutate(motif_activity=motif_dat2L[motif0, sampleID],
                      treats="LPS",
                      gene=motif0)#, 

    #cv2$gene <- motif0
    df_geneall <- rbind(cv2C, cv2L)


      fname=paste(outDir, '/00.', motif0, "_combined_treat_C", cluster, "_", myvar, "_significant_gene_", motif0, ".png", sep="")
      #pdf(fname)
      png(fname, width=600, height=600, res=120)        

      #plot combined
        #geneSymbol <- grch38$symbol[grch38$ensgene == geneID]
     p3 <- ggplot(df_geneall, aes(x = get(myvar), y = motif_activity, color = treats)) +
          geom_point(position=position_jitter(w=0.1, h=0)) +
          #geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
          geom_smooth(method = lm, se = TRUE) +
            scale_color_manual(values = c("CTRL" = "#00BFC4",  
                                          "LPS"  = "#F8766D")) + 
          theme_minimal() +
          #ggtitle(geneID) +
          ggtitle(paste0(motif0, " - C", cluster, " - CTRL and LPS"))+#, 
            #subtitle = paste0("Control: r=", cor_valueC, ", p=", p_valueC, 
            # "\nLPS: r=", cor_valueL, ", p=", p_valueL)) +
            #subtitle = paste0(
             # "Control: LFC=", LFC_C, ", padj=", padj_C, 
             # "\nLPS: LFC=", LFC_L, ", padj=", padj_L)) + #, 
              #"\nInteraction: LFC=", LFC_I, ", padj=", padj_I)) +
            #subtitle = paste0("r=", cor_valueL, ", p=", p_valueL)) +
          labs(x=var_name,
               y="Normalized Motif Activity") +
          theme(plot.title = element_text(hjust=0.5, size = rel(1.3)))
        print(p3)
        dev.off() 


      }
    }

  }





















  #if(nrow(resSig) > 0 ) {
    #DEGs <- rownames(resSig)
    #stopifnot(all(DEGs %in% names(dds)))
    #if(length(DEGs) > 10) {
     # DEGs <- DEGs[DEGs == mygene]
    #}

    ########################### ctrl 
    DEGs <- mygene
    # normalized counts for the genes
    #tcountsC <- t(log2((counts(ddsC[DEGs, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    
    #tcountsC <- t(normC[DEGs, ]) #%>%
      #merge(colData(ddsC), ., by="row.names") %>%
      #gather(gene, expression, (ncol(.)-length(DEGs)+1):ncol(.))
 
      #tcounts <- setNames(as.data.frame(normcnt[DEGs, ]), DEGs) %>%
      #merge(colData(ddsC), ., by="row.names") %>% as.data.frame() %>%
      #gather(mygene, expression, (ncol(.)-length(DEGs)+1):ncol(.))

      tcountsC <- t(normC[DEGs, , drop = FALSE]) %>%
      merge(colData(ddsC), ., by = "row.names") #%>%
      #rename(expression = all_of(DEGs))
      colnames(tcountsC)[ncol(tcountsC)] <- "expression"
      tcountsC$gene <- DEGs


    fname=paste(outDir, '/', mygene, "_CTRL_C", cluster, "_", myvar, "_significant_gene_", mygene, ".png", sep="")
    #pdf(fname)
    png(fname, width=600, height=600, res=120)        

    for(j in DEGs) {
      #geneID <- substr(j, 0, 15)
      geneID <- j
        df_geneC <- filter(tcountsC, gene == j)
              # Calculate Pearson correlation

        #regressedC <- lm(df_geneC$expression ~ df_geneC$age + df_geneC$sex_alph + df_geneC$BATCH + df_geneC$PC1 + df_geneC$PC2)
        #df_geneC$residual <- residuals(regressedC)
        #df_geneC$expression <-  df_geneC$residual


        cor_testC <- cor.test(df_geneC[[myvar]], df_geneC$expression)
        cor_valueC <- round(cor_testC$estimate, 2)  # Pearson correlation
        p_valueC <- format(cor_testC$p.value, format = "f", digits = 4)  # P-value

      #geneSymbol <- grch38$symbol[grch38$ensgene == geneID]
      p1 <- ggplot(df_geneC, aes(get(myvar), expression)) +
        geom_point(position=position_jitter(w=0.1, h=0)) +
        geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
        theme_minimal() +
        #ggtitle(geneID) +
        ggtitle(paste0(geneID, " - C", cluster, " - CTRL"), subtitle = paste0("r=", cor_valueC, ", p=", p_valueC)) +
        labs(x=var_name,
             y="Normalized Gene Expression") +
        theme(plot.title = element_text(hjust=0.5, size = rel(1.3)))
      print(p1)
    }     
    dev.off() 

    ########################### LPS 
    DEGs <- mygene
    # normalized counts for the genes
    #tcountsL <- t(log2((counts(ddsL[DEGs, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
    #  merge(colData(ddsL), ., by="row.names") %>%
    #  gather(gene, expression, (ncol(.)-length(DEGs)+1):ncol(.))

      tcountsL <- t(normL[DEGs, , drop = FALSE]) %>%
      merge(colData(ddsL), ., by = "row.names") #%>%
      #rename(expression = all_of(DEGs))
      colnames(tcountsL)[ncol(tcountsL)] <- "expression"
      tcountsL$gene <- DEGs


    fname=paste(outDir, '/', mygene, "_LPS_C", cluster, "_", myvar, "_significant_gene_", mygene, ".png", sep="")
    #pdf(fname)
    png(fname, width=600, height=600, res=120)        

    for(j in DEGs) {
      #geneID <- substr(j, 0, 15)
      geneID <- j
        df_geneL <- filter(tcountsL, gene == j)
              # Calculate Pearson correlation
 
        #regressedL <- lm(df_geneL$expression ~ df_geneL$age + df_geneL$sex_alph + df_geneL$BATCH + df_geneL$PC1 + df_geneL$PC2)
        #df_geneL$residual <- residuals(regressedL)
        #df_geneL$expression <-  df_geneL$residual


        cor_testL <- cor.test(df_geneL[[myvar]], df_geneL$expression)
        cor_valueL <- round(cor_testL$estimate, 2)  # Pearson correlation
        p_valueL <- format(cor_testL$p.value, format = "f", digits = 4)  # P-value

      #geneSymbol <- grch38$symbol[grch38$ensgene == geneID]
      p2 <- ggplot(df_geneL, aes(get(myvar), expression)) +
        geom_point(position=position_jitter(w=0.1, h=0)) +
        geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
        theme_minimal() +
        #ggtitle(geneID) +
        ggtitle(paste0(geneID, " - C", cluster, " - LPS"), subtitle = paste0("r=", cor_valueL, ", p=", p_valueL)) +
        labs(x=var_name,
             y="Normalized Gene Expression") +
        theme(plot.title = element_text(hjust=0.5, size = rel(1.3)))
      	print(p2)
    }     
    dev.off() 
  #}
  
  res_dataC %>% filter(identifier == "CCL20")
  ################ combined and do it in one plot: 
  		df_geneC$treatment <- "CTRL"
  		df_geneL$treatment <- "LPS"

      #df_geneC$out[df_geneC %in% sids] <- "outlier_CTRL"
      #df_geneL$out[df_geneL %in% sids] <- "outlier_LPS"

  		LFC_C <- round(res_dataC$logFC[res_dataC$identifier==DEGs], 4)
  		padj_C <- round(res_dataC$padj[res_dataC$identifier==DEGs], 4)

  		LFC_L <- round(res_dataL$logFC[res_dataL$identifier==DEGs], 4)
  		padj_L <- round(res_dataL$padj[res_dataL$identifier==DEGs], 4)

  		#LFC_I <- round(res_dataI$logFC[res_dataI$identifier==DEGs], 4)
  		#padj_I <- round(res_dataI$padj[res_dataI$identifier==DEGs], 4)

  		df_geneall <- rbind(df_geneC, df_geneL)

	    fname=paste(outDir, '/00.', mygene, "_combined_treat_C", cluster, "_", myvar, "_significant_gene_", mygene, ".png", sep="")
	    #pdf(fname)
	    png(fname, width=600, height=600, res=120)        

  		#plot combined
	      #geneSymbol <- grch38$symbol[grch38$ensgene == geneID]
		 p3 <- ggplot(df_geneall, aes(x = get(myvar), y = expression, color = treatment)) +
	        geom_point(position=position_jitter(w=0.1, h=0)) +
	        #geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
	        geom_smooth(method = lm, se = TRUE) +
            scale_color_manual(values = c("CTRL" = "#00BFC4",  
                                          "LPS"  = "#F8766D")) + 
	        theme_minimal() +
	        #ggtitle(geneID) +
	        ggtitle(paste0(geneID, " - C", cluster, " - CTRL and LPS"), 
	        	#subtitle = paste0("Control: r=", cor_valueC, ", p=", p_valueC, 
	        	#	"\nLPS: r=", cor_valueL, ", p=", p_valueL)) +
	        	subtitle = paste0(
              "Control: LFC=", LFC_C, ", padj=", padj_C, 
	        		"\nLPS: LFC=", LFC_L, ", padj=", padj_L)) + #, 
	        		#"\nInteraction: LFC=", LFC_I, ", padj=", padj_I)) +
	        	#subtitle = paste0("r=", cor_valueL, ", p=", p_valueL)) +
	        labs(x=var_name,
	             y="Normalized Gene Expression") +
	        theme(plot.title = element_text(hjust=0.5, size = rel(1.3)))
	      print(p3)
	      dev.off() 

	      # get the top 5 highly expressed ctrl genes
	      #subdfC <- df_geneall %>% arrange(desc(expression)) %>% filter(treatment == "CTRL") %>% slice_head(n = 5)
		    #subdfL <- df_geneall %>% filter(Sample_ID %in% subdfC$Sample_ID) %>% filter(treatment == "LPS")
		    #subdfall <- rbind(subdfC, subdfL)  
          #write table
          #fname=paste(statsDir, '/00.', mygene, "_combined_treat_C", cluster, "_", myvar, "_top5CTRL_", mygene, ".txt", sep="")
          #write.table(subdfall, fname, sep="\t", row.names=F, quote=F)

            }, 
            error = function(e) {
             #Handle the error, e.g., log it or print a message
            print(paste("first mess:Error at iteration for variable ", myvar ," in cell type ", cluster, ":", e$message))
            },
            finally = {
            # Optional: code to run after each iteration, regardless of success or error
            print(paste("first mess: Completed iteration for variable ", myvar ," in cell type ", cluster))
            })
  }
 }
}


















fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

resdegs <- resall %>% filter(padj < 0.1)

clusters <- c("0", "1", "2", "3", "4")#, "5")#, "6")
celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")#, "R6 DC")
vars <- c("ISEL_Mean", "PSS_all_mean")

resdegs$celltype <- NA
resdegs$celltype[resdegs$cluster == "C0"] <- "R0 T CD4+"
resdegs$celltype[resdegs$cluster == "C1"] <- "R1 T CD8+"
resdegs$celltype[resdegs$cluster == "C2"] <- "R2 NK"
resdegs$celltype[resdegs$cluster == "C3"] <- "R3 Monocyte"
resdegs$celltype[resdegs$cluster == "C4"] <- "R4 B"
#resdegs$celltype[resdegs$cluster == "C5"] <- "R6 DC"

resdegs <- resdegs %>% filter(var %in% vars) %>% filter(!is.na(celltype))


res_sub <- resdegs %>%
  filter(var %in% vars)

opposite_degs <- res_sub %>%
  select(identifier, celltype, var, logFC, padj) %>%
  group_by(identifier, celltype) %>%
  filter(n_distinct(var) == 2) %>%  # must have both ISEL and PSS
  summarise(
    logFC_ISEL = logFC[var == "ISEL_Mean"],
    padj_ISEL  = padj[var == "ISEL_Mean"],
    logFC_PSS  = logFC[var == "PSS_all_mean"],
    padj_PSS   = padj[var == "PSS_all_mean"],
    .groups = "drop"
  ) %>%
  filter(sign(logFC_ISEL) != sign(logFC_PSS)) %>% # opposite direction
  mutate(max_abs_logFC = pmax(abs(logFC_ISEL), abs(logFC_PSS))) %>%
  arrange(desc(max_abs_logFC)) %>%
  select(-max_abs_logFC)   # optional helper column removal

dim(opposite_degs)
opptop <- head(opposite_degs, n=10) %>% pull(identifier)

fname=paste("/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/stats/", "ISEL_DEGs_opposite_sign.txt", sep="")
write.table(opposite_degs, fname, sep="\t", row.names=F, quote=F)


#resdegs <- resdegs %>% filter(var == "ISEL_Mean") %>% filter(!is.na(celltype))
#resdegs <- resdegs %>% filter(var == "PSS_all_mean") %>% filter(!is.na(celltype))

#resdegs <- resdegs %>% arrange(padj)
resdegs <- resdegs %>% arrange(desc(logFC))


psstop <- head(resdegs, n=10) %>% pull(identifier) %>% unique()


iseltop <- head(resdegs, n=10) %>% pull(identifier)

#padjtop <- head(resdegs, n=10) %>% pull(identifier)







iseltop <- head(resdegs, n=10) %>% pull(identifier) %>% unique()



head(resdegs, n=10) %>% select(identifier, padj, logFC, cluster, celltype, var)
   identifier       padj      logFC cluster    celltype       var
1      IFI44L 0.04991480 -1.0641925      C2       R2 NK ISEL_Mean
2       IFIT3 0.03154633 -0.9580942      C0   R0 T CD4+ ISEL_Mean
3      IFI44L 0.03754693 -0.9173664      C0   R0 T CD4+ ISEL_Mean
4       IFIT3 0.03349860 -0.8857745      C1   R1 T CD8+ ISEL_Mean
5       USP18 0.02527458 -0.8554276      C0   R0 T CD4+ ISEL_Mean
6       CMPK2 0.07835700 -0.8546501      C1   R1 T CD8+ ISEL_Mean
7       CMPK2 0.03909400 -0.8293099      C0   R0 T CD4+ ISEL_Mean
8        CD38 0.09365548 -0.8196971      C1   R1 T CD8+ ISEL_Mean
9      P2RY14 0.05890934 -0.8112191      C3 R3 Monocyte ISEL_Mean
10      IFIT3 0.08412461 -0.7947082      C2       R2 NK ISEL_Mean



head(resdegs, n=10) %>% select(identifier, padj, logFC, cluster, celltype, var)
   identifier        padj    logFC cluster    celltype          var
1  AC022816.1 0.001435024 1.769451      C3 R3 Monocyte PSS_all_mean
2        G0S2 0.006696371 1.284763      C3 R3 Monocyte PSS_all_mean
3     ZNF385D 0.009934604 1.236329      C3 R3 Monocyte PSS_all_mean
4      CXCL11 0.021971479 1.227355      C3 R3 Monocyte PSS_all_mean
5        EBI3 0.006547446 1.210538      C3 R3 Monocyte PSS_all_mean
6    MIR155HG 0.059426349 1.155267      C2       R2 NK PSS_all_mean
7        IDO1 0.010250821 1.129564      C3 R3 Monocyte PSS_all_mean
8     SLC39A8 0.000342077 1.078637      C3 R3 Monocyte PSS_all_mean
9        MT1F 0.001708990 1.046995      C3 R3 Monocyte PSS_all_mean
10      TNIP3 0.006044895 1.023855      C3 R3 Monocyte PSS_all_mean


head(opposite_degs, n=10)
 identifier celltype    logFC_ISEL padj_ISEL logFC_PSS padj_PSS
   <chr>      <chr>            <dbl>     <dbl>     <dbl>    <dbl>
 1 STAC       R3 Monocyte     -0.783    0.0802     1.02   0.00334
 2 IFIT3      R0 T CD4+       -0.958    0.0315     0.920  0.0379
 3 USP18      R0 T CD4+       -0.855    0.0253     0.698  0.0402
 4 CMPK2      R0 T CD4+       -0.829    0.0391     0.696  0.0701
 5 P2RY14     R3 Monocyte     -0.811    0.0589     0.581  0.0955
 6 CXCL10     R0 T CD4+       -0.753    0.0391     0.803  0.0379
 7 RSAD2      R0 T CD4+       -0.759    0.0375     0.674  0.0484
 8 IFIT1      R0 T CD4+       -0.753    0.0387     0.598  0.0852
 9 CD40       R3 Monocyte     -0.653    0.0414     0.731  0.00439
10 OAS1       R0 T CD4+       -0.681    0.0436     0.611  0.0662










    fname = paste0(normDir, "/C0_CTRL_rna.normal.no_voom.rds")
    fname = paste0(normDir, "/C0_LPS_rna.normal.no_voom.rds")
    fname = paste0(normDir, "/../C0_LPS_rna.normal.rds")

    norm <- readRDS(fname)

C0_CTRL_rna.normal.rds
C0_CTRL_rna.normal.no_voom.rds

C0_LPS_rna.normal.no_voom.rds

