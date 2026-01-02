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

dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

#normDir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data/alltreat"
normDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/norm_counts/alltreats"

timestamp()


platePrefix <- "ALL.0.1.13." 
plateSufixC <- "-RNA-CTRL.noCombat_DESeq" 
plateSufixL <- "-RNA-LPS.noCombat_DESeq" 


#countSufix <- ".SES_PCs_sex_age_and_treats_adjusted"
deseqSuffix <- ".noCombat_DESeq"
deseqDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/"


variables <- c("PSS_all_mean", "ISEL_Mean")#, "cytocomp")
variable_names <- c("Perceived Stress",  "Social Support")#, "Cytokine Composite")

# top CTRL genes
#mygenes <- c("STAT1", "STAT2", "STAT3", "JAK1", "JAK2", "IRF2", "IRF9", "INFA1", "IFNB1", "IFNGR1", "IFNGR2", "IFIT1", "IFITM2", "IFNA21", "IFI44", "ISG15")
#mygenes <- c("STAT2", "STAT1", "LIPA", "ISG15", "PARP12", "IFIT3", "TAP1", "IFI6", "ITGAX", "ICAM1", "CCL4L2", "DDX58")
#mygenes <- c("ITGAX", "ICAM1", "CCL4L2")
#mygenes <- c("DDX58")

#top LPS DEGs
# PSS C0: CXCL11 and TUBB8, C3: RUNDC1
# ISEL C0 CXCL11, C2: TNFRSF1B C12orf75 MYO1F, C4: ITGAM AC003991.1 TNFRSF18
mygenes <- c("CXCL11", "TUBB8", "RUNDC1", "TNFRSF1B", "C12orf75", "MYO1F", "ITGAM", "AC003991.1", "TNFRSF18")


#top treat variable interactions
# PSS C4: ICAM1. more C4: CCL3L1, ICAM1, PPM1L, C5: IGLV6-57
# ISEL C0: RSAD2, USP18, C2: STAT1, ISG15, C3: IFIT3, C4: HCAR2
#mygenes <- c("ICAM1", "RSAD2", "USP18", "STAT1", "ISG15", "IFIT3", "HCAR2", "CCL3L1", "ICAM1", "TNIP3", "PPM1L")


# intereactons only
# PSS C3: CCL20, C4: RAD51B, JUNB, CCL2, WTAP
# isel C0: RSAD2, USP18. C1: CCL3. C4: USP18, PTX3, IL27, BATF3
mygenes <- c("CCL20", "RAD51B", "JUNB", "CCL2", "WTAP", "RSAD2", "USP18", "CCL3", "IL27", "BATF3")

# intereactons + CTRL
# PSS C4: SOD2, CXCL3, CD44, IL27, ICAM1, ISG20, IRAK2
# isel C2: PARP9, STAT1, IFI1, ISG15, IRF7 C3: STAT1, ISG15, IFI6, DDX58, GBP1. C4: ISG15, CD38, ISG20, IFIT2, CD40



mygenes <- c("CXCL11", "TUBB8")


mygenes <- c("CCL20")


mygenes <- c("SOD2", "CXCL3", "CD44", "IL27", "ICAM1", "ISG20", "IRAK2", "PARP9", "STAT1", "ISG15", "IRF7", "CD40", "DDX58")


mygenes <- c("CXCL11", "RUNDC1", "TNFRSF18", "BAZ2B", "ARHGAP10", "FOXO3", "SYTL3", "AC003991.1", "AC024382.1", "AL451123.1",
            "LINC02732", "FGD4", "ITGAM", "AC004540.2", "TRAV6", "BID")


#top 10 lowest pvalues, PSS C3
mygenes <- padjtop

mygenes <- iseltop
mygenes <- psstop


mygenes <- opptop


fname="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/barplot/12.3.all_LPS_ressigs.txt"
lpssig <- read.table(fname, sep="\t", header=T, quote='"', comment="")

mygenes <- c("CXCL11", "RUNDC1", "TNFRSF18", "BAZ2B", "ARHGAP10", "FOXO3", "SYTL3", "AC003991.1", "AC024382.1", "AL451123.1",
            "LINC02732", "FGD4", "ITGAM", "AC004540.2", "TRAV6", "BID")


mygenes <- c("CXCL10", "RSAD2", "CMPK2", "IRF4", "CXCL11", "MT1F", "P2RY14")
mygenes <- opptop


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
  for (cluster in c(0,3)){
    #tryCatch({

    cat("######## loading CTRL for cluster ", cluster, "\n")
    # Read in the data CTRL
    fname <- paste(dataDir, platePrefix, "C", cluster, ".deseqres_", myvar, plateSufixC, ".txt", sep = "")
    #fname <- paste(dataDir, platePrefix, "C", clust, ".deseqres_", myvar, plateSufix, ".txt", sep = "")

    res_dataC <- read.table(fname, sep="\t", header=T, quote='"', comment="")
    resSigC <- res_dataC %>% filter(padj < 0.1)
    rownames(resSigC) <- resSigC$identifier
    resSigC <- resSigC[order(resSigC$padj), ]

    # Store each data frame in the list with the name "resX"
    dfvarC[[paste0(cluster)]] <- res_dataC
    dfsigC[[paste0(cluster)]] <- res_dataC %>% filter(padj < 0.1)

    #cat("## numb DEGs present in cluster", cluster, "for", var_name,":", nrow(dfsigC[[paste0(cluster)]]),  "\n")

    # read count data
    #ALL.ComBat_seq.C1.SES_PCs_sex_age_and_treats_adjusted.RData
    #cname <- paste(dataDir, platePrefix, ".ComBat_seq", ".C", cluster, countSufix, ".RData", sep="")
    #load(cname)
    #count_data <- adjusted_counts

    # read DESeq object for each cluster
    #deseqvar <- readRDS("ALL.DESeq_output-RNA-CTRL-iselC0.SES_PCs_sex_age_and_treats_adjusted.RDS")
    #ALL.0.15.13.DESeq_output-RNA-CTRL-PSS_all_meanC3.noCombat_DESeq.RData
    dname <- paste(deseqDir, platePrefix, "DESeq_output-RNA-CTRL-",myvar, "C", cluster, deseqSuffix, ".RData", sep="")
    load(dname)
    ddsC <- dds
    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS
    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS    
    cat("## total genes present in cluster", cluster, "for", var_name,":", nrow(dfvarC[[paste0(cluster)]]),  "\n")
    cat("## numb DEGs present in cluster", cluster, "for", var_name,":", nrow(dfsigC[[paste0(cluster)]]),  "\n")
    #cat("## number of genes for ", cluster, " count data:", nrow(counts),  "\n")
    cat("## dimention of DEseq count ", cluster, "for", myvar, dim(counts(ddsC)),  "\n")

    # load normalized count data
    #fname = paste0(normDir, "/C", cluster, "_CTRL_rna.normal.rds")
    #fname = paste0(normDir, "/C", cluster, "_CTRL_rna.normal.no_voom.rds")
    fname = paste0(normDir, "/C", cluster, "_rna.normal.rds")
    norm <- readRDS(fname)

    normC <- norm[, grep("CTRL", colnames(norm))]
    normL <- norm[, grep("LPS", colnames(norm))]

    #normC <- readRDS(fname)

    #C0_CTRL_rna.normal.rds
    cat("################# CTRL LOADED ##################\n")

    cat("######## loading LPS for cluster ", cluster, "\n")
    # Read in the data LPS
    fname <- paste(dataDir, platePrefix, "C", cluster, ".deseqres_", myvar, plateSufixL, ".txt", sep = "")
    
    res_dataL <- read.table(fname, sep="\t", header=T, quote='"', comment="")
    resSigL <- res_dataL %>% filter(padj < 0.1)
    rownames(resSigL) <- resSigL$identifier
    resSigL <- resSigL[order(resSigL$padj), ]

    # Store each data frame in the list with the name "resX"
    dfvarL[[paste0(cluster)]] <- res_dataL
    dfsigL[[paste0(cluster)]] <- res_dataL %>% filter(padj < 0.1)

    # read count data
    #ALL.ComBat_seq.C1.SES_PCs_sex_age_and_treats_adjusted.RData
    #cname <- paste(dataDir, platePrefix, ".ComBat_seq", ".C", cluster, countSufix, ".RData", sep="")
    #load(cname)
    #count_data <- adjusted_counts

    # read DESeq object for each cluster
    #deseqvar <- readRDS("ALL.DESeq_output-RNA-CTRL-iselC0.SES_PCs_sex_age_and_treats_adjusted.RDS")
    dname <- paste(deseqDir, platePrefix, "DESeq_output-RNA-LPS-",myvar, "C", cluster, deseqSuffix, ".RData", sep="")
    load(dname)
    ddsL <- dds

    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS
    #/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.DESeq_output-RNA-CTRL-pr_compC0.SES_PCs_sex_age_and_treats_adjusted.RDS    
    cat("## total genes present in cluster", cluster, "for", var_name,":", nrow(dfvarL[[paste0(cluster)]]),  "\n")
    cat("## numb DEGs present in cluster", cluster, "for", var_name,":", nrow(dfsigL[[paste0(cluster)]]),  "\n")
    cat("## number of genes for ", cluster, " count data:", nrow(counts),  "\n")
    cat("## dimention of DEseq count ", cluster, "for", myvar, dim(counts(ddsL)),  "\n")

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
          tryCatch({
        mygene <- mygenes[g]
         cat("gene ", mygene, "\n")

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
            # "\nLPS: r=", cor_valueL, ", p=", p_valueL)) +
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













############################### get top opposite ###########################################################

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
