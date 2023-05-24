require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)

require(BiocParallel)


setwd("/wsu/home/hd/hd37/hd3768/piquelab/SCAIP_2022/Ali/09-ALOFT-SCAIP-batch-removed/01-ALOFT-SCAIP-pseudobulk-LPS")

# create a table to save all stats, need modification for each treatment!
table <- data.frame(matrix(ncol=18, nrow=1))
x <- c("symb", "variable", "LPS_NK_numb", "LPS_NK_uniq_ID", "LPS_NK_gene_numb", "LPS_NK_DEGs", "LPS_T_numb", "LPS_T_uniq_ID", "LPS_T_gene_numb", "LPS_T_DEGs",
        "LPS_B_numb", "LPS_B_uniq_ID", "LPS_B_gene_numb", "LPS_B_DEGs", "LPS_M_numb", "LPS_M_uniq_ID", "LPS_M_gene_numb", "LPS_M_DEGs")
colnames(table) <- x
   
fname=paste("stats_all_cell_types", "_numbers_genes_degs", "_LPS", ".txt", sep="")
write.table(table,file=fname,quote=F,sep="\t", row.names=F)


timestamp()

cores <- 8 

platePrefix <- "SCAIP_DESeq_measured_variables_4-27-22"
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## To run DESeq2 in parallel, using the
## BiocParallel library
register(MulticoreParam(cores))


## Gene counts: this is our data for anlaysis
load("/wsu/home/hd/hd37/hd3768/piquelab/SCAIP_2022/Ali/YtX_sel_37.comb.RData")
data <- YtX_sel_37

#Covariates file
#cov.file <- "ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_Wave-a-added_4-8-2021.txt"
#cov.file <- "/wsu/home/hd/hd37/hd3768/piquelab/SCAIP_2022/Ali/SCAIP_ALOFT_covariates_with_treatments_7-25-22.txt"

#Covariates with LSI variables added:
cov.file <- "/wsu/home/groups/piquelab/Shreya/SCAIP-ALOFT/covariates/SCAIP_ALOFT_covariate_all_treatments_LSI_added_8-11-22.txt"
covariates <- read.table(cov.file, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
covariates$sample.ID <- str_c(covariates$sample.ID, "-", covariates$Wave)

aloft_cov <- read.table("/wsu/home/groups/piquelab/Shreya/SCAIP-ALOFT/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_Wave-a-added_4-8-2021.txt", sep="\t", header=T, quote='"', comment="")

covariates$log2FC.IL13 <- NA
for(i in covariates$sample.ID) {
    if(i %in% aloft_cov$sample.ID)
    covariates$log2FC.IL13[covariates$sample.ID == i] <- aloft_cov$log2FC.IL13[aloft_cov$sample.ID == i]
}

# [1] "pedu"    "pincme"  "psesl"   "pnsi"    "cpwm"    "ppeqcf"  "cpeqcm"  "pcesdm"  "criskf" 
#[10] "cditsm"  "IL5_co"  "IL13_co" "IFNG_co" "baso_av" "eosi_av" "lymp_av" "mono_av" "neut_av"
#[19] "aBPFAM"  "aBPFPM"  "FEVPP"   "FVCPP"   "FFPP"    "csasgm"  "ctasfq"  "ctasev"  "cslpefc"
 
# [1] "parental education" "Parental income" "Neighborhood stress" "parental warmth" "PR Parent child confilct" "YR Parent child conflict" "parental depression" "Risky family"
# [10] "Youth depression" "Stimulated IL-5" "Stimulated IL-13" "Stimulated IFNy" "Basophils" "esosinophils" "Lymphocytes" "Monocytes" "Neutrophils"
# [19] "Peak flow AM" "Peak flow PM" "FEV1 percent predicted" "FVC precent predicted" "FEV1/FVC precent predicted" "Nightly asthma" "Asthma frequency" "Asthma severity" "Sleep efficiency"

#variables <- c("cssdhm", "cpeqcm", "csasgm", "FEVPP", "neut_av", "ctasev", "IL5_co", "log2FC.IL13", "cddstfm", "psesl")
#variable_names <- c("Sleep Duration", "Mother-child Conflict", "Nightly Asthma", "FEV %1 Predicted", "Neutrophil Fraction", "Asthma Severity", "Stimulated IL5", "IL13 GC Resistance", "Self Disclosure", "Subjective SES")
# test with 3 variables
#variables <- c("cpeqcm", "pincme", "cditsm")
#variable_names <- c("YR Parent child conflict", "Parental income", "Youth depression")

 
#variables <- c("LSI_C1", "LSI_C2", "LSI_C3", "LSI_C4", "LSI_PC1", "LSI_PC2", "LSI_PC3")
#variable_names <- c("YR Child-Family relationship", "YR School Functioning", "YR Marital Relationship", "YR Family Health", "YRPR Child-Family relationship", "YRPR School Functioning", "YRPR Parent Intimate Relationship")

# removing batches 2 and 3 from the covairate file
table(covariates$LibBatch)
#cv <- filter(cv, LibBatch != "SCAIP2")
covariates <- covariates %>% filter(!LibBatch %in% c('SCAIP2', 'SCAIP3'))


# all variables and LSI variables added
variables <- c("pedu", "pincme", "psesl", "pnsi", "cpwm", "ppeqcf", "cpeqcm", "pcesdm", "criskf", 
                "cditsm", "IL5_co", "IL13_co", "IFNG_co", "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP","csasgm", "ctasfq", "ctasev", "cslpefc",
                "cssdhm", "log2FC.IL13",  "cddstfm",
                "LSI_C1", "LSI_C2", "LSI_C3", "LSI_C4", "LSI_PC1", "LSI_PC2", "LSI_PC3"

                )
 
variable_names <- c("Parental Education", "Parental Income","Subjective SES", "Neighborhood Stress", "Parental Warmth", "PR Parent Child Confilct", "YR Parent Child Conflict", "Parental Depression", "Risky Family",
                "Youth Depression", "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", "Nightly Asthma", "Asthma Frequency", "Asthma Severity", "Sleep Efficiency",
                "Sleep Duration", "IL13 GC Resistance", "Self Disclosure",
                "YR Child-Family relationship", "YR School Functioning", "YR Marital Relationship", "YR Family Health", "YRPR Child-Family relationship", "YRPR School Functioning", "YRPR Parent Intimate Relationship"
                
                )



#variables <- c("pedu", "pincme", "psesl"
                
 #              )
 

#variable_names <- c("Parental Education", "Parental Income","Subjective SES"
                 
#              )



#writing R output to this textfile
sink(paste0(platePrefix,"_results_LPS_treatment.txt"))

timestamp()

for(i in 1:length(variables)){
    myvar <- variables[i]
    var_name <- variable_names[i]
    cat("######################################################################\n")
    cat("##Processing: ", var_name, "\n")
    cv <- covariates
    cv <- cv[!is.na(cv[, myvar]),]
    #isolate data for each cell type and make the first column match the columns in data set
    cv <- cv %>% select(bti, everything())
    cvNK <- cv[cv$MCls == "NKcell" & cv$treats == "LPS-EtOH",]
    cvT <- cv[cv$MCls == "Tcell" & cv$treats == "LPS-EtOH",]
    cvB <- cv[cv$MCls == "Bcell" & cv$treats == "LPS-EtOH",]
    cvM <- cv[cv$MCls == "Monocyte" & cv$treats == "LPS-EtOH",]
    
    #add treatment to the subset
    #lps_cvNK <- cv[cv$MCls == "NKcell" & cv$treats == "LPS-EtOH",]

    ## assign variables, load data, and load experiment information
    topDirectory <- 'out_data_'
    outDir <- paste(topDirectory, platePrefix,"/",myvar, sep='')
    system(paste("mkdir -p",outDir))
    ##
    plotsDir <- paste(outDir, '/plots', sep='')
    system(paste("mkdir -p", plotsDir))
    ##
    statsDir <- paste(outDir, '/stats', sep='')
    system(paste("mkdir -p", statsDir))
    ##
    dataDir <- paste(outDir, '/data_objects', sep='')
    system(paste("mkdir -p", dataDir))

    ##################################################################
    ## Manual conversion to R factor objects:

    cvNK$Wave <- factor(cvNK$Wave)
    NKWaveLevels <- levels(cvNK$Wave)
    cat("#", NKWaveLevels, "\n")

    cvNK$LibBatch <- factor(cvNK$LibBatch)
    NKLibBatchLevels <- levels(cvNK$LibBatch)
    cat("#", NKLibBatchLevels, "\n")

    cvNK$Sex <- factor(cvNK$Sex)
    NKSexLevels <- levels(cvNK$Sex)
    cat("#", NKSexLevels, "\n")

    cvT$Wave <- factor(cvT$Wave)
    TWaveLevels <- levels(cvT$Wave)
    cat("#", TWaveLevels, "\n")

    cvT$LibBatch <- factor(cvT$LibBatch)
    TLibBatchLevels <- levels(cvT$LibBatch)
    cat("#", TLibBatchLevels, "\n")

    cvT$Sex <- factor(cvT$Sex)
    TSexLevels <- levels(cvT$Sex)
    cat("#", TSexLevels, "\n")

    cvB$Wave <- factor(cvB$Wave)
    BWaveLevels <- levels(cvB$Wave)
    cat("#", BWaveLevels, "\n")

    cvB$LibBatch <- factor(cvB$LibBatch)
    BLibBatchLevels <- levels(cvB$LibBatch)
    cat("#", BLibBatchLevels, "\n")

    cvB$Sex <- factor(cvB$Sex)
    BSexLevels <- levels(cvB$Sex)
    cat("#", BSexLevels, "\n")

    cvM$Wave <- factor(cvM$Wave)
    MWaveLevels <- levels(cvM$Wave)
    cat("#", MWaveLevels, "\n")

    cvM$LibBatch <- factor(cvM$LibBatch)
    MLibBatchLevels <- levels(cvM$LibBatch)
    cat("#", MLibBatchLevels, "\n")

    cvM$Sex <- factor(cvM$Sex)
    MSexLevels <- levels(cvM$Sex)
    cat("#", MSexLevels, "\n")
    ##################################################################
    ## Preparing data for DEseq:
    ## Combine processed data into a DESeqDataSet
    ## & filter to keep genes expressed in 50% of samples at least

    cat("Now testing ", var_name, " for NK Cells \n")
    ##NK Cells
    dataNK <- data[,colnames(data) %in% cvNK$bti]
    n.barcodes <- dim(dataNK)[2]
    cvNK<- (cvNK[which(cvNK$bti %in% colnames(dataNK)),]) #So cv and data match
    cvNK <- cvNK[order(cvNK$bti),]
    dataNK <- dataNK[,order(colnames(dataNK))]
    colnames(dataNK)==cvNK$bti #Check that files are in order
    cat("Number of NK cells: ", nrow(cvNK)," \n")
    cat("Number of individuals - NK cells: ", length(unique(cvNK$sample.ID)), "\n")


    filtered_dataNK <- dataNK > 0
    table(rowSums(filtered_dataNK))
    keepNK <- rowSums(filtered_dataNK) >= (ncol(dataNK) / 4)
    dataNK <- dataNK[keepNK, ]
    cat("# of genes tested for NK cells - ", myvar, " :", nrow(dataNK)," \n")

    ddsNK <- DESeqDataSetFromMatrix(
        countData = round(dataNK),
        colData = cvNK,
        design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    colnames(ddsNK) <- cvNK$bti
    ## Fit the model on the whole plate
    system.time(ddsNK <- DESeq(ddsNK,parallel=TRUE))

    ##save
    save(ddsNK, file=paste(dataDir, '/DESeq2_', platePrefix,"_",myvar, '_NKcell.RData', sep=''))
    resNK <- results(ddsNK, parallel=TRUE)
    summary(resNK)
    resNK <- as.data.frame(resNK)
    fname=paste(statsDir, '/', platePrefix,"_",myvar, "_NKcells_", "DESeq_results", ".txt", sep="")
    write.table(resNK,file=fname,quote=F,sep="\t")
    resNK_sig <- resNK %>% filter(padj <0.1)
    resNK_sig <- resNK_sig[order(resNK_sig$padj),]
    cat("# of DEGs for NK cells -", myvar, "(FDR < 0.1) :", nrow(resNK_sig),"\n")
    fname=paste(statsDir, '/', platePrefix,"_",myvar,"_NKcells_", "significant_DESeq_results", ".txt", sep="")
    write.table(resNK_sig,file=fname,quote=F,sep="\t")


    ##T Cells
    cat("Now testing ", var_name, " for T Cells \n")
    dataT <- data[,colnames(data) %in% cvT$bti]
    n.barcodes <- dim(dataT)[2]
    cvT<- (cvT[which(cvT$bti %in% colnames(dataT)),]) #So cv and data match
    cvT <- cvT[order(cvT$bti),]
    dataT <- dataT[,order(colnames(dataT))]
    colnames(dataT)==cvT$bti #Check that files are in order
    cat("Number of T cells: ", nrow(cvT)," \n")
    cat("Number of individuals - T cells: ", length(unique(cvT$sample.ID)), "\n")

    filtered_dataT <- dataT > 0
    table(rowSums(filtered_dataT))
    keepT <- rowSums(filtered_dataT) >= (ncol(dataT) / 4)
    dataT <- dataT[keepT, ]
    cat("# of genes tested for T cells - ", myvar, " :", nrow(dataT)," \n")


    ddsT <- DESeqDataSetFromMatrix(
        countData = round(dataT),
        colData = cvT,
        design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    colnames(ddsT) <- cvT$bti

    ## Fit the model on the whole plate
    system.time(ddsT <- DESeq(ddsT,parallel=TRUE))

    ##save
    save(ddsT, file=paste(dataDir, '/DESeq2_', platePrefix,"_",myvar, '_Tcell.RData', sep=''))
    resT <- results(ddsT, parallel=TRUE)
    summary(resT)
    resT <- as.data.frame(resT)
    fname=paste(statsDir, '/', platePrefix,"_",myvar, "_Tcells_", "DESeq_results", ".txt", sep="")
    write.table(resT,file=fname,quote=F,sep="\t")
    resT_sig <- resT %>% filter(padj <0.1)
    resT_sig <- resT_sig[order(resT_sig$padj), ]
    cat("# of DEGs for T cells - ", myvar, " (FDR < 0.1) :", nrow(resT_sig)," \n")
    fname=paste(statsDir, '/', platePrefix,"_",myvar, "_Tcells_", "significant_DESeq_results", ".txt", sep="")
    write.table(resT_sig,file=fname,quote=F,sep="\t")

    ##B Cells
    cat("Now testing", var_name, "for B Cells\n")
    dataB <- data[,colnames(data) %in% cvB$bti]
    n.barcodes <- dim(dataB)[2]
    cvB<- (cvB[which(cvB$bti %in% colnames(dataB)),]) #So cv and data match
    cvB <- cvB[order(cvB$bti),]
    dataB <- dataB[,order(colnames(dataB))]
    colnames(dataB)==cvB$bti #Check that files are in order
    cat("Number of B cells: ", nrow(cvB)," \n")
    cat("Number of individuals - B cells: ", length(unique(cvB$sample.ID)), "\n")

    filtered_dataB <- dataB > 0
    table(rowSums(filtered_dataB))
    keepB <- rowSums(filtered_dataB) >= (ncol(dataB) / 4)
    dataB <- dataB[keepB, ]
    cat("# of genes tested for B cells - ", myvar, " :", nrow(dataB)," \n")


    ddsB <- DESeqDataSetFromMatrix(
        countData = round(dataB),
        colData = cvB,
        design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    colnames(ddsB) <- cvB$bti

    ## Fit the model on the whole plate
    system.time(ddsB <- DESeq(ddsB,parallel=TRUE))

    ##save
    save(ddsB, file=paste(dataDir, '/DESeq2_', platePrefix, "_",myvar,'_Bcell.RData', sep=''))
    resB <- results(ddsB, parallel=TRUE)
    summary(resB)
    resB <- as.data.frame(resB)
    fname=paste(statsDir, '/', platePrefix,"_",myvar, "_Bcells_", "DESeq_results", ".txt", sep="")
    write.table(resB,file=fname,quote=F,sep="\t")
    resB_sig <- resB %>% filter(padj <0.1)
    resB_sig <- resB_sig[order(resB_sig$padj), ]
    cat("# of DEGs for B cells - ", myvar, " (FDR < 0.1) :", nrow(resB_sig)," \n")
    fname=paste(statsDir, '/', platePrefix,"_",myvar, "_Bcells_", "significant_DESeq_results", ".txt", sep="")
    write.table(resB_sig,file=fname,quote=F,sep="\t")

    ##Monocytes
    cat("Now testing ", var_name, " for Monocytes\n")
    dataM <- data[,colnames(data) %in% cvM$bti]
    n.barcodes <- dim(dataM)[2]
    cvM<- (cvM[which(cvM$bti %in% colnames(dataM)),]) #So cv and data match
    cvM <- cvM[order(cvM$bti),]
    dataM <- dataM[,order(colnames(dataM))]
    colnames(dataM)==cvM$bti #Check that files are in order
    cat("Number of Monocytes: ", nrow(cvM)," \n")
    cat("Number of individuals - Monocytes: ", length(unique(cvM$sample.ID)), "\n")

    filtered_dataM <- dataM > 0
    table(rowSums(filtered_dataM))
    keepM <- rowSums(filtered_dataM) >= (ncol(dataM) / 4)
    dataM <- dataM[keepM, ]
    cat("# of genes tested for Monocytes - ", myvar, " :", nrow(dataM)," \n")

    ddsM <- DESeqDataSetFromMatrix(
        countData = round(dataM),
        colData = cvM,
        design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    colnames(ddsM) <- cvM$bti

    ## Fit the model on the whole plate
    system.time(ddsM <- DESeq(ddsM,parallel=TRUE))

    ##save
    save(ddsM, file=paste(dataDir, '/DESeq2_', platePrefix, "_",myvar,'_Monocyte.RData', sep=''))
    resM <- results(ddsM, parallel=TRUE)
    summary(resM)
    resM <- as.data.frame(resM)
    fname=paste(statsDir, '/', platePrefix, "_",myvar,"_Monocytes_", "DESeq_results", ".txt", sep="")
    write.table(resM,file=fname,quote=F,sep="\t")
    resM_sig <- resM %>% filter(padj <0.1)
    resM_sig <- resM_sig[order(resM_sig$padj), ]
    cat("# of DEGs for Monocytes - ", myvar, " (FDR < 0.1) :", nrow(resM_sig)," \n")
    fname=paste(statsDir, '/', platePrefix, "_",myvar,"_Monocytes_", "significant_DESeq_results", ".txt", sep="")
    write.table(resM_sig,file=fname,quote=F,sep="\t")

    

    table <- data.frame(matrix(ncol=18, nrow=1))
    x <- c("symb", "variable", "LPS_NK_numb", "LPS_NK_uniq_ID", "LPS_NK_gene_numb", "LPS_NK_DEGs", "LPS_T_numb", "LPS_T_uniq_ID", "LPS_T_gene_numb", "LPS_T_DEGs",
        "LPS_B_numb", "LPS_B_uniq_ID", "LPS_B_gene_numb", "LPS_B_DEGs", "LPS_M_numb", "LPS_M_uniq_ID", "LPS_M_gene_numb", "LPS_M_DEGs")
    colnames(table) <- x
    table$symb <- paste(myvar)
    table$variable <- paste(var_name)

    table$LPS_NK_numb <- paste(nrow(cvNK))
    table$LPS_NK_uniq_ID <- paste(length(unique(cvNK$sample.ID)))
    table$LPS_NK_gene_numb <- paste(nrow(dataNK))
    table$LPS_NK_DEGs <- paste(nrow(resNK_sig))

    table$LPS_T_numb <- paste(nrow(cvT))
    table$LPS_T_uniq_ID <- paste(length(unique(cvT$sample.ID)))
    table$LPS_T_gene_numb <- paste(nrow(dataT))
    table$LPS_T_DEGs <- paste(nrow(resT_sig))

    table$LPS_B_numb <- paste(nrow(cvB))
    table$LPS_B_uniq_ID <- paste(length(unique(cvB$sample.ID)))
    table$LPS_B_gene_numb <- paste(nrow(dataB))
    table$LPS_B_DEGs <- paste(nrow(resB_sig))

    table$LPS_M_numb <- paste(nrow(cvM))
    table$LPS_M_uniq_ID <- paste(length(unique(cvM$sample.ID)))
    table$LPS_M_gene_numb <- paste(nrow(dataM))
    table$LPS_M_DEGs <- paste(nrow(resM_sig))

    # save it as a csv file

    #fname=paste(statsDir, '/', platePrefix,"_",myvar,"_all_cell_types_stats", "_numbers_genes_degs", "_LPS", ".txt", sep="")
    fname=paste("stats_all_cell_types", "_numbers_genes_degs", "_LPS", ".txt", sep="")

    #write.table(table,file=fname,quote=F,sep="\t", row.names=F)
    write.table(table, file = fname, sep = "\t", append = T, quote = F, col.names = F, row.names = F)


}


