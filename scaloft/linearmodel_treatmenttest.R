require(broom) 
library(plyr)   
library(tidyverse)
library(edgeR)
library(limma)
library(annotables)
library(data.table)

#load in original counts info for clusters and treats used
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.RData")
load(opfn)
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/"

# load annotation
u_eigenvec2 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt")
cv <- unique(u_eigenvec2)
##check for individuals with more than one sample
cv %>% count(Sample_ID) %>% dplyr::filter(n>1)
eigenvec2 <- u_eigenvec2 %>% dplyr::filter(SCAIP7_18 == 1 | SCAIP1_6 == 1)
eigenvec2$Sample_ID <- eigenvec2$dbgap.ID
psychvarstorun <- colnames(eigenvec2[,24:59])

combatrun="income_PCs_sex_age_and_treats_adjusted"
combatfolder=paste0(baseoutFolder,combatrun,"/")
    contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA-DEX"))

lmcoefs <- ldply(lapply(names(counts_ls),function(cluster){
    #cluster="C0"
        cat("running ", cluster, "\n")
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    df <- ldply(lapply(1:length(contrastdf$control),function(x){
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        ## Normalization of data
        # make edgeR object
        dge <- DGEList(counts=cluster_counts_t)
        #Transform counts to counts per million
        dge <- calcNormFactors(dge)
        cpm <- log2(cpm(dge)+1)
        rownames(cpm) <- rownames(cluster_counts_t)
        samples <- dim(cluster_counts_t)[2]
        design <- model.matrix(~ cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
        model <- lm(t(cpm) ~ design, data = data.frame(t(cpm), cv_d))
        tmodel <- as.data.frame(tidy(model))
        tmodel <- transform(tmodel, cluster=cluster, term=ifelse(grepl("treats", term), contrast, term))
        tmodeltreat <- subset(tmodel, term==contrast)
        return(tmodeltreat)
    }), data.frame)
    return(df)
}), data.frame)

lmcoefspha <- subset(lmcoefs, term=="PHA_vs_CTRL")
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

var="pnsi"
lmcoefs <- ldply(lapply(names(counts_ls),function(cluster){
    #cluster="C0"
        cat("running ", cluster, "\n")
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    df <- ldply(lapply(1:length(contrastdf$control),function(x){
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_metadata_var <- cluster_metadata_var[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        ## Normalization of data
        # make edgeR object
        dge <- DGEList(counts=cluster_counts_t)
        #Transform counts to counts per million
        dge <- calcNormFactors(dge)
        cpm <- log2(cpm(dge)+1)
        rownames(cpm) <- rownames(cluster_counts_t)
        samples <- dim(cluster_counts_t)[2]
        design <- model.matrix(~ as.numeric(cv_d$pnsi) + cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
        model <- lm(t(cpm) ~ design, data = data.frame(t(cpm), cv_d))
        tmodel <- as.data.frame(tidy(model))
        tmodel <- transform(tmodel, cluster=cluster, term=ifelse(grepl("treats", term), contrast, term))
        tmodel <- transform(tmodel, cluster=cluster, term=ifelse(grepl(var, term), var, term))
        tmodeltreat <- subset(tmodel, term==var)
        return(tmodeltreat)
    }), data.frame)
    return(df)
}), data.frame)

#lmcoefspha <- subset(lmcoefs, term=="PHA_vs_CTRL")
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

var="pnsi"
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
allvars <- psychvarstorun[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]
firstrunvars <- allvars[c(1:10)]
secondrunvars <- c(allvars[c(1:20)],"cage1")

lmcoefs <- ldply(lapply(names(counts_ls),function(cluster){
    #cluster="C0"
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    vars <- ldply(lapply(secondrunvars,function(var){
        cat("running ", cluster, var, "\n")
    df <- ldply(lapply(1:length(contrastdf$control),function(x){
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_metadata_var <- cluster_metadata_var[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        #cdf <- subset(plyr::count(cluster_metadata_var$Sample_ID),freq<2)
        #cluster_metadata_var <- subset(cluster_metadata_var, !Sample_ID %in% cdf$x)
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        cv_d <- transform(cv_d, treats=factor(treats))
        ## Normalization of data
        # make edgeR object
        dge <- DGEList(counts=cluster_counts_t)
        #Transform counts to counts per million
        dge <- calcNormFactors(dge)
        cpm <- log2(cpm(dge)+1)
        rownames(cpm) <- rownames(cluster_counts_t)
        samples <- dim(cluster_counts_t)[2]
        varcol=cv_d[,var]
        design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
        model <- lm(t(cpm) ~ design, data = data.frame(t(cpm), cv_d))
        tmodel <- as.data.frame(tidy(model))
        tmodel <- transform(tmodel, cluster=cluster, contrast=contrast, variable=var, term=ifelse(grepl("treats", term) & grepl(var, term), paste0(contrast,":",var), ifelse(grepl("treats", term),contrast, ifelse(grepl(var, term), var, term))))
        tmodeltreat <- subset(tmodel, term==paste0(contrast,":",var) | term==contrast | term==var)
        tmodeltreata <- ldply(lapply(unique(tmodeltreat$term),function(t){
            dft <- subset(tmodeltreat, term==t)
            names(dft)[1] <- "identifier"
            dft <- transform(dft, padj=p.adjust(p.value,method="BH"))
            fwrite(dft, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".",var,".lm_",t,".treatvarint.txt"))
            table <- data.frame(variable= var, cluster=cluster,term=t)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts_t))
    table$tested_genes <- paste(nrow(dft))
    table$DEGs_FDR <- paste(length(which(dft$padj<0.05)))
    table$DEGs_FDR_10 <- paste(length(which(dft$padj<0.1)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".",var,".lmstats_",t,".treatvarint.txt"))
        return(dft)
        }),data.frame)
        return(tmodeltreata)
    }), data.frame)
    return(df)
    }), data.frame)
    return(vars)
}), data.frame)

lmcoefspha <- subset(lmcoefs, term=="PHA_vs_CTRL" & variable==var)
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

lmcoefspha <- subset(lmcoefs, term=="PHA_vs_CTRL:pnsi")
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

lmcoefspha <- subset(lmcoefs, term=="pnsi" & contrast=="PHA_vs_CTRL")
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

    ldply(lapply(names(counts_ls),function(cl){
        ldply(lapply(1:length(contrastdf$control),function(x){
            con=contrastdf[x,]
            co=paste0(con$treatment,"_vs_",con$control)
            df <- subset(lmcoefs, cluster==cl & contrast==co)
            var=unique(df$variable)

    })