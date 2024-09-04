require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(tidyverse)
library(data.table)
library(plyr)
library(cowplot)
library(parallel)
library(future)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_converted_dbgapID_key_n183_AR.txt","ALL","fastdemux") #for testing
base <- args[1]

cov_file=fread(args[2]) #this is the psych cov file

#if(!is.na(args[3])){
#removeind <- fread(paste0(base,args[3]),header=F)$V1
#project <- sapply(strsplit(args[3],"_"),function(y)y[1])
#}
project=args[3]
method=args[4]

outFolder=paste0(base,method,"_pseudobulk_ctrl/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 30 * 1024 ^ 3)

opfn <- paste0(outFolder,project,".DESeq_countlists.RData")
load(opfn)

for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[7]
        cat("running ", cluster, "\n")
    var="treat_SES_int"
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    cluster_metadata <- transform(cluster_metadata, sex_alph=as.factor(sex_alph))
    cluster_metadata <- within(cluster_metadata, sex_alph <- relevel(sex_alph, ref = "Male"))

    cluster_metadata_var <- cluster_metadata[,c("Sample_ID","treats","PC1","PC2","BATCH","sex_alph","age","SES")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

    cluster_counts <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts) == rownames(cluster_metadata_var))

    cat("running treatment deseq \n")
    # Create DESeq2 object        
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata_var, 
                                  design = ~ treats + BATCH + age + PC1 + PC2 + sex_alph + SES + treats*SES)
    #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".DESeq_output-",var,cluster,".RData")
    save(dds, file=opfn)

    res <- results(dds)
    sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
    names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
    sub.table <- sub.table[!is.na(sub.table$padj), ]
    cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    sub.table$var=var
    sub.table$cluster=cluster
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,".txt"))
    fdr=0.05
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,".txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_var))
    table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".txt"))
}

for (var in c("treat_SES_int")){
    cat("running ",var,"\n")
sub.table <- ldply(lapply(names(counts_ls), function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,".txt"),header=T)
    }),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,".txt"))
sub.table <- ldply(lapply(names(counts_ls), function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,".txt"),header=T)
}),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,".txt"))
sub.table <- ldply(lapply(names(counts_ls), function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".txt"),header=T)
}),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,".txt"))
}

