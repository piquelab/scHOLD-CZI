require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
library(annotables)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(plyr); library(dplyr)
library(parallel)
library("AnnotationHub")
library(ggseurat)
library(cowplot)
library(sva)
library(ggpubr)
library(flextable)
library(ftExtra)
library(rlist)
library(officer)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 15 * 1024 ^ 3)

# running treatment DEGs
zingeRrun=FALSE
job="CZI"

if(job=="ALOFT"){
cat("need to input")
combat="withCOMBAT"
run=paste0("treatment_",combat)
contrastdf <- data.frame(control=c("CTRL"),treatment=c("LPS"))
} else if (job=="CZI"){
    args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.15) #for testing
    base <- args[1]
    project=args[3]
    method=args[4]
    dimset=args[5]
    resset=args[6]
    baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
    combat="noCOMBAT"
    run=paste0("treatment_",combat)
    #combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
    opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
    load(opfn)
    contrastdf <- data.frame(control=c("RNA-CTRL"),treatment=c("RNA-LPS"))
}

outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

contrastdf <- transform(contrastdf, contrast=paste0(treatment,"_vs_",control))

lapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
#    if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt")) > 0)){
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    if(combat=="withCOMBAT"){
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
        if(job=="ALOFT"){
            adjusted_counts <- adjusted
        }
    } else {
        adjusted_counts <- cluster_counts
    }

    if(job=="ALOFT"){
    contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA-DEX"))
    } else if(job=="CZI"){
    contrastdf <- data.frame(control=c("RNA-CTRL"),treatment=c("RNA-LPS"))
    }
    mclapply(1:nrow(contrastdf),function(x) {
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cat("running",contrast,"\n")
        if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt")) > 0)){
        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","treats")]
        cluster_metadata_var <- subset(cluster_metadata_var, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        indcount <- plyr::count(cluster_metadata_var,"Sample_ID")
        cluster_metadata_var <- subset(cluster_metadata_var, Sample_ID %in% subset(indcount, freq>1)$Sample_ID)
        design <-  paste0("~ Sample_ID + treats")
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        genestoremove <- data.frame(genes=unique(rownames(which(cluster_counts_t >= .Machine$integer.max, arr.ind = TRUE))))
        cluster_counts_t <- cluster_counts_t[!rownames(cluster_counts_t) %in% genestoremove$genes,]

        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        if(zingeRrun){
        designb = model.matrix( as.formula(design) , cluster_metadata_var)
        weights <- zeroWeightsLS(counts=cluster_counts_t, design=designb, maxit=200, normalization="DESeq2_poscounts", colData=cluster_metadata_var, designFormula=as.formula(design))
        assays(dds)[["weights"]]=weights
        dds = DESeq2::estimateSizeFactors(dds, type="poscounts")
        dds = estimateDispersions(dds)
        dds = nbinomWaldTest(dds, betaPrior=FALSE, useT=TRUE, df=rowSums(weights)-2) #betaPrior=FALSE should be used for designs with interactions >> note: orignial script is betaPrior=TRUE
        } else {
        dds <- DESeq(dds,parallel=TRUE) #was parallel=TRUE causing issue?
        }
        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",cluster,".",contrast,".",run,".RData")
        save(dds, file=opfn)

        res <- results(dds, contrast=c("treats",con$treatment,con$control))
        sub.table <- data.frame(res@rownames, res$'baseMean', res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier','baseMean','padj','pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$cluster=cluster
        sub.table$contrast=contrast
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres.",contrast,".",run,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        table <- data.frame(cluster=cluster, contrast=contrast)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$tested_genes <- paste(nrow(sub.table))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt"))
        }
        })    
    #}
})

for (x in c(1:nrow(contrastdf))){
    con=contrastdf[x,]
    contrast=con$contrast
    cat("running ",contrast,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres.",contrast,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres.",contrast,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types.",contrast,".",run,".txt"))
}

subvars <- ldply(lapply(1:nrow(contrastdf), function(x){
    con=contrastdf[x,]
    contrast=con$contrast
    cat("running",contrast,"\n")
    if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types.",contrast,".",run,".txt"))){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types.",contrast,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newtable <- data.frame(contrast=contrast,list.cbind(newls))
            return(newtable)
        }
    }
}),data.frame)

dft <- subvars %>% flextable() %>% span_header()  

dft <- align(dft, i = 1:2, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=max)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header() %>% set_caption(caption = run) 
cols <- seq(1,length(subsubvars),by=3)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1:2, j = c(2:length(subsubvars)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(subsubvars)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(cols), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("numind|ngenes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() %>% split_header()
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = c(1), border = border, part = "all")
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(subsubvars_degonly)), align = "center", part = "body")

#dft <- align(dft, i = 1, j = NULL, align = "left", part = "body")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".ndegonly.png"))

library(EnhancedVolcano)
deseqres <- fread(paste0(baseoutFolder,run,"/deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
for (c in names(counts_ls)){
    cat("running ", c, "\n")
    sub.table <- subset(deseqres, cluster==c)
    sub.table <- sub.table[order(sub.table$padj,-abs(sub.table$logFC)), ]
    xmax=max(sub.table$logFC)+0.5
    xmin=min(sub.table$logFC)-0.5
    topp1 <- min(subset(sub.table, padj > quantile(padj, prob = 1 -99/100,na.rm=T))$pvalue)
    tolab <- c(unique(head(sub.table,n=20)$identifier))
    png(width = 8, height = 8, file=paste0(outFolder,"figures/dge_volcano-",run,"-",c,".png"), pointsize=12, 
          bg = "transparent", units = "in", res = 1200)
    p <- EnhancedVolcano(sub.table,
      lab = sub.table$identifier,
      x = 'logFC',
      y = 'pvalue',
      title=paste0(c),
      subtitle = NULL,
      xlim = c(xmin, xmax),
      pCutoffCol= 'padj',
      pCutoff = 0.1,
      FCcutoff = 0.5,
      labSize = 3.0,
      #pointSize = c(ifelse(sub.table$gene_symbol %in% tolab, 5,1.5)),
      col=c('black', 'black', 'blue', 'red3'),
      hline = c(topp1),
      hlineCol = c('green'),
      selectLab = tolab,
      legendPosition = 'none',
      drawConnectors = TRUE)
    print(p)
    dev.off()
}
