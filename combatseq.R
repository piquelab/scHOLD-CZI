library(sva)
require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(tidyverse)
library(data.table)
library(plyr)
library(cowplot)
library(parallel)
library(ggpubr)
future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)


args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_converted_dbgapID_key_n183_AR_18newvar_updated_08-05.txt","ALL","fastdemux") #for testing
base <- args[1]

cov_file=fread(args[2]) #this is the psych cov file

#if(!is.na(args[3])){
#removeind <- fread(paste0(base,args[3]),header=F)$V1
#project <- sapply(strsplit(args[3],"_"),function(y)y[1])
#}
project=args[3]
method=args[4]

outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip")
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]

opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".DESeq_countlists.RData")
load(opfn)

for (cluster in names(counts_ls)){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))

	adjusted <- ComBat_seq(cluster_counts, batch=cluster_metadata$BATCH, group=cluster_metadata$treats)
    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".adjusted.RData")
    save(adjusted, file=opfn)

    cluster_metadata_var <- cluster_metadata[,c("PC1","PC2","sex_alph","age","treats")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

    adjusted <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
    all(colnames(adjusted) == rownames(cluster_metadata_var))

    cat("running treatment deseq \n")
    # Create DESeq2 object        
    dds <- DESeqDataSetFromMatrix(adjusted, 
                                  colData = cluster_metadata_var, 
                                  design = ~ PC1 + PC2 + sex_alph + age + treats)
    #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".DESeq_output-treats",cluster,".adjusted.RData")
    save(dds, file=opfn)

    var="LPSvsCTRL"
    cluster_metadata_t <- subset(cluster_metadata, treats %in% c("RNA-LPS","RNA-CTRL"))
    res <- results(dds, contrast=c("treats","RNA-LPS","RNA-CTRL"))
    sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
    names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
    sub.table <- sub.table[!is.na(sub.table$padj), ]
    cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    sub.table$var=var
    sub.table$cluster=cluster
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,".adjusted.txt"))
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,".adjusted.txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_t))
    table$number_individuals <- paste(length(unique(cluster_metadata_t$Sample_ID)))
    table$gene_number <- paste(nrow(adjusted))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_logfc2 <- paste(nrow(subset(sub.table,padj<fdr & abs(logFC>0.5))))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".adjusted.txt"))
}

var="LPSvsCTRL"
clustersrun <- ldply(lapply(names(counts_ls), function(cluster){
    if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,".adjusted.txt"))){
        return(cluster)
    }
}),data.frame)
sub.table <- ldply(lapply(clustersrun$X..i.., function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,".adjusted.txt"),header=T)
    }),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,".adjusted.txt"))
sub.table <- ldply(lapply(clustersrun$X..i.., function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,".adjusted.txt"),header=T)
}),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,".adjusted.txt"))
sub.table <- ldply(lapply(clustersrun$X..i.., function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".adjusted.txt"),header=T)
}),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,".adjusted.txt"))

filenames <- list.files(outFolder) #file list from directory
filenames1 <- filenames[grep("DESeq_output-", filenames)] #pick specific files from list
x <- na.omit(sapply(strsplit(filenames1,"(?<=.)(?=C[0-9])",perl=TRUE),function(y)y[2]))
clusters <- unique(sapply(strsplit(x,"[.]"),function(z)z[1]))

for (cluster in clusters){
    #cluster=clusters[1]
    cat("running ", cluster, "\n")
    var="LPSvsCTRL"
    load(paste0(outFolder,project,".DESeq_output-treats",cluster,".adjusted.RData"))

    vsd <- vst(dds, blind=FALSE)
    # Plot PCA
    fig0 <- DESeq2::plotPCA(vsd, intgroup = c("treats"))+ggtitle("treats")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    #fig1 <- DESeq2::plotPCA(vsd, intgroup = c("PC1"))+ggtitle("PC1")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    #fig2 <- DESeq2::plotPCA(vsd, intgroup = c("PC2"))+ggtitle("PC2")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig3 <- DESeq2::plotPCA(vsd, intgroup = c("sex_alph"))+ggtitle("sex")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig4 <- DESeq2::plotPCA(vsd, intgroup = c("age"))+ggtitle("age")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig <-plot_grid(fig0,fig3,fig4, nrow=5, ncol=1, align="h")
    png(paste0(figuredir,project,".deseqPCA-treats","-",cluster,".adjusted.png"), width=1000, height=1000, res=120)
    print(fig)
    dev.off()
}
cluster=clusters[1]
for(cov in c("treats","sex_alph","age")){
	cat("running ",cov, "\n")
    fig <- DESeq2::plotPCA(vsd, intgroup = c(cov))+ggtitle(cov)+theme(plot.margin=margin(b=-0.8,unit="cm"))
    png(paste0(figuredir,project,".deseqPCA-treats","-",cluster,cov,".adjusted.png"), width=1000, height=1000, res=120)
    print(fig)
    dev.off()
}



for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".adjusted.RData")
    load(opfn)

            #run covariates separately for each treatment condition
        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        mclapply(psychvarstorun,function(var){
        #var="SES"
        #var <- colnames(cluster_metadata_t[con])
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(nrow(cluster_metadata_var))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".adjusted.txt"))
        })
    }
}

#this is separate as CVD risk was calculated taking into account sex and age so need those out of hte model
for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    var="CVDRISK"

    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","treats","PC1","PC2","sex_alph","age",var)]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

    adjusted <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age",var)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".",var,"_sex_age_and_treats_adjusted.RData")
    save(adjusted, file=opfn)

            #run covariates separately for each treatment condition
        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        #var <- colnames(cluster_metadata_t[con])
        design <-  paste0("~ PC1 + PC2 + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".",var,"_sex_age_and_treats_adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(nrow(cluster_metadata_var))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
    }
}

#this is separate to have CRP as a factor (known way to hande the variable) and remove individuals >10
for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    var="HS_CRP"
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","treats","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age",var)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".",var,"_sex_age_and_treats_adjusted.RData")
    save(adjusted, file=opfn)

            #run covariates separately for each treatment condition
        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- subset(cluster_metadata_t, HS_CRP<10) #advised to remove as likely an infection
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph), factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))
        cluster_metadata_t <- within(cluster_metadata_t, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))

        var="factor_HS_CRP" #we are factorizing the variable here so saving separately 
        #var <- colnames(cluster_metadata_t[con])
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".",var,"_sex_age_and_treats_adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
    }
}



cluster="C7"
for (var in c("sex","age","sex_age_int",psychvarstorun)){
    #for (var in c("SES"))
    for (i in unique(metadata_ls[[1]]$treats)){
    	        if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".adjusted.txt"))){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".adjusted.txt"),header=T)
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".adjusted.txt"))
    }
}
}

l <- ldply(lapply(c("sex","age","sex_age_int",psychvarstorun), function(var){
    ml<-ldply(lapply(unique(metadata_ls[[1]]$treats),function(i){
                if(file.exists(paste0(outFolder,"adjusted/",project,".",cluster,".deseqres_",var,"-",i,".adjusted.txt"))){
        cat("running ",var," ",i,"\n")
    df <- fread(paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".txt"))
    df_adj <- fread(paste0(outFolder,"adjusted/",project,".stats_all_cell_types-",var,"-",i,".adjusted.txt"))
    m <- merge(df[,c("variable","cluster","DEGs_FDR_10")],df_adj[,c("variable","cluster","DEGs_FDR_10")],by=c("variable","cluster"))
    colnames(m) <- c("variable","cluster","DEGs_FDR_10","COMBAT_DEGs_FDR_10")
    m <- transform(m, treat=i,diff=DEGs_FDR_10-COMBAT_DEGs_FDR_10)
    m_s <- subset(m, COMBAT_DEGs_FDR_10<DEGs_FDR_10)
    m_s <- subset(m_s, DEGs_FDR_10>10 & diff>10)
    return(m_s)
        }
    }),data.frame)
    return(ml)
}), data.frame)

l <- ldply(lapply(c("sex","age","sex_age_int",psychvarstorun), function(var){
    ml<-ldply(lapply(unique(metadata_ls[[1]]$treats),function(i){
                if(file.exists(paste0(outFolder,"adjusted/",project,".",cluster,".deseqres_",var,"-",i,".adjusted.txt"))){
        cat("running ",var," ",i,"\n")
    df <- fread(paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".txt"))
    df_adj <- fread(paste0(outFolder,"adjusted/",project,".stats_all_cell_types-",var,"-",i,".adjusted.txt"))
    m <- merge(df[,c("variable","cluster","DEGs_FDR_10")],df_adj[,c("variable","cluster","DEGs_FDR_10")],by=c("variable","cluster"))
    colnames(m) <- c("variable","cluster","DEGs_FDR_10","COMBAT_DEGs_FDR_10")
    m <- transform(m, treat=i,diff=DEGs_FDR_10-COMBAT_DEGs_FDR_10)
    m_s <- subset(m, COMBAT_DEGs_FDR_10>DEGs_FDR_10)
    m_s <- subset(m_s, COMBAT_DEGs_FDR_10>10 & diff<(-10))
    return(m_s)
        }
    }),data.frame)
    return(ml)
}), data.frame)

#combat vs no combat DEGs, coloring by variable
l <- ldply(lapply(c("sex","age","sex_age_int",psychvarstorun), function(var){
    ml<-ldply(lapply(unique(metadata_ls[[1]]$treats),function(i){
                if(file.exists(paste0(outFolder,"adjusted/",project,".",cluster,".deseqres_",var,"-",i,".adjusted.txt"))){
        cat("running ",var," ",i,"\n")
    df <- fread(paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".txt"))
    df_adj <- fread(paste0(outFolder,"adjusted/",project,".stats_all_cell_types-",var,"-",i,".adjusted.txt"))
    m <- merge(df[,c("variable","cluster","DEGs_FDR_10")],df_adj[,c("variable","cluster","DEGs_FDR_10")],by=c("variable","cluster"))
    colnames(m) <- c("variable","cluster","DEGs_FDR_10","COMBAT_DEGs_FDR_10")
    m <- transform(m, treat=i)
    return(m)
        }
    }),data.frame)
    return(ml)
}), data.frame)

p <- ggplot(l, aes(x=DEGs_FDR_10, y=COMBAT_DEGs_FDR_10)) +
  theme_bw()+
  geom_point(aes(color=variable),size=3)+ #
  geom_vline(xintercept = 0)+
  #facet_wrap(.~treat,scales="free")+
  geom_smooth(method = "lm", fill = NA)+
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))

png(width = 12, height = 12, file=paste0(figuredir,project,".combatvsnocombat.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

#now just focusing on the more DEGs variables
#combat vs no combat DEGs, coloring by variable
l <- ldply(lapply(c(c("pr_comp","isel","SES","Chol_HDL","PSS_all_mean","SNI_NoP","DED_all_mean","BPd_avg","Trig")), function(var){
    ml<-ldply(lapply("RNA-CTRL",function(i){
                if(file.exists(paste0(outFolder,"adjusted/",project,".",cluster,".deseqres_",var,"-",i,".adjusted.txt"))){
        cat("running ",var," ",i,"\n")
    df <- fread(paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".txt"))
    df_adj <- fread(paste0(outFolder,"adjusted/",project,".stats_all_cell_types-",var,"-",i,".adjusted.txt"))
    m <- merge(df[,c("variable","cluster","DEGs_FDR_10")],df_adj[,c("variable","cluster","DEGs_FDR_10")],by=c("variable","cluster"))
    colnames(m) <- c("variable","cluster","DEGs_FDR_10","COMBAT_DEGs_FDR_10")
    m <- transform(m, treat=i)
    return(m)
        }
    }),data.frame)
    return(ml)
}), data.frame)

p <- ggplot(l, aes(x=DEGs_FDR_10, y=COMBAT_DEGs_FDR_10)) +
  theme_bw()+
  geom_point(aes(color=variable),size=3)+ #
  geom_vline(xintercept = 0)+
  facet_wrap(.~cluster,scales="free")+
  geom_smooth(method = "lm", fill = NA)+
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))

png(width = 12, height = 12, file=paste0(figuredir,project,".combatvsnocombat_sub.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

varrun <- psychvarstorun[!psychvarstorun %in% c("CVDRISK","SES","pr_comp","isel","SNI_NoP","Trig","BPd_avg","PSS_all_mean")]
for (cluster in names(counts_ls)){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))

        mclapply(varrun,function(var){
#var="SES"
        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","treats","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age",var)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".",var,"_sex_age_and_treats_adjusted.RData")
    save(adjusted_counts, file=opfn)
        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".",var,"_sex_age_and_treats_adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
}
})
}

cluster="C7"
for (var in c("sex","age","sex_age_int",psychvarstorun)){
    #for (var in c("factor_HS_CRP","CVDRISK")){
    for (i in unique(metadata_ls[[1]]$treats)){
                if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"),header=T)
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
    }
}

noncorvars <- psychvarstorun[!psychvarstorun %in% c("cytocomp","DSES_07","DSES_09","LogCRP","NAIscr2019","StressCount","NII_fam","Chol_LDL")]
mclapply(names(counts_ls)[c(7:8)],function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))

        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","treats","sex_alph","age",psychvarstorun)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age",noncorvars)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".","allvars_sex_age_and_treats_adjusted.RData")
    save(adjusted_counts, file=opfn)
})
#take 2
shortvars <- psychvarstorun[!psychvarstorun %in% c(removed,removed2,removed3,ilremove,cholremove,other)]
mclapply(names(counts_ls),function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))

        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","treats","sex_alph","age",psychvarstorun)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age",shortvars)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".","shortvars_sex_age_and_treats_adjusted.RData")
    save(adjusted_counts, file=opfn)
})

mclapply(names(counts_ls),function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".","shortvars_sex_age_and_treats_adjusted.RData")
    load(opfn)

        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        mclapply(psychvarstorun[!psychvarstorun %in% c("SES","pr_comp","isel","DED_all_mean","Chol_HDL","SNI_NoP","Trig","BPd_avg","PSS_all_mean")],function(var){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".","shortvars_sex_age_and_treats_adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))
})
    }
})


cluster="C7"
for (var in c("sex","age","sex_age_int",psychvarstorun)){
    #for (var in c("SES"))
    for (i in unique(metadata_ls[[1]]$treats)){
                if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"),header=T)
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".",var,"_sex_age_and_treats_adjusted.txt"))
    }
}