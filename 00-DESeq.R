library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(plyr)
library(parallel)
library("AnnotationHub")
library(ggseurat)

ah <- AnnotationHub()
if(length(ah["AH98047"]) == 0) {
  edb <- ah[["AH75011"]]
} else {
  edb <- ah[["AH98047"]]
}
geneIDs <- genes(edb) %>%
  as.data.frame() %>% 
  setDT(keep.rownames = "ensembl_gene_id") %>%
  .[, c("ensembl_gene_id","entrezid","symbol","seqnames","start","end","strand","gene_biotype", "description")]
names(geneIDs)[c(1,2,4,8)] <- c("ensgene","entrez","chr","biotype")
geneIDs.sex <- subset(geneIDs, chr=="Y" | chr=="X")
geneIDs.male <- subset(geneIDs, chr=="Y" )
geneIDs.female <- subset(geneIDs, chr=="X")

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_converted_dbgapID_key_n183_AR_18newvar_updated_08-05.txt","ALL","fastdemux",0.3) #for testing
base <- args[1]

cov_file=fread(args[2]) #this is the psych cov file

#if(!is.na(args[3])){
#removeind <- fread(paste0(base,args[3]),header=F)$V1
#project <- sapply(strsplit(args[3],"_"),function(y)y[1])
#}
project=args[3]
method=args[4]
resset <- args[5]

#uncoment this section if no longer loading in eigenvec pc file (contains all cov info already)
#cov_file=args[2]
#read in samples file (just list of samples to run, each sample on newline)
#samples=read.table(args[3],header=F)
#samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
#}
#exp <- read.table(cov_file, row.names=NULL,header=T)
#exp <- transform(exp, Sample_ID=dbgap.ID)
#exp$Batch <- gsub("HOLD0","HOLD",exp$Batch)
#if(!is.na(args[3])){
#  exp <- exp %>% dplyr::filter(Batch %in% samples$Batch)
#}

outFolder=paste0(base,method,"_pseudobulk_ctrl/")
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

#opfn_i <- file.info(dir(paste0(base,"5b_IdenCelltype_",method,"/"), full.names=T, pattern=paste0(project,".",resset,".seuratObj-.harmony-sctype-")))
opfn_i <- file.info(dir(paste0(base,"5b_IdenCelltype_",method,"/"), full.names=T, pattern=paste0(".seuratObj-.harmony-sctype-")))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

#from https://www.biostars.org/p/9482789/

sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
colOrd <- data.frame(NEW_BARCODE=names(sc$orig.ident))
coldata_ex <- merge(as.data.frame(sc@meta.data),eigenvec2,by="Sample_ID",all.x=T) #first 2 columns are dbgapid ie. sampleid and batch
RcolData <- merge(colOrd, coldata_ex, by="NEW_BARCODE",all.x=T)
RcolData <- RcolData[match(rownames(sc@meta.data), RcolData$NEW_BARCODE), ]
sc@meta.data <-cbind(sc@meta.data,RcolData[,which(!colnames(RcolData) %in% colnames(sc@meta.data))])

#if(!is.na(args[3])){
#sc@meta.data$keep <- ifelse(sc@meta.data$Sample_ID %in% removeind, FALSE, TRUE)
#sc <- subset(sc, subset = keep == TRUE)
#sc$keep <- NULL
#}

sc[["sex_male"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.male$symbol])
sc[["sex_female"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.female$symbol])

cellcount <- as.data.frame(table(sc@meta.data$letter_clusters, sc@meta.data$orig.ident))
lowcell <- subset(cellcount, Freq<3000)
if(dim(lowcell)[1]==0){cat( "no low cell counts")}else{sc <-subset(x = sc, subset = letter_clusters == lowcell$Var1, invert = TRUE)}

p <- ggplot(data = sc) +
  cowplot::theme_cowplot()
p1 <- p + geom_boxplot(aes(x=sex_alph, y=sex_male,colour=sex_alph))
p2 <- p + geom_boxplot(aes(x=sex_alph, y=sex_female,colour=sex_alph))
fig <-plot_grid(p1, p2, nrow=2, ncol=1, align="h")
png(paste0(figuredir,project,".box_ncount_sex.png"), width=3000, height=1000, res=120)
print(fig)
dev.off()

fig0 <- VlnPlot(sc, features = c("sex_male", "sex_female"), ncol = 2,pt.size = FALSE) 
png(paste0(figuredir,project,".violin_ncount_sex.png"), width=3000, height=1000, res=120)
print(fig0)
dev.off()

mean_sex <- ddply(sc@meta.data, "Sample_ID", plyr::summarize,
    sex_male=mean(sex_male, na.rm=T),
    sex_female=mean(sex_female, na.rm=T),
    sex_alph=unique(sex_alph))

p <- ggplot(data = mean_sex) +
  cowplot::theme_cowplot()
p1 <- p + geom_boxplot(aes(x=sex_alph, y=sex_male,colour=sex_alph))
p2 <- p + geom_boxplot(aes(x=sex_alph, y=sex_female,colour=sex_alph))
fig <-plot_grid(p1, p2, nrow=2, ncol=1, align="h")
png(paste0(figuredir,project,".box_ncount_sex_averageind.png"), width=3000, height=1000, res=120)
print(fig)
dev.off()

sce <- as.SingleCellExperiment(sc)
#seurat_clusters, treats, BATCH, Library

#/ aggregate by cluster,library info and covariate of interest:
sum_by <- c("letter_clusters", "BATCH","Sample_ID","treats")
summed <- aggregateAcrossCells(sce, id=colData(sce)[,sum_by])

#/ add rownames using the information from the colData:
colnames(summed) <- apply(colData(summed)[,sum_by], 1, function(x) paste(x, collapse="_"))
head(assay(summed, "counts"))
raw <- assay(summed, "counts")

counts_ls <- lapply(unique(summed$letter_clusters), function(i){
    cat("running ",i)
    #i <- unique(summed$letter_clusters)[1]
    cell_idx <- which(tstrsplit(colnames(summed), "_")[[1]] == i)
    keep <- which(colnames(raw) %in% colnames(summed[, cell_idx]))
    data <- raw[, keep]
    #filtered_data <- data > 0 
    filtered_data <- data
    filtered_data[filtered_data < 0] <- NA
    if(is.matrix(filtered_data)){
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(data) / 4) #filter to keep genes expressed in 50% of samples
    data <- unlist(data[keep, ])
    data <- data[!rownames(data) %in% geneIDs.sex$symbol, ]
    summed_filt <- summed[rownames(summed) %in% rownames(data), cell_idx]
    #table(rownames(summed_filt) %in% geneIDs.sex$symbol)
    return(summed_filt)
    } else {
    NULL
    }
  #counts_ls[[i]] <- summed[, cell_idx]
  #names(counts_ls) <-
    })
names(counts_ls) <-unique(summed$letter_clusters)
counts_ls[sapply(counts_ls, is.null)] <- NULL

# Number of cells per sample and cluster
t <- table(colData(sce)$Sample_ID,
           colData(sce)$letter_clusters)

metadata_ls <- lapply(counts_ls, function(i){
    #i <- counts_ls[[1]]
    df <- data.frame(cluster_sample_id = colnames(i)) ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$letter_clusters <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$BATCH  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    df$Sample_ID  <- tstrsplit(df$cluster_sample_id, "_")[[3]]
    df$treats  <- tstrsplit(df$cluster_sample_id, "_")[[4]]
    metadata <- as.data.frame(colData(sce))
    test <- metadata %>% dplyr::select(-c(sex_male,sex_female,orig.ident,NEW_BARCODE,percent.mt,nCount_RNA,nFeature_RNA,NUM.READS,NUM.SNPS,RNA_snn_res.0.3))

    df_n <- merge(df, unique(test), all.x=T)

    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df_n$letter_clusters))
    cell_counts <- t[, idx]
    ## Remove samples with zero cell contributing to the cluster
    cell_counts <- cell_counts[cell_counts > 0]
    ## Match order of cell_counts and sample_ids
    sample_order <- match(df_n$Sample_ID, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    
    ## Append cell_counts to data frame
    df_n$cell_count <- cell_counts
## Join data frame (capturing metadata specific to cluster) to generic metadata
    #df_n <- plyr::join(df_n, metadata, 
    #                 by = intersect(names(df_n), names(metadata)))
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df_n) <- df_n$cluster_sample_id
    #df_n <- df_n[complete.cases(df_n), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    return(df_n)
})

all(names(counts_ls) == names(metadata_ls))

opfn <- paste0(outFolder,project,".DESeq_countlists.RData")
save(counts_ls,metadata_ls, file=opfn)
#load(opfn)

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

    cluster_metadata_var <- cluster_metadata[,c("treats","Sample_ID")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

    cluster_counts <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts) == rownames(cluster_metadata_var))

    cat("running treatment deseq \n")
    # Create DESeq2 object        
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata_var, 
                                  design = ~ treats + Sample_ID)
    #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".DESeq_output-treat",cluster,".RData")
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
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,".txt"))
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,".txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_t))
    table$number_individuals <- paste(length(unique(cluster_metadata_t$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_logfc2 <- paste(nrow(subset(sub.table,padj<fdr & abs(logFC>0.5))))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".txt"))

    cluster_metadata_t <- subset(cluster_metadata, treats %in% c("RNA-LPS-DEX","RNA-LPS"))
    cluster_metadata_var <- cluster_metadata_t[,c("treats","Sample_ID")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

    cat("running treatment deseq \n")
    # Create DESeq2 object        
    dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                  colData = cluster_metadata_var, 
                                  design = ~ treats + Sample_ID)
    #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".DESeq_output-treat",cluster,".RData")
    save(dds, file=opfn)

    var="DEXvsLPS"
    res <- results(dds, contrast=c("treats","RNA-LPS-DEX","RNA-LPS"))
    sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
    names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
    sub.table <- sub.table[!is.na(sub.table$padj), ]
    cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    sub.table$var=var
    sub.table$cluster=cluster
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,".txt"))
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,".txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_t))
    table$number_individuals <- paste(length(unique(cluster_metadata_t$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts_t))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_logfc2 <- paste(nrow(subset(sub.table,padj<fdr & abs(logFC>0.5))))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".txt"))

    #var="DEXvsCTRL"
    #cluster_metadata_t <- subset(cluster_metadata, treats %in% c("RNA-LPS-DEX","RNA-CTRL"))
    #res <- results(dds, contrast=c("treats","RNA-LPS-DEX","RNA-CTRL"))
    #sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
    #names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
    #sub.table <- sub.table[!is.na(sub.table$padj), ]
    #cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    #sub.table$var=var
    #sub.table$cluster=cluster
    #fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,".txt"))
    #sigDEGs <- subset(sub.table,padj<fdr)
    #sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    #fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,".txt"))
    #table <- data.frame(symb=var, variable= var, cluster=cluster)
    #table$number_samples <- paste(nrow(cluster_metadata_t))
    #table$number_individuals <- paste(length(unique(cluster_metadata_t$Sample_ID)))
    #table$gene_number <- paste(nrow(cluster_counts))
    #table$DEGs_FDR <- paste(nrow(sigDEGs))
    #table$DEGs_FDR_logfc2 <- paste(nrow(subset(sub.table,padj<fdr & abs(logFC>0.5))))
    #fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".txt"))

        #run covariates separately for each treatment condition
        mclapply(unique(cluster_metadata_sce$treats),function(i){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","BATCH","sex_alph","age")]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        var="sex"

        cat("running deseq ",var,cluster," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = ~ PC1 + PC2 + BATCH + sex_alph + age)
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds, contrast=c("sex_alph","Female","Male"))
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_logfc2 <- paste(nrow(subset(sub.table,padj<fdr & abs(logFC>0.5))))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".txt"))

        ####
        var="age"
        cat("getting deseq results for ",var,cluster," \n")

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_logfc2 <- paste(nrow(subset(sub.table,padj<fdr & abs(logFC>0.5))))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".txt"))

        })
}

#   - Interaction term deseq (treat*sex, treat*age, sex*age)

for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[7]
        cat("running ", cluster, "\n")
            fdr=0.1
    var="treat_sex_int"
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    cluster_metadata <- transform(cluster_metadata, sex_alph=as.factor(sex_alph))
    cluster_metadata <- within(cluster_metadata, sex_alph <- relevel(sex_alph, ref = "Male"))

    cluster_metadata_var <- cluster_metadata[,c("Sample_ID","treats","PC1","PC2","BATCH","sex_alph","age")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

    cluster_counts <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts) == rownames(cluster_metadata_var))

    cat("running treatment deseq \n")
    # Create DESeq2 object        
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata_var, 
                                  design = ~ treats + BATCH + age + PC1 + PC2 + sex_alph + treats*sex_alph)
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
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,".txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_var))
    table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<fdr)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".txt"))

    var="treat_age_int"
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata_var, 
                                  design = ~ treats + BATCH + age + PC1 + PC2 + sex_alph + treats*age)
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
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,".txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_var))
    table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<fdr)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,".txt"))

    mclapply(unique(cluster_metadata_sce$treats),function(i){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","treats","PC1","PC2","BATCH","sex_alph","age")]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

    var="sex_age_int"
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata_var, 
                                  design = ~ BATCH + age + PC1 + PC2 + sex_alph + sex_alph*age)
    #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".RDS")
    save(dds, file=opfn)

    res <- results(dds)
    sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
    names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
    sub.table <- sub.table[!is.na(sub.table$padj), ]
    cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    sub.table$var=var
    sub.table$cluster=cluster
    sub.table$treats =i
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".txt"))
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_var))
    table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<fdr)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".txt"))
    })
}

#psychosocial variables
#int <- sapply(cov_file,function(x) is.integer(x))
#cov_int <- cov_file[,..int]

for (cluster in names(counts_ls)[-c(1:5)]){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    cluster_counts <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata))]

    all(colnames(cluster_counts) == rownames(cluster_metadata))
            #run covariates separately for each treatment condition
        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        var="SES"
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        design <-  paste0("~ BATCH + sex_alph + age + ",var)
        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".noPC.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".noPC.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".noPC.txt"))

        mclapply(psychvarstorun[c(27:28)],function(var){
        #var="SES"
        #var <- colnames(cluster_metadata_t[con])
        design <-  paste0("~ PC1 + PC2 + BATCH + sex_alph + age + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","BATCH","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(nrow(cluster_metadata_var))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".txt"))
        })
    }
}

#this is separate to have CRP as a factor (known way to hande the variable) and remove individuals >10
for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    cluster_counts <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata))]

    all(colnames(cluster_counts) == rownames(cluster_metadata))

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

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".RDS")
    save(dds, file=opfn)

    res <- results(dds)
    sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
    names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
    sub.table <- sub.table[!is.na(sub.table$padj), ]
    cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    sub.table$var=var
    sub.table$cluster=cluster
    sub.table$treats =i
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".txt"))
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_var))
    table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<fdr)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".txt"))
    }
}

#this is separate as CVD risk was calculated taking into account sex and age so need those out of hte model
for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    cluster_counts <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata))]

    all(colnames(cluster_counts) == rownames(cluster_metadata))

            #run covariates separately for each treatment condition
        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

var="CVDRISK"
        #var <- colnames(cluster_metadata_t[con])
        design <-  paste0("~ PC1 + PC2 + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".RDS")
    save(dds, file=opfn)

    res <- results(dds)
    sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
    names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
    sub.table <- sub.table[!is.na(sub.table$padj), ]
    cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    sub.table$var=var
    sub.table$cluster=cluster
    sub.table$treats =i
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".txt"))
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
    table$number_samples <- paste(nrow(cluster_metadata_var))
    table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<fdr)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".txt"))
    }
}

#running SES per each batch separately to see if different #individuals is causing difference
for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[1]
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    cluster_counts <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata))]

    all(colnames(cluster_counts) == rownames(cluster_metadata))
            #run covariates separately for each treatment condition
        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        lapply(unique(cluster_metadata_t$BATCH),function(b){
        var="SES"
        #b="HOLD1"
        #var <- colnames(cluster_metadata_t[con])
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)
        cluster_metadata_t <- subset(cluster_metadata_t, BATCH==b)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","BATCH","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i,b," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".",b,".RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$BATCH=b
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".",b,".deseqres_",var,"-",i,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".",b,".sigDEGs_",var,"-",i,".txt"))
        table <- data.frame(variable= var, cluster=cluster, BATCH=b)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".",b,".stats_all_cell_types-",var,"-",i,".txt"))
        })
    }
}

var="SES"
cluster="C7"
b="HOLD1"
    for (i in unique(metadata_ls[[1]]$treats)){
        if(file.exists(paste0(outFolder,project,".",cluster,".",b,".deseqres_",var,"-",i,".txt"))){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        cluster_metadata<- metadata_ls[[cluster]]
        ldply(lapply(unique(cluster_metadata$BATCH),function(b){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".",b,".deseqres_",var,"-",i,".txt"),header=T)
        }),data.frame)}),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".perBatch.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        cluster_metadata<- metadata_ls[[cluster]]
        ldply(lapply(unique(cluster_metadata$BATCH),function(b){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".",b,".sigDEGs_",var,"-",i,".txt"),header=T)
    }),data.frame)}),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".perBatch.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        cluster_metadata<- metadata_ls[[cluster]]
        ldply(lapply(unique(cluster_metadata$BATCH),function(b){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".",b,".stats_all_cell_types-",var,"-",i,".txt"),header=T)
    }),data.frame)}),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".perBatch.txt"))
    }
    }
l <- fread(paste0(outFolder,project,".stats_all_cell_types-",var,"-RNA-CTRL.perBatch.txt"))
p <- ggplot(l, aes(x=number_samples, y=DEGs_FDR_10)) +
  theme_bw()+
  geom_point(aes(color=BATCH),size=3)+ #
  facet_wrap(.~cluster,scales="free")+
  geom_smooth(method = "lm", fill = NA)+
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))

png(width = 12, height = 12, file=paste0(figuredir,project,".perBatch.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off() 


#combine clusters into 1 output
for (var in c("LPSvsCTRL","DEXvsLPS","DEXvsCTRL","treat_sex_int","treat_age_int")){
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

for (var in c("sex","age","sex_age_int",psychvarstorun)){
        #for (var in c("factor_HS_CRP","CVDRISK")){

    for (i in unique(metadata_ls[[1]]$treats)){
        if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".txt"))){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".txt"),header=T)
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".txt"))
    }
    }
}
 

var="SES"
    for (i in unique(metadata_ls[[1]]$treats)){
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".noPC.txt"),header=T)
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".noPC.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".noPC.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".noPC.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".noPC.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".noPC.txt"))
    }


