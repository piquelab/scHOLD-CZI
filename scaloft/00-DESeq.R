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
library(sva)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

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
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_n521_uniq-n265_04-03-2024.txt","ALL","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt",0.2) 

base <- args[1]

cov_file=fread(args[2]) #this is the psych cov file

#if(!is.na(args[3])){
#removeind <- fread(paste0(base,args[3]),header=F)$V1
#project <- sapply(strsplit(args[3],"_"),function(y)y[1])
#}
project=args[3]
method=args[4]
sample_batch <- args[5]
resset <- args[6]
dimset=13

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
#for resolution 0.2 and including V2 chem
outFolder=paste0(outFolder,"lessfilt/")

if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
#eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
#eigenvec2 <- merge(eigenvec2_o,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)

exp=fread(paste0(base,"filtered.incV2.",gsub(".*/","",sample_batch)))
eigenvec2 <- merge(exp,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)

#notrun_var <- c("DSES_01","DSES_03","PWaist","PHip")
#colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
#psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
#psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
psychvarstorun <- colnames(eigenvec2[,22:57])

outdir=paste0(base,"5b_IdenCelltype_",method,"/")
opfn <- paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
sc <- read_rds(opfn)

#from https://www.biostars.org/p/9482789/

sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
colOrd <- data.frame(NEW_BARCODE=names(sc$orig.ident))
coldata_ex <- merge(as.data.frame(sc@meta.data),eigenvec2,by.x=c("Sample_ID","BATCH"),by.y=c("Sample_ID","Batch2"),all.x=T) #first 2 columns are dbgapid ie. sampleid and batch
RcolData <- merge(colOrd, coldata_ex, by="NEW_BARCODE",all.x=T)
RcolData <- RcolData[match(rownames(sc@meta.data), RcolData$NEW_BARCODE), ]
sc@meta.data <-cbind(sc@meta.data,RcolData[,which(!colnames(RcolData) %in% colnames(sc@meta.data))])
sc@meta.data$orig.ident <- "scaloft_comb"

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
p1 <- p + geom_boxplot(aes(x=Sex, y=sex_male,colour=Sex))
p2 <- p + geom_boxplot(aes(x=Sex, y=sex_female,colour=Sex))
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
    Sex=unique(Sex))

p <- ggplot(data = mean_sex) +
  cowplot::theme_cowplot()
p1 <- p + geom_boxplot(aes(x=Sex, y=sex_male,colour=Sex))
p2 <- p + geom_boxplot(aes(x=Sex, y=sex_female,colour=Sex))
fig <-plot_grid(p1, p2, nrow=2, ncol=1, align="h")
png(paste0(figuredir,project,".box_ncount_sex_averageind.png"), width=3000, height=1000, res=120)
print(fig)
dev.off()

sce <- as.SingleCellExperiment(sc)
#seurat_clusters, treats, BATCH, Library
opfn <- paste0(outFolder,project,".",resset,".",dimset,".singlecellexp.RData")
save(sce, file=opfn)
#load(opfn)

library(biomaRt)  
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",1:22), mart = mart)

#/ aggregate by cluster,library info and covariate of interest:
sum_by <- c("letter_clusters", "BATCH","Sample_ID","treats")
summed <- aggregateAcrossCells(sce, id=colData(sce)[,sum_by])

#/ add rownames using the information from the colData:
colnames(summed) <- apply(colData(summed)[,sum_by], 1, function(x) paste(x, collapse="_"))
#head(assay(summed, "counts"))
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
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(data) / 4) #filter to keep genes expressed in 25% of samples
    data <- unlist(data[keep, ])
    data <- data[rownames(data) %in% genes$hgnc_symbol, ]
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
    test <- metadata %>% dplyr::select(-c(sex_male,sex_female,orig.ident,NEW_BARCODE,percent.mt,nCount_RNA,nFeature_RNA,NUM.READS,NUM.SNPS,RNA_snn_res.0.2))
    #need to remove all RNA_snn_res columns so one per each resolution
    #test <- test %>% dplyr::select(-c(RNA_snn_res.0.3,RNA_snn_res.0.4))
    test <- unique(transform(test,cluster_sample_id=paste(letter_clusters,BATCH,Sample_ID,treats,sep="_")))
    #test[test$cluster_sample_id=="C0_SCAIP10_AL-021_CTRL",c(1:20)]
    test <- test[!duplicated(test$cluster_sample_id),]
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
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df_n) <- df_n$cluster_sample_id
    #df_n <- df_n[complete.cases(df_n), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    return(df_n)
})

all(names(counts_ls) == names(metadata_ls))

opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists.RData")
save(counts_ls,metadata_ls, file=opfn)
#load(opfn)

mclapply(names(counts_ls), function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)

    #I commented these out here, as it caused COMBAT to return the confounded error
    #cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    #cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))

    x <- subset(as.data.frame(table(cluster_metadata$BATCH)),Freq>0)
    cluster_metadata <- subset(cluster_metadata,BATCH %in% unique(x$Var1))
    cluster_metadata_var <- unique(cluster_metadata[,c("Sample_ID","BATCH","treats","Sex","cage1")])
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
    
    adjusted <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","Sex","cage1")])
    opfn <- paste0(outFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".sex_age_and_treats_adjusted.RData")
    save(adjusted, file=opfn)
})

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
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))

    opfn <- paste0(outFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".sex_age_and_treats_adjusted.RData")
    load(opfn)

    cluster_metadata_var <- cluster_metadata[,c("Sample_ID","treats")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts) == rownames(cluster_metadata_var))

    cat("running treatment deseq \n")
    # Create DESeq2 object        
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata_var, 
                                  design = ~ Sample_ID + treats)
    #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-treats",cluster,".sex_age_and_treats_adjusted.RData")
    save(dds, file=opfn)

    var="LPSvsCTRL"
    cluster_metadata_t <- unique(subset(cluster_metadata[,c("Sample_ID","treats")], treats %in% c("LPS","CTRL")))
    cluster_metadata_var <- cluster_metadata_t[complete.cases(cluster_metadata_t),] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    res <- results(dds, contrast=c("treats","LPS","CTRL"))
    sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
    names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
    sub.table <- sub.table[!is.na(sub.table$padj), ]
    cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    sub.table$var=var
    sub.table$cluster=cluster
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,".sex_age_and_treats_adjusted.txt"))
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,".sex_age_and_treats_adjusted.txt"))
    table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(adjusted))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_logfc2 <- paste(nrow(subset(sub.table,padj<fdr & abs(logFC>0.5))))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,".sex_age_and_treats_adjusted.txt"))
}

var="LPSvsCTRL"
clustersrun <- ldply(lapply(names(counts_ls), function(cluster){
    if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,".txt"))){
        return(cluster)
    }
}),data.frame)
sub.table <- ldply(lapply(clustersrun$X..i.., function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",resset,".",dimset,".",cluster,".deseqres_",var,".txt"),header=T)
    }),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".deseqres_",var,".txt"))
sub.table <- ldply(lapply(clustersrun$X..i.., function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,".txt"),header=T)
}),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".sigDEGs_",var,".txt"))
sub.table <- ldply(lapply(clustersrun$X..i.., function(cluster){
    sub.table <- fread(paste0(outFolder,project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,".txt"),header=T)
}),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".stats_all_cell_types-",var,".txt"))

filenames <- list.files(outFolder) #file list from directory
filenames1 <- filenames[grep("DESeq_output-", filenames)] #pick specific files from list
x <- na.omit(sapply(strsplit(filenames1,"(?<=.)(?=C[0-9])",perl=TRUE),function(y)y[2]))
clusters <- unique(sapply(strsplit(x,"[.]"),function(z)z[1]))

#for (cluster in clusters){
cluster=clusters[1]
    cat("running ", cluster, "\n")
    var="LPSvsCTRL"
    load(paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-treats",cluster,".RData"))

    vsd <- vst(dds, blind=FALSE)
    # Plot PCA
    fig0 <- DESeq2::plotPCA(vsd, intgroup = c("treats"))+ggtitle("treats")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    #fig1 <- DESeq2::plotPCA(vsd, intgroup = c("PC1"))+ggtitle("PC1")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    #fig2 <- DESeq2::plotPCA(vsd, intgroup = c("PC2"))+ggtitle("PC2")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig3 <- DESeq2::plotPCA(vsd, intgroup = c("Sex"))+ggtitle("sex")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig4 <- DESeq2::plotPCA(vsd, intgroup = c("cage1"))+ggtitle("age")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig <-plot_grid(fig0,fig3,fig4, nrow=3, ncol=1, align="h")
    png(paste0(figuredir,project,".",resset,".",dimset,".deseqPCA-treats","-",cluster,".png"), width=1000, height=1000, res=120)
    print(fig)
    dev.off()
#}


puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))
    opfn <- paste0(outFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".sex_age_and_treats_adjusted.RData")
    load(opfn)
            #run covariates separately for each treatment condition
        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, Sex=as.factor(Sex))
        if(all(c("0","1") %in% levels(cluster_metadata_t$Sex))){
            cluster_metadata_t <- within(cluster_metadata_t, Sex <- relevel(Sex, ref = "0"))
        }
        mclapply(psychvarstorun[c(11:length(psychvarstorun))][!psychvarstorun[c(11:length(psychvarstorun))] %in% puberty],function(var){
        #var="SES"
            if(!isTRUE(file.size(paste0(outFolder,project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".sex_age_and_treats_adjusted.txt")) > 0)){

        design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + cage1 + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","genPC1","genPC2","genPC3","Sex","cage1",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".sex_age_and_treats_adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".sex_age_and_treats_adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".sex_age_and_treats_adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".sex_age_and_treats_adjusted.txt"))
            }
        })
    }
}


clustersrun <- ldply(lapply(names(counts_ls), function(cluster){
    i <- tail(unique(metadata_ls[[1]]$treats),n=1)
    var=tail(c(psychvarstorun)[c(1:3)],n=1)
    if(file.exists(paste0(outFolder,project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".sex_age_and_treats_adjusted.txt"))){
        return(cluster)
    }
}),data.frame)

for (var in c(psychvarstorun)[c(11:length(psychvarstorun))][!psychvarstorun[c(11:length(psychvarstorun))] %in% puberty]){
    for (i in unique(metadata_ls[[1]]$treats)){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".sex_age_and_treats_adjusted.txt"),header=T)
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".deseqres_",var,"-",i,".sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".sex_age_and_treats_adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".sigDEGs_",var,"-",i,".sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".sex_age_and_treats_adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".sex_age_and_treats_adjusted.txt"))
}
}


#puberty variables run separately as they run per each sex
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")[c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds") %in% colnames(cov_file)]
for (cluster in names(counts_ls)){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats), Sex=as.factor(ifelse(Sex==1,"female","male")))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))
    cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "male"))
    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".adjusted.RData")
    load(opfn)

            #run covariates separately for each treatment condition
        for(i in unique(cluster_metadata$treats)){
            for(s in c("male","female")[2]){
            #lapply(list.df, subset, B!=2)
            #i <- "CTRL"

        lapply(puberty,function(var){
        #var="SES"
        #var <- colnames(cluster_metadata_t[con])
        cluster_metadata_t <- subset(cluster_metadata, treats==i & Sex==s)
        design <-  paste0("~ genPC1 + genPC2 + genPC3 + cage1 + ",var)
        cat("running ",var,cluster,i,s," \n")

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","genPC1","genPC2","genPC3","cage1",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        if(dim(cluster_metadata_var)[1]>0) {

        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
        cat("running deseq ",var,cluster,i," sex= ",s," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,"_",s,".adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        sub.table$Sex =s
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,"_",s,".adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,"_",s,".adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,"_",s,".adjusted.txt"))
        }
        })
    }
    }
}


clustersrun <- ldply(lapply(names(counts_ls), function(cluster){
    i <- tail(unique(metadata_ls[[1]]$treats),n=1)
    var=tail(c(puberty),n=1)
    if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,"_",s,".adjusted.txt"))){
        return(cluster)
    }
}),data.frame)

for (var in puberty){
    for (i in unique(metadata_ls[[1]]$treats)){
        for (s in c("male","female")){
    if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,"_",s,".adjusted.txt"))){
        cat("running ",var," ",i," ",s,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,"_",s,".adjusted.txt"),header=T)
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,"_",s,".adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,"_",s,".adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,"_",s,".adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,"_",s,".adjusted.txt"),header=T)
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,"_",s,".adjusted.txt"))
}
}
}
}
