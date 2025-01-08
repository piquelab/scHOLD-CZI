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
library(Seurat)
library(scuttle)
library(tidyverse)
library(data.table)
library(plyr)
library(parallel)
library("AnnotationHub")
library(ggseurat)
library(cowplot)
library(sva)
library(ggpubr)
future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_SNI_NoP_LogDED-Trig_chronic-sum-categ_nii-mean_11_25_2024.txt","ALL","fastdemux") #for testing
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
eigenvec2 <- merge(eigenvec2_o[,-c("sex","sex_alph","age")],cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip","SNI_NoP","age","sex","sex_alph") #SNI_NoP is the only variable that should be excluded based on the observed issues in score distributions that don’t make sense (negative values and extreme outliers)
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
focus_vars <-c("factor_HS_CRP","HS_CRP","PSQI_total","HVS_mean","EDS_mean","SES","LogCRP","cytocomp","SLS","PSS_all_mean","DED_all_mean","Chol_HDL","Chol_LDL","BPs_avg","BPd_avg")

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

#opfn_i <- file.info(dir(paste0(base,"5b_IdenCelltype_",method,"/"), full.names=T, pattern=paste0(project,".seuratObj-.harmony-sctype-")))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)
sc <- read_rds("/rs/rs_grp_schold/CZI/RNA/analysis/5b_IdenCelltype_fastdemux/ALL.seuratObj-.harmony-sctype-2024-05-15.rds")
#from https://www.biostars.org/p/9482789/

sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
colOrd <- data.frame(NEW_BARCODE=names(sc$orig.ident))
coldata_ex <- merge(as.data.frame(sc@meta.data),eigenvec2,by="Sample_ID") #first 2 columns are dbgapid ie. sampleid and batch
RcolData <- merge(colOrd, coldata_ex, by="NEW_BARCODE",all.x=T)
RcolData <- RcolData[match(rownames(sc@meta.data), RcolData$NEW_BARCODE), ]
sc@meta.data <-cbind(sc@meta.data,RcolData[,which(!colnames(RcolData) %in% colnames(sc@meta.data))])
#sc@meta.data <- sc@meta.data[,-c(27:29)]
#names(sc@meta.data)[22:24] <- c("age","sex","sex_alph")
#if(!is.na(args[3])){
#sc@meta.data$keep <- ifelse(sc@meta.data$Sample_ID %in% removeind, FALSE, TRUE)
#sc <- subset(sc, subset = keep == TRUE)
#sc$keep <- NULL
#}

sc[["sex_male"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.male$symbol])
sc[["sex_female"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.female$symbol])

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

cellcount <- as.data.frame(table(sc@meta.data$letter_clusters, sc@meta.data$orig.ident))
lowcell <- subset(cellcount, Freq<3000)
if(dim(lowcell)[1]==0){cat( "no low cell counts")}else{sc <-subset(x = sc, subset = letter_clusters %in% unique(lowcell$Var1), invert = TRUE)}

sce <- as.SingleCellExperiment(sc)
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".SingleCellExperiment.RDS")
write_rds(sce, opfn)

#seurat_clusters, treats, BATCH, Library
#/ aggregate by cluster,library info and covariate of interest:
sum_by <- c("letter_clusters", "BATCH","Sample_ID","treats")
summed <- aggregateAcrossCells(sce, id=colData(sce)[,sum_by])

#/ add rownames using the information from the colData:
colnames(summed) <- apply(colData(summed)[,sum_by], 1, function(x) paste(x, collapse="_"))
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

#opfn <- paste0(outFolder,project,".DESeq_countlists.RData")
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".DESeq_countlists.RData")
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
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
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
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
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
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
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
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
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
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
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
shortvars <- c("FCDEM_12", "smoke", "CVDRISK", "PSQI_total", "SES", "cytocomp", "pr_comp", "StressSev", "SNI_NoP", "DED_all_mean", "Chol", "BPd_avg")
mclapply(names(counts_ls),function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)

        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","treats","sex_alph","age",shortvars)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age",shortvars)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".","shortvars_sex_age_and_treats_adjusted.RData")
    save(adjusted_counts, file=opfn)
})
#[!psychvarstorun %in% c(removed,removed2,removed3,ilremove,cholremove,other)]

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

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".shortvars_sex_age_and_treats_adjusted.RData")
    load(opfn)

        for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-LPS"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        mclapply(psychvarstorun,function(var){
            #var="SES"
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
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
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
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

#cluster="C6"
for (var in c("sex","age","sex_age_int",psychvarstorun)){
    #for (var in psychvarstorun){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".shortvars_sex_age_and_treats_adjusted.txt"))
    }
}

#filter: ALL.stats_all_cell_types*.shortvars*

vars=fread(file=paste0(base,"COMBAT_modelvars.txt"),header=F)$V1
mclapply(names(counts_ls),function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)

        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","treats","sex_alph","age",vars)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age",vars)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".vars_sex_age_and_treats_adjusted.RData")
    save(adjusted_counts, file=opfn)
})
#[!psychvarstorun %in% c(removed,removed2,removed3,ilremove,cholremove,other)]

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

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".vars_sex_age_and_treats_adjusted.RData")
    load(opfn)

    for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-LPS"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        mclapply(psychvarstorun,function(var){
            #var="SES"
            if(!file.exists(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)

        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        cat("running deseq ",var,cluster,i," \n")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".","vars_sex_age_and_treats_adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))
            }
        })
    }
})

#cluster="C6"
for (var in c("sex","age","sex_age_int",psychvarstorun)){
    #for (var in psychvarstorun){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".vars_sex_age_and_treats_adjusted.txt"))
    }
}
#filter: ALL.stats_all_cell_types*.vars*

#only using vars that have no missing data
nomissvars=fread(file=paste0(base,"COMBAT_modelvars_nomiss.txt"),header=F)$V1
mclapply(names(counts_ls),function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)

        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","treats","sex_alph","age",nomissvars)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age",nomissvars)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".nomissvars_sex_age_and_treats_adjusted.RData")
    save(adjusted_counts, file=opfn)
})
#[!psychvarstorun %in% c(removed,removed2,removed3,ilremove,cholremove,other)]

firstrunvars=c("SES","pr_comp","isel","SNI_NoP","PSS_all_mean","BPd_avg")
outFolder=paste0(outFolder,"combat_nomissvars/")
mclapply(names(counts_ls),function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[5]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".nomissvars_sex_age_and_treats_adjusted.RData")
    load(opfn)

    for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        mclapply(psychvarstorun,function(var){
            #var="BPd_avg"
            if(!isTRUE(file.size(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt")) > 0)){
                cat("running deseq ",var,cluster,i," \n")
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".nomissvars_sex_age_and_treats_adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))
            }
        })
    }
})

#cluster="C6"
for (var in c("sex","age","sex_age_int",psychvarstorun)){
    #for (var in firstrunvars){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".nomissvars_sex_age_and_treats_adjusted.txt"))
    }
}
#filter: ALL.stats_all_cell_types*.nomissvars*



#only using vars that have no missing data, also accounting for PC1 and PC2
nomissvars=fread(file=paste0(base,"COMBAT_modelvars_nomiss.txt"),header=F)$V1
mclapply(names(counts_ls),function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)

        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","PC1","PC2","treats","sex_alph","age",nomissvars)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age","PC1","PC2",nomissvars)])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".nomissvars_PCs_sex_age_and_treats_adjusted.RData")
    save(adjusted_counts, file=opfn)
})

outFolder=paste0(outFolder,"combat_nomissvars_pcs/")
firstrunvars=c("SES","pr_comp","isel","SNI_NoP","PSS_all_mean","BPd_avg")
lapply(names(counts_ls),function(cluster){
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
    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".nomissvars_PCs_sex_age_and_treats_adjusted.RData")
    load(opfn)

    for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        mclapply(psychvarstorun,function(var){
            #var="SNI_NoP"
            if(!isTRUE(file.size(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt")) > 0)){
                cat("running deseq ",var,cluster,i," \n")
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".nomissvars_PCs_sex_age_and_treats_adjusted.RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))
            }
        })
    }
})

#cluster="C6"
for (var in c("sex","age","sex_age_int",psychvarstorun)){
    #for (var in firstrunvars){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".nomissvars_PCs_sex_age_and_treats_adjusted.txt"))
    }
}
#filter: ALL.stats_all_cell_types*.nomissvars_PCs*

#only using SES var, also accounting for PC1 and PC2
mclapply(names(counts_ls),function(cluster){
 #[c(7:9)] #to run clusters if job terminates part way through
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)

        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","PC1","PC2","treats","sex_alph","age","SES")]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age","PC1","PC2","SES")])

    opfn <- paste0(outFolder,project,".ComBat_seq.",cluster,".SES_PCs_sex_age_and_treats_adjusted.RData")
    save(adjusted_counts, file=opfn)
})

run="SES_PCs_sex_age_and_treats_generem"
combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
firstrunvars=c("SES","pr_comp","isel","SNI_NoP","PSS_all_mean","BPd_avg")
updated_vars <- c("SNI_NumPeople_r", "chronic_sum_categ", "nii_mean", "LogDED", "LogTrig")
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid/"))) dir.create(paste0(outFolder,"sampleid/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid_all/"))) dir.create(paste0(outFolder,"sampleid_all/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

lapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    opfn <- paste0(baseoutFolder,project,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)

    d <- data.frame(number_ind=length(unique(cluster_metadata[,c("Sample_ID")])))
    fwrite(d, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid_all/",project,".",cluster,".sampleid_all.",run,".txt"))

    for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        mclapply(c(psychvarstorun,"age","factor_HS_CRP"),function(var){
            #c(psychvarstorun,"factor_HS_CRP",updated_vars)
            #var="SNI_NoP" ,"factor_HS_CRP"
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
                cat("running deseq ",var,cluster,i," \n")
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)
        if(var=="CVDRISK"){
        design <-  paste0("~ PC1 + PC2 + ",var)
        }
        if(var=="factor_HS_CRP"){
        cluster_metadata_t <- subset(cluster_metadata_t, HS_CRP<10) #advised to remove as likely an infection
        cluster_metadata_t <- transform(cluster_metadata_t, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
        cluster_metadata_t <- within(cluster_metadata_t, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
        }
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        if(var=="age"){
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph",var)]
        }
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".",run,".RDS")
        saveRDS(dds, file=opfn)

        fwrite(data.frame(Sample_ID=cluster_metadata_var[,c("Sample_ID")]), sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid/",project,".",cluster,".sampleid_",var,"-",i,".",run,".txt"))

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$tested_genes <- paste(nrow(sub.table))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
            }
        })
    }
})


for (var in c("sex","age","sex_age_int",psychvarstorun,"factor_HS_CRP")){
    #for (var in c("factor_HS_CRP","CVDRISK")){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".deseqres_",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".deseqres_",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".sigDEGs_",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".sigDEGs_",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"))){
        sub.table <- fread(paste0(outFolder,project,".",cluster,".stats_all_cell_types-",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".stats_all_cell_types-",var,"-",i,".SES_PCs_sex_age_and_treats_adjusted.txt"))
    }
}
#filter: ALL.stats_all_cell_types*.SES_PCs_sex*

run="SES_PCs_SES_sex_age_and_treats_adjusted"
combatrun="SES_PCs_sex_age_and_treats_adjusted"
firstrunvars=c("SES","pr_comp","isel","SNI_NoP","PSS_all_mean","BPd_avg")
updated_vars <- c("SNI_NumPeople_r", "chronic_sum_categ", "nii_mean", "LogDED", "LogTrig")
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
combatfolder=paste0(baseoutFolder,"combat_sessexagetreatpcs/")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid/"))) dir.create(paste0(outFolder,"sampleid/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid_all/"))) dir.create(paste0(outFolder,"sampleid_all/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

lapply(names(counts_ls),function(cluster){
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
    opfn <- paste0(combatfolder,project,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)

    d <- data.frame(number_ind=length(unique(cluster_metadata[,c("Sample_ID")])))
    fwrite(d, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid_all/",project,".",cluster,".sampleid_all.",run,".txt"))


    for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        lapply(updated_vars,function(var){
            #var="SNI_NoP" ,"factor_HS_CRP"
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
                cat("running deseq ",var,cluster,i," \n")
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + SES + ",var)
        if(var=="SES"){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)
        }
        if(var=="CVDRISK"){
        design <-  paste0("~ PC1 + PC2 + SES + ",var)
        }
        if(var=="factor_HS_CRP"){
        cluster_metadata_t <- subset(cluster_metadata_t, HS_CRP<10) #advised to remove as likely an infection
        cluster_metadata_t <- transform(cluster_metadata_t, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
        cluster_metadata_t <- within(cluster_metadata_t, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
        }
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","SES",var)]
        if(var=="SES"){
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        }
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".DESeq_output-",i,"-",var,cluster,".",run,".RDS")
        saveRDS(dds, file=opfn)

        fwrite(data.frame(Sample_ID=cluster_metadata_var[,c("Sample_ID")]), sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid/",project,".",cluster,".sampleid_",var,"-",i,".",run,".txt"))

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"))
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$tested_genes <- paste(nrow(sub.table))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
            }
        })
    }
})

#combine output per cluster
sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sampleid_all/",project,".",cluster,".sampleid_all.",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sampleid_all/",project,".",cluster,".sampleid_all.",run,".txt"),header=T)
        }
        sub.table <- transform(sub.table,cluster=cluster)
        }),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid_all/",project,".sampleid_all_combtreat.",run,".txt"))

for (var in c("age",psychvarstorun,"factor_HS_CRP")){ #"sex","sex_age_int"
    #for (var in c(firstrunvars[c(5:10)],"age")){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sampleid/",project,".",cluster,".sampleid_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sampleid/",project,".",cluster,".sampleid_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid/",project,".sampleid_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",cluster,".deseqres_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".deseqres_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sigDEGs/",project,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sigDEGs/",project,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".sigDEGs_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    }
}
#filter: ALL.stats_all_cell_types*.SES_PCs_SES_*

################## UNDER CONSTRUCTION ##########################
#################################################################

#automated table making
#- top is cluster; sub top is num ind. | num genes | num DEG
#- side is variable
#one per each treatment
# assign variabeles
# assign variabeles
variables <- c("DSES_07","DSES_09","FCDEM_12","chronic_sum","smoke","CVDRISK",
                "HS_CRP","il6","PSQI_total","HVS_mean","EDS_mean","SES","LogCRP","Logifny","Logil10","Logil12",
                "Logil13","Logil1b","Logil2","Logil4","Logil6","Logil8","Logtnfa","cytocomp","NAIscr2019","NAIscr2017",
                "LQ2019","StressSev","StressCount","NII_mean","NII_fam","SLS","LivingAlone","SNI_HCG","SNI_NoP","PSS_all_mean",
                "DED_all_mean","Chol","Trig","Chol_HDL","Chol_Ratio","Chol_LDL","HbA1C","BPs_avg","BPd_avg","isel", "pr_comp",
                "SNI_NumPeople_r","LogDED","LogTrig","chronic_sum_categ","nii_mean",
                "age","factor_HS_CRP"
                )
 
variable_names <- c("Participant's highest level of education","Pre-tax household income","Are you currently taking any prescription medications",
                    "Total number of chronic conditions","Smoking status","Cardiovascular disease risk",
                    "High Sensitivity CRP (mg/L) raw values","IL-6 raw values","Total sleep score (PSQI)","Vigilance Scale - Mean","Everyday Discrimination Scale - Mean",
                    "SES composite [mean(ZDSES_07, ZDSES_09)]","LogCRP","Logifny","Logil10","Logil12",
                    "Logil13","Logil1b","Logil2","Logil4","Logil6",
                    "Logil8","Logtnfa","Cytokines composite [mean(ZLogtnfa, ZLogi6, ZLogifny)]","Neighborhood Adversity Score","Neighborhood Adversity Score",
                    "Neighborhood Segregation Score","Core: Total Severity of Stressors","Core: Total Count of Stressors","Negative interactions mean","Negative interactions with family","Avg on loneliness scale",
                    "Dichotomous variable measuring living alone (1 = living alone)","High contact groups: number of categories of social groups p interacted with at least once every two weeks",
                    "This is the mean score representing the number of people with whom the respondent has regular contact","Mean of daily pss measures","Daily experience of discrimination cumulative average",
                    "Cholesterol in mg/dL","Triglycerides in mg/dL","HDL Cholesterol in mg/dL","CHOL/HDL ratio","LDL Cholesterol in mg/dL",
                    "Glycated HGB Affinity HPLC HbA1C","Average systolic blood pressure (wonky BP3_s measure corrected)","Average diastolic blood pressure","social support","psychological resources composite",
                    "number of pople in contact with","Log of Daily Discrimination", "LogTrig", "Chronic conditions transformed", "Negative interactions mean",
                    "age","High Sensitivity CRP (mg/L) factorized 1 to 5"
                    )




# Variables to loop over
variablesL <- c("DSES_09","chronic_sum","SES", "Logifny", "Logil4", "Logtnfa", "cytocomp",
                "PSS_all_mean", "Chol_HDL", "BPd_avg", "isel", "pr_comp", "SNI_NumPeople_r", "LogDED", "LogTrig", "chronic_sum_categ", "nii_mean") 



variable_namesL <- c("Income", "Chronic conditions","SES composite", "Logifny", "Logil4", "Logtnfa", "Cytokines composite", 
                    "Perceived Stress Scale", "HDL Cholesterol", "Diastolic blood pressure","social support","psychological resources composite", 
                    "number of pople in contact with", "Log of Daily Discrimination", "LogTrig", "Chronic conditions transformed", "Negative interactions mean")


variables_df <- data.frame(variable=variables, description=variable_names)

library(flextable)
library(ftExtra)
library(rlist)
#[!psychvarstorun %in% puberty]
for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
subvars <- ldply(lapply(variables_df$variable, function(var){
    cat("running",i,var,"\n")
    if(isTRUE(file.size(paste0(outFolder,"stats/",project,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
    stats <- fread(paste0(outFolder,"stats/",project,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".deseqres_",var,"-",i,".",run,".txt"),select=c(1,7))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,list.cbind(newls))
            return(newtable)
        }
    } 
}),data.frame)

dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",i,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
#subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=max)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(names(which(apply(subvars[,degcols],MARGIN=2,FUN=max)<50)),"[.]"),function(y)y[1]),NA),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",i,".",run,".50degs.png"))
}

