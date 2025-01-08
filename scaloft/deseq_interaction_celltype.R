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
library(ggseurat)
library(sva)
library(data.table)
library(flextable)
library(ftExtra)
library(rlist)
library(officer)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_n521_uniq-n265_04-03-2024.txt","ALL","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt",0.2) 

base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
sample_batch <- args[5]
resset <- args[6]
dimset=50

outFolder=paste0(base,method,"_pseudobulk_ctrl/")
#for resolution 0.2 and including V2 chem
outFolder=paste0(outFolder,"lessfilt/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

exp <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt")
eigenvec2 <- merge(exp,cov_file,by="dbgap.ID",all.x=T)
#Ali shared that not all waves were run for singlecell and to make sure everything is grabbed correctly this filter is needed
eigenvec2 <- eigenvec2 %>% dplyr::filter(SCAIP7_18 == 1 | SCAIP1_6 == 1)
eigenvec2$Sample_ID <- eigenvec2$dbgap.ID
psychvarstorun <- colnames(eigenvec2[,24:59])
variables <- c("pedu", "pincme", 
                "psesl", 
                "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasg", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq", 
                "genPC1", "genPC2", "genPC3",  "Sex", "cage1", "ceth1", "cwght1", "chght1", "csex1",
                "cgpd5", "cgpd", "cbpd"
                )
 
variable_names <- c("Parental Education", "Parental Income",
                "Subjective SES", 
                "Neighborhood Stress", "Self-disclosure", 
                "Perceived responsiveness", "Parental Warmth", "YR Parent Child Conflict",  "Risky Family", "Youth Depression", 
                "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Stimulated IL-5 cortisol trt", "Stimulated IL-13 cortisol trt", "Stimulated IFNy cortisol trt",    
                "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", 
                "DD Asthma symptoms", "Nightly Asthma", "Asthma Frequency", "Asthma Severity", 
                "Nightly Awakenings",  "Sleep Duration", "Sleep quality",
                "New Genotype PC1", "New Genotype PC2", "New Genotype PC3","Sex_alph", "Age", "Ethnicity", "Weight", "Height", "sex",
                "Female Menarche Status", "Female Puberty Score", "Male puberty score"
                )
variables_df <- data.frame(variable=variables, description=variable_names)

> nrow(cluster_counts_t)
[1] 142636

#trial to inc expression threshold
opfn <- paste0(outFolder,project,".",resset,".",dimset,".wavefilt.singlecellexp.RData")
load(opfn)
library(biomaRt)  
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",1:22), mart = mart)
sum_by <- c("letter_clusters", "BATCH","Sample_ID","treats")
summed <- aggregateAcrossCells(sce, id=colData(sce)[,sum_by])
# Number of cells per sample and cluster
t <- table(colData(sce)$Sample_ID,
           colData(sce)$letter_clusters)
metadata <- as.data.frame(colData(sce))
#/ add rownames using the information from the colData:
colnames(summed) <- apply(colData(summed)[,sum_by], 1, function(x) paste(x, collapse="_"))
raw <- assay(summed, "counts")

counts_ls <- lapply(unique(summed$letter_clusters), function(i){
    cat("running ",i,"\t")
    #i <- unique(summed$letter_clusters)[1]
    cell_idx <- which(tstrsplit(colnames(summed), "_")[[1]] == i)
    keep <- which(colnames(raw) %in% colnames(summed[, cell_idx]))
    data <- raw[, keep]
    #filtered_data <- data > 0 
    filtered_data <- data
    filtered_data[filtered_data < 0] <- NA
    cat(dim(data),"\t")
    if(is.matrix(filtered_data)){
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(data) / 3) #filter to keep genes expressed in 33% of samples
    data <- unlist(data[keep, ])
    cat(dim(data),"\t")
    data <- data[rownames(data) %in% genes$hgnc_symbol, ]
    cat(dim(data),"\n")
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

metadata_ls <- lapply(counts_ls, function(i){
    #i <- counts_ls[[1]]
    cat("running ",unique(i$letter_clusters),"\t")
    df <- data.frame(cluster_sample_id = colnames(i)) ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$letter_clusters <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$BATCH  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    df$Sample_ID  <- tstrsplit(df$cluster_sample_id, "_")[[3]]
    df$treats  <- tstrsplit(df$cluster_sample_id, "_")[[4]]
    test <- metadata %>% dplyr::select(-c(sex_male,sex_female,orig.ident,NEW_BARCODE,percent.mt,nCount_RNA,nFeature_RNA,NUM.READS,NUM.SNPS,RNA_snn_res.0.2))
    #need to remove all RNA_snn_res columns so one per each resolution
    test <- dplyr::select(test, -contains("RNA_snn_res"))
    #test <- test %>% dplyr::select(-c(RNA_snn_res.0.3,RNA_snn_res.0.4))
    test <- unique(transform(test,cluster_sample_id=paste(letter_clusters,BATCH,Sample_ID,treats,sep="_")))
    #test[test$cluster_sample_id=="C0_SCAIP10_AL-021_CTRL",c(1:20)]
    #test <- test[!duplicated(test$cluster_sample_id),]
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

opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt_33.RData")
save(counts_ls,metadata_ls, file=opfn)
#load(opfn)






baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt_33.RData")
load(opfn)

all_cluster_metadata <- ldply(lapply(metadata_ls,function(i){
    i$rows <- rownames(i)
    return(i)
    }), data.frame)[,-1]
rownames(all_cluster_metadata) <- all_cluster_metadata$rows

combatrun="income_PCs_sex_age_and_treats_adjusted_33"
run="treatxcelltype_income_PCs_sex_age_and_treats_adjusted_withWave_33"
combatfolder=paste0(baseoutFolder,combatrun,"/")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleidwave/"))) dir.create(paste0(outFolder,"sampleidwave/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid_all/"))) dir.create(paste0(outFolder,"sampleid_all/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

combatrun="cluster_income_PCs_sex_age_and_treats_adjusted_allclus_33"

all_counts_raw_ls <- lapply(unique(all_cluster_metadata$letter_clusters),function(cluster){
    #cluster=unique(all_cluster_metadata$letter_clusters)[1]
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    genes <- rownames(cluster_counts)
    cluster_counts <- cbind(genes, as.data.frame(cluster_counts))
    return(cluster_counts)
})
all_counts_raw <- ldply(all_counts_raw_ls,data.frame)
rownames(all_counts_raw) <- make.names(all_counts_raw$genes, unique = TRUE)
all_counts_raw[is.na(all_counts_raw)] <- 0
all_counts_raw <- as.matrix(all_counts_raw[,-1])
colnames(all_counts_raw) <- gsub("[.]","-",colnames(all_counts_raw))

cluster_metadata_var <- unique(all_cluster_metadata[,c("Sample_ID","BATCH","treats","Sex","cage1","genPC1","genPC2","genPC3","pincme","letter_clusters")])
cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
cluster_counts_t <- all_counts_raw[,which(colnames(all_counts_raw) %in% rownames(cluster_metadata_var))]
all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
adjusted <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","Sex","cage1","genPC1","genPC2","genPC3","pincme","letter_clusters")])

#took too long to run- timed out. will try running per treat
cluster_metadata_var <- unique(all_cluster_metadata[,c("Sample_ID","BATCH","treats","Sex","cage1","genPC1","genPC2","genPC3","pincme","letter_clusters")])
cluster_metadata_var <- subset(cluster_metadata_var, treats=="CTRL")
cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
cluster_counts_t <- all_counts_raw[,which(colnames(all_counts_raw) %in% rownames(cluster_metadata_var))]
all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

adjusted <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("Sex","cage1","genPC1","genPC2","genPC3","pincme","letter_clusters")])
opfn <- paste0(outFolder,project,".",resset,".",dimset,".","CTRL",".ComBat_seq.",combatrun,".RData")
save(adjusted, file=opfn)

max=unique(rownames(which(adjusted >= .Machine$integer.max, arr.ind = TRUE)))
max

#all_counts_ls <- lapply(unique(all_metadata$letter_clusters),function(cluster){
#    #cluster=unique(all_metadata$letter_clusters)[1]
#    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
#    load(opfn)
#    genes <- rownames(adjusted)
#    adjusted <- cbind(genes, as.data.frame(adjusted))
#    return(adjusted)
#})
#all_counts <- ldply(all_counts_ls,data.frame)
#rownames(all_counts) <- make.names(all_counts$genes, unique = TRUE)
#all_counts[is.na(all_counts)] <- 0
#all_counts <- as.matrix(all_counts[,-1])
#colnames(all_counts) <- gsub("[.]","-",colnames(all_counts))

runage=TRUE
withWave=TRUE
withCOMBAT=TRUE
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
firstrunvars <- psychvarstorun[c(1:2)]

cluster_metadata <- transform(all_cluster_metadata, treats=as.factor(treats))
cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))
cluster_metadata <- transform(cluster_metadata, Sex=as.factor(Sex))
if(all(c("0","1") %in% levels(cluster_metadata$Sex))){
    cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "0"))
}

#for(i in unique(cluster_metadata$treats)){
    for(i in c("CTRL")){
    #i <- "CTRL"
    cluster_metadata_t <- subset(cluster_metadata, treats==i)
    opfn <- paste0(outFolder,project,".",resset,".",dimset,".","CTRL",".ComBat_seq.",combatrun,".RData")
    load(opfn)

    lapply(firstrunvars[c(1:length(firstrunvars))][!firstrunvars[c(1:length(firstrunvars))] %in% puberty],function(var){
        #var="pnsi"
        contrastdf <- data.frame(control_cluster="C0",contrast_cluster=unique(cluster_metadata$letter_clusters)[!unique(cluster_metadata$letter_clusters)=="C0"])
        contrastdf <- transform(contrastdf, contrast1=paste0(var,".","letter_clusters",control_cluster),contrast2=paste0(var,".","letter_clusters",contrast_cluster))
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","letter_clusters",var)]
            mclapply(1:nrow(contrastdf),function(x) {
                #x=1
                con=contrastdf[x,]
                contrast=paste0(con$contrast2,"_vs_",con$contrast1)
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",contrast,".",run,".txt")) > 0)){
            cat("running deseq",var,i,contrast," \n")
            cluster_metadata_var <- subset(cluster_metadata_var, letter_clusters %in% c(con$control_cluster,con$contrast_cluster)) #removing due to low ind counts

        if(withWave==TRUE){
            design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + cage1 + letter_clusters + ", var, " + letter_clusters:",var)
        } else {
            design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + cage1 + letter_clusters + ", var, " + letter_clusters:",var)
        }
        if(runage==TRUE & var=="cage1"){
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","letter_clusters",var)]
            if(withWave==TRUE){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + letter_clusters + ",var, " + letter_clusters:",var)
            } else {
                design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + letter_clusters + ",var, " + letter_clusters:",var)
            }
        }

        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        fwrite(cluster_metadata_var[,c("Sample_ID","Wave")], sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".sampleidwave_",var,"-",i,".",contrast,".",run,".txt"))
        if(withCOMBAT){
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        } else{
        cluster_counts_t <- all_counts_raw[,which(colnames(all_counts_raw) %in% rownames(cluster_metadata_var))]
        }
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                  colData = cluster_metadata_var, 
                  design = as.formula(design))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",var,"-",i,".",contrast,".",run,".RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds, name =con$contrast2)
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$treats =i
        sub.table$contrast=contrast
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",contrast,".",run,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".sigDEGs_",var,"-",i,".",contrast,".",run,".txt"))
        table <- data.frame(symb=var, description=variables_df[variables_df$variable==var,]$description, treats=i,contrast=contrast)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$tested_genes <- paste(nrow(sub.table))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",contrast,".",run,".txt"))
        }
        })
    })
}



for (var in c(psychvarstorun,"cage1")[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]){
cat("running",var,"\n")
myDir <- paste0(outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0(var,".treats"),filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.", "", data_names)
data_names <- gsub("stats_all_cell_types-", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
stats$contrast <- gsub(paste0(var,".treats"),"",stats$contrast)
fwrite(stats, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".stats_all_cell_types-",var,".",run,".treatinteraction.txt"))

myDir <- paste0(outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0(var,".treats"),filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.", "", data_names)
data_names <- gsub("deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
deseqres$contrast <- gsub(paste0(var,".treats"),"",deseqres$contrast)
fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".deseqres_",var,".",run,".treatinteraction.txt"))

myDir <- paste0(outFolder,"sigDEGs/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0(var,".treats"),filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.", "", data_names)
data_names <- gsub("sigDEGs_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
df <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
df$contrast <- gsub(paste0(var,".treats"),"",df$contrast)
fwrite(df, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".sigDEGs_",var,".",run,".treatinteraction.txt"))
}

#Work in progress

myDir <- paste0(outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".treatinteraction.txt",filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.treatinteraction.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#colnames(deseqres) <- c("identifier","padj","pvalue","logFC","SE","var","cluster","contrast")

myDir <- paste0(outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".treatinteraction.txt",filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.treatinteraction.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.stats_all_cell_types-", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
colnames(stats) <- c("symb","description","cluster","contrast","number_samples","number_individuals","gene_number","tested_genes","DEGs_FDR","sigDEGs")
stats <- transform(stats, symb_contrast=paste(symb,contrast,sep=":"))

var50 <- subset(stats, sigDEGs>50)
stats50 <- subset(stats, symb_contrast %in% var50$symb_contrast)
wstats50 <- reshape(stats50[,-c(7,9,11)], idvar = c("symb","description","contrast"), timevar = "cluster", v.names=c("number_individuals","tested_genes","sigDEGs"), direction = "wide",sep=":")

 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

degcols <- grep("sigDEGs",colnames(wstats50))
subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit(colnames(wstats50[,degcols[ apply(wstats50[,degcols],MARGIN=2,FUN=my.max)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(wstats50))]
dft <- subsubvars[,-4] %>% flextable() %>% span_header(sep=":")
cols <- seq(3,36,by=3)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1:2, j = c(4:length(subsubvars[,-4])), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars[,-4])), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(cols), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".treatinteraction.",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("number_individuals|tested_genes",colnames(subsubvars))]
dft <- subsubvars_degonly[,-4] %>% flextable() %>% split_header(sep=":")
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = c(3), border = border, part = "all")
dft <- align(dft, i = 2, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars_degonly[,-4])), align = "center", part = "body")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".treatinteraction.",run,".ndegonly.png"))



    set_caption(caption = "Table 8.1")  




