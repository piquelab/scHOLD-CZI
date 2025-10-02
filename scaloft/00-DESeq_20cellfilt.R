library(DESeq2)
library(qvalue)
library(annotables)
library(tidyr)
library(tidyverse)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(plyr);library(dplyr)
library(parallel)
library(ggseurat)
library(sva)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_uniq-n265_psesl-a2_fixed_12-19-2024.txt","ALL","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt",0.2) 

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
#eigenvec2 <- merge(exp,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
eigenvec2 <- merge(exp,cov_file,by="dbgap.ID",all.x=T)
#Ali shared that not all waves were run for singlecell and to make sure everything is grabbed correctly this filter is needed
eigenvec2 <- eigenvec2 %>% dplyr::filter(SCAIP7_18 == 1 | SCAIP1_6 == 1)
eigenvec2$Sample_ID <- eigenvec2$dbgap.ID
eigenvec2 <- transform(eigenvec2, SCAIP1_6=ifelse(is.na(SCAIP1_6),0,SCAIP1_6)) #weird case where some are NA ... Ali not sure why
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
treatmentsfirst=c("CTRL","PHA")
treatments=c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")

#opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.RData")
opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
#opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")

load(opfn)

#this checks what combos are affected by the filter
cluster="C6"
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  cluster_metadatacell <- cluster_metadata[,c("comb","letter_clusters","cell_count","treats","BATCH","SCAIP1_6")]
cluster_metadatacellctrl <- unique(subset(cluster_metadatacell,treats=="CTRL"))

table(cluster_metadatacellctrl$BATCH)

allmeta <- ldply(lapply(names(counts_ls),function(cluster){
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  cluster_metadatacell <- cluster_metadata[,c("comb","letter_clusters","cell_count","treats","BATCH","SCAIP1_6")]
cluster_metadatacellctrl <- unique(subset(cluster_metadatacell,treats=="CTRL"))
return(cluster_metadatacellctrl)
    }),data.frame)
table(allmeta$BATCH,allmeta$letter_clusters)

allmetaL <- ldply(lapply(names(counts_ls),function(cluster){
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  cluster_metadatacell <- cluster_metadata[cluster_metadata$cell_count<20,c("comb","letter_clusters","cell_count","treats","BATCH","SCAIP1_6")]
cluster_metadatacellctrl <- unique(subset(cluster_metadatacell,treats=="CTRL"))
return(cluster_metadatacellctrl)
    }),data.frame)
table(allmetaL$BATCH,allmetaL$letter_clusters)

allmeta <- ldply(lapply(names(counts_ls),function(cluster){
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  cluster_metadatacell <- cluster_metadata[,c("comb","letter_clusters","cell_count","treats","BATCH","SCAIP1_6")]
cluster_metadatacellctrl <- unique(subset(cluster_metadatacell,treats=="LPS"))
return(cluster_metadatacellctrl)
    }),data.frame)
table(allmeta$BATCH,allmeta$letter_clusters)

allmetaL <- ldply(lapply(names(counts_ls),function(cluster){
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  cluster_metadatacell <- cluster_metadata[cluster_metadata$cell_count<20,c("comb","letter_clusters","cell_count","treats","BATCH","SCAIP1_6")]
cluster_metadatacellctrl <- unique(subset(cluster_metadatacell,treats=="LPS"))
return(cluster_metadatacellctrl)
    }),data.frame)
table(allmetaL$BATCH,allmetaL$letter_clusters)



lcf <- ldply(lapply(names(counts_ls), function(cluster){
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  lowcell <- cluster_metadata[cluster_metadata$cell_count<20,c("comb","letter_clusters","cell_count","treats","BATCH","SCAIP1_6")]
  return(lowcell)
  }), data.frame)
table(lcf$letter_clusters)
table(lcf$letter_clusters,lcf$treats)
lcfctrl <- subset(lcf,treats=="CTRL")
table(lcfctrl$letter_clusters,lcfctrl$SCAIP1_6)

hcf <- ldply(lapply(names(counts_ls), function(cluster){
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  highcell <- cluster_metadata[cluster_metadata$cell_count>=20,c("comb","letter_clusters","cell_count","treats","BATCH","SCAIP1_6")]
  return(highcell)
  }), data.frame)
table(hcf$letter_clusters,hcf$treats)
hcfctrl <- subset(hcf,treats=="CTRL")
table(hcfctrl$letter_clusters,hcfctrl$SCAIP1_6)
hcfctrl6 <- subset(hcfctrl,letter_clusters=="C6")
table(hcfctrl6$BATCH)
png(width = 8, height = 8, file=paste0(baseoutFolder,"figures/",project,".",resset,".",dimset,".C6.ctrl.cellcounts_batch_hist.png"), pointsize=12, 
  bg = "transparent", units = "in", res = 1200)
ggplot(data=unique(hcfctrl6), aes(cell_count)) + 
theme_bw()+
facet_wrap(.~BATCH,ncol=3,scales="free")+
geom_histogram(bins=30)
dev.off()


cellcount <- ldply(lapply(names(counts_ls), function(cluster){
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  cell <- cluster_metadata[cluster_metadata$cell_count,c("comb","letter_clusters","cell_count","treats","BATCH","SCAIP1_6")]
  return(cell)
  }), data.frame)
c6count <- unique(cellcount[cellcount$letter_clusters=="C6",])
png(width = 8, height = 8, file=paste0(baseoutFolder,"figures/",project,".",resset,".",dimset,".cellcounts_cluster_hist.png"), pointsize=12, 
  bg = "transparent", units = "in", res = 1200)
ggplot(data=unique(cellcount[!is.na(cellcount$letter_clusters),]), aes(cell_count)) + 
theme_bw()+
facet_wrap(.~letter_clusters,ncol=3,scales="free")+
geom_histogram(bins=50)
dev.off()

combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
if (!file.exists(baseoutFolder)) dir.create(baseoutFolder, showWarnings=F)
if (!file.exists(paste0(baseoutFolder,"figures/"))) dir.create(paste0(baseoutFolder,"figures/"), showWarnings=F)

mclapply(names(counts_ls), function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata_bf <- data.frame(cluster_metadata_sce)
    highcell <- cluster_metadata_bf[cluster_metadata_bf$cell_count>=20,c("comb","letter_clusters","cell_count")]
    cluster_metadata <- cluster_metadata_bf[rownames(cluster_metadata_bf) %in% rownames(highcell),]
    x <- subset(as.data.frame(table(cluster_metadata$BATCH)),Freq>1)
    cluster_metadata <- subset(cluster_metadata,BATCH %in% unique(x$Var1))
    cluster_metadata_var <- unique(cluster_metadata[,c("Sample_ID","BATCH","treats","Sex","cage1","genPC1","genPC2","genPC3","pincme")])
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
    adjusted <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","Sex","cage1","genPC1","genPC2","genPC3","pincme")])
    #cat(length(rownames(adjusted)),"\n")    
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    save(adjusted, file=opfn)
})

genestoremove <- ldply(lapply(names(counts_ls), function(cluster){
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
max=unique(rownames(which(adjusted >= .Machine$integer.max, arr.ind = TRUE)))
return(data.frame(genes=max))
}),data.frame)
genestoremove
#none so move on

combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
run="income_PCs_sex_age_and_treats_adjusted_withWave"
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

runage=TRUE
withWave=TRUE
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
allvars <- psychvarstorun[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]
firstrunvars <- allvars[c(1:10)]
secondrunvars <- c(allvars[c(1:20)])
for (cluster in names(counts_ls)){
    #cluster="C0"
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    for(i in treatmentsfirst){
            #i <- "CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, Sex=as.factor(Sex))
        if(all(c("0","1") %in% levels(cluster_metadata_t$Sex))){
            cluster_metadata_t <- within(cluster_metadata_t, Sex <- relevel(Sex, ref = "0"))
        }
        if(runage==TRUE){
        if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-age-",i,".",run,".txt")) > 0)){
            var="age"
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1")]
            cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
            cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
            cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
            all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
            if(withWave==TRUE){
            design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + cage1")
            } else {
            design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + cage1")
            }
            cat("running deseq ",var,cluster,i," \n")
            dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                          colData = cluster_metadata_var, 
                                          design = as.formula(design))
            dds <- DESeq(dds,parallel=TRUE)
            opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RDS")
            saveRDS(dds, file=opfn)

            res <- results(dds)
            sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
            names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
            sub.table <- sub.table[!is.na(sub.table$padj), ]
            cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
            sub.table$var=var
            sub.table$cluster=cluster
            sub.table$treats =i
            fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))
            sigDEGs <- subset(sub.table,padj<fdr)
            sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
            table <- data.frame(symb=var, variable= var, cluster=cluster)
            table$number_samples <- paste(nrow(cluster_metadata_var))
            table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
            table$gene_number <- paste(nrow(cluster_counts_t))
            table$DEGs_FDR <- paste(nrow(sigDEGs))
            table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
        }
        }
        mclapply(secondrunvars,function(var){
        #var="pnsi"
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1",var)]
                cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
                cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
                cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
                all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
                if(withWave==TRUE){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + cage1 + ",var)
                } else {
                design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + cage1 + ",var)
                }
                cat("running deseq ",var,cluster,i," \n")
                dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                              colData = cluster_metadata_var, 
                                              design = as.formula(design))
                dds <- DESeq(dds,parallel=TRUE)
                opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RDS")
                saveRDS(dds, file=opfn)

                res <- results(dds)
                sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
                names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
                sub.table <- sub.table[!is.na(sub.table$padj), ]
                cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
                sub.table$var=var
                sub.table$cluster=cluster
                sub.table$treats =i
                fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))
                sigDEGs <- subset(sub.table,padj<fdr)
                sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
                table <- data.frame(symb=var, variable= var, cluster=cluster)
                table$number_samples <- paste(nrow(cluster_metadata_var))
                table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
                table$gene_number <- paste(nrow(cluster_counts_t))
                table$tested_genes <- paste(nrow(sub.table))
                table$DEGs_FDR <- paste(nrow(sigDEGs))
                table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
                fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
            }
        })
    }
}


for (var in c("age",psychvarstorun)[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]){
    #for (var in c(secondrunvars,"age")){
    for (i in treatmentsfirst){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,".sampleidwave_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,".sampleidwave_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".sampleidwave_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".sigDEGs_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    }
}
#filter: ALL.stats_all_cell_types*.nomissvars_PCs*

variables <- c("pedu", "pincme", 
                "psesl", 
                "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasg", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq", 
                "genPC1", "genPC2", "genPC3",  "Sex", "age", "ceth1", "cwght1", "chght1", "csex1",
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

library(flextable)
library(ftExtra)
library(rlist)
library(officer)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)

#[!psychvarstorun %in% puberty]
for (i in treatmentsfirst){
#subvars <- ldply(lapply(psychvarstorun, function(var){
subvars <- ldply(lapply(c(secondrunvars,"age"), function(var){
    cat("running",i,var,"\n")
    if(isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
    #if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"),select=c(1,7))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newls <- newls[sapply(newls, nrow)>0]
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,list.cbind(newls))
            return(newtable)
        }
    }
}),data.frame)

dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=max,na.rm=T)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("numind|ngenes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() 
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = 2, border = border, part = "all")
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(3:length(subsubvars_degonly)), align = "center", part = "body")

#dft <- align(dft, i = 1, j = NULL, align = "left", part = "body")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".ndegonly.png"))
}

puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")[c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds") %in% colnames(cov_file)]
for (cluster in names(counts_ls)){
#for (cluster in c("C8","C9")){
    #cluster="C0"
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))
    cluster_metadata <- transform(cluster_metadata, Sex=as.factor(ifelse(Sex==0,"female","male")))
    if(all(c("male","female") %in% levels(cluster_metadata$Sex))){
        cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "male"))
    }
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
            #run covariates separately for each treatment condition
    for(i in treatmentsfirst){
    #i <- "CTRL"
        for(s in c("male","female")){
        if (s=="female"){
        puberty_sex <- c("cgpd5","cgpd")
            }else {
        puberty_sex <- c("cbpd")
            }
        cluster_metadata_t <- subset(cluster_metadata, treats==i & Sex==s)
            mclapply(puberty_sex, function(var){
        #var="pnsi"
        missing=all(is.na(cluster_metadata_t[,var]))
            if(isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,"_",s,".stats_all_cell_types-",var,"-",i,".",run,".txt")) >0)){
                cat(var,cluster,i,s, "has run already \n")
                next
                } else if (missing){
                    cat("var is all missing information")
                    next
                    } else{
                cat("running")
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3",var)]
                cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
                cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
                cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
                all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
                if(withWave==TRUE){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + ",var)
                } else {
                design <-  paste0("~ genPC1 + genPC2 + genPC3 + ",var)
                }
                cat("running deseq ",var,cluster,i,s," \n")
                dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                              colData = cluster_metadata_var, 
                                              design = as.formula(design))
                dds <- DESeq(dds,parallel=TRUE)
                opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,"_",s,".",run,".RDS")
                saveRDS(dds, file=opfn)

                res <- results(dds)
                sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
                names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
                sub.table <- sub.table[!is.na(sub.table$padj), ]
                cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
                sub.table$var=var
                sub.table$cluster=cluster
                sub.table$treats =i
                sub.table$sex =s
                fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,"_",s,".deseqres_",var,"-",i,".",run,".txt"))
                sigDEGs <- subset(sub.table,padj<fdr)
                sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
                table <- data.frame(symb=var, variable= var, cluster=cluster,sex=s)
                table$number_samples <- paste(nrow(cluster_metadata_var))
                table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
                table$gene_number <- paste(nrow(cluster_counts_t))
                table$tested_genes <- paste(nrow(sub.table))
                table$DEGs_FDR <- paste(nrow(sigDEGs))
                table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
                fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,"_",s,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
                }
            })
        }
    }
}

for (var in puberty){
    #for (var in c(firstrunvars[c(5:10)],"age")){
    for (i in treatmentsfirst){
        for (s in c("male","female")){
        cat("running ",var," ",i,s,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,"_",s,".sampleidwave_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,"_",s,".sampleidwave_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".sampleidwave_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,"_",s,".deseqres_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,"_",s,".deseqres_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,"_",s,".sigDEGs_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,"_",s,".sigDEGs_",var,"-",i,".",run,".txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".sigDEGs_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,"_",s,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,"_",s,".stats_all_cell_types-",var,"-",i,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    }
}
}

zoom="puberty"
for (i in treatmentsfirst){
#subvars <- ldply(lapply(psychvarstorun, function(var){
subvars <- ldply(lapply(puberty, function(var){
    cat("running",i,var,"\n")
    if(isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
    #if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"),select=c(1,7))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newls <- newls[sapply(newls, nrow)>0]
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,list.cbind(newls))
            return(newtable)
        }
    }
}),data.frame)

if(dim(subvars)[1]<1){
    cat("no DEGs>50 \n")
    } else {
dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",zoom,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=max,na.rm=T)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",zoom,".",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("numind|ngenes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() 
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = 2, border = border, part = "all")
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(3:length(subsubvars_degonly)), align = "center", part = "body")

#dft <- align(dft, i = 1, j = NULL, align = "left", part = "body")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",zoom,".",run,".ndegonly.png"))
}
}
