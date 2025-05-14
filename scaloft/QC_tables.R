library(DESeq2)
library(plyr)
library(dplyr)
library(parallel)
library(data.table)

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
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
combatrun="income_PCs_sex_age_and_treats_adjusted"
run="income_PCs_sex_age_and_treats_adjusted_withWave"
outFolder=paste0(baseoutFolder,run,"/")
cluster_celltype <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.RData")
load(opfn)

#Num testable genes across clusters
testablegenes <- ldply(lapply(cluster_celltype$cluster, function(cluster){
    #cluster="C0"
        cat("running ", cluster, "\n")
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    testablegenes <- dim(adjusted)[1]
    df <- data.frame(cluster=cluster,genecount=testablegenes)
    return(df)
}), data.frame)
fwrite(testablegenes, sep='\t', quote=F, row.names=F, col.names=T, paste0(baseoutFolder,project,".",resset,".",dimset,".",combatrun,".testablegenes.txt"))

#Number ind per treatment
all_cluster_metadata <- ldply(lapply(metadata_ls,function(i){
    i$rows <- rownames(i)
    return(i)
    }), data.frame)[,-1]
rownames(all_cluster_metadata) <- all_cluster_metadata$rows

indpertreat <- ddply(all_cluster_metadata,"treats",plyr::summarize,
    numind=length(unique(Sample_ID)))
fwrite(indpertreat, sep='\t', quote=F, row.names=F, col.names=T, paste0(baseoutFolder,project,".",resset,".",dimset,".numindpertreat.txt"))

#Num final cells
paste0(baseoutFolder,project,".",resset,".",dimset,".cellcount.txt")

#num ind per cluster in control
all_cluster_metadata_ctrl <- subset(all_cluster_metadata, treats=="CTRL")
indperclus <- ddply(all_cluster_metadata_ctrl,"letter_clusters",plyr::summarize,
    numind=length(unique(Sample_ID)))
fwrite(indperclus, sep='\t', quote=F, row.names=F, col.names=T, paste0(baseoutFolder,project,".",resset,".",dimset,".numindpercluster_CTRL.txt"))

#post combat ind counts
for (i in c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")){
testableind <- ldply(lapply(cluster_celltype$cluster, function(cluster){
    #cluster="C0"
        cat("running ", cluster, i, "\n")
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    testableind <- length(unique(colnames(adjusted)[grep(paste0(i,"$"),colnames(adjusted))]))
    df <- data.frame(cluster=cluster,indcount=testableind)
    return(df)
}), data.frame)
fwrite(testableind, sep='\t', quote=F, row.names=F, col.names=T, paste0(baseoutFolder,project,".",resset,".",dimset,".",combatrun,".postcombat_numindpercluster_",i,".txt"))
}

