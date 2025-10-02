R
library(data.table)
library(qvalue)
library(ggplot2)
library(plyr)
library(tidyverse)

cluster="C4"
treat="CTRL"
PC=2
ci=0.95
resset=0.1
dimset=50
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"

#baseoutFolder=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/")
#opfn <- paste0(baseoutFolder,"ALL.0.2.50.DESeq_countlists_wavefilt.icfilt.RData")
#load(opfn)

#outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
#if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

filter <- "CTRLonly" #ALOFT used
outFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
opfn <- paste0(outFolder,"ALL.",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
load(opfn)

outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/"
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

#############
#loop over clusters
if (!file.exists(paste0(outFolder,"results/"))) dir.create(paste0(outFolder,"results/"), showWarnings=F)

bestPCtable <- lapply(names(counts_ls),function(cluster){
#  lapply(treatments,function(treat){
    cat("running ",cluster,treat,"\n")
myDir <- outFolder #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".cis_qtl.txt.gz", filenames)] #pick specific files from list
filenames <- filenames[grepl(paste0(cluster,".",treat), filenames)] #pick specific files from list
data_names <- gsub("_tensorqtlr.cis_qtl.txt.gz", "", filenames) #remove file ending
PCnames <- paste0("PC",gsub(paste0(cluster,".",treat,".*PC"), "", data_names)) #remove file ending

for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = T)[,PCs:=gsub(".*PC","",data_names[i])]) #read in specific files and set the df object names. can dro0p unwanted columns
#combining data
all_PCs <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(all_PCs) <- PCnames
all_PCs_r <- tapply(all_PCs, names(all_PCs), dplyr::bind_rows)

res <- ldply(lapply(all_PCs_r, function(i){
    df <- data.frame(PCs=as.numeric(unique(i$PCs)),teseted_genes=dim(i)[1],eGenes=sum(i$qval<0.1,na.rm=T))
    return(df)
}), data.frame)
names(res)[1] <- "PCname"

best.PCs <- res[res$eGenes==max(res$eGenes),"PCs"]
best.index <- res[res$eGenes==max(res$eGenes),"PCname"]
res <- res[order(res$PCs),]
if(length(best.PCs)>1){
    best.index <- paste0(cluster,".",treat,".PC1-",1)
    res <- res[order(res$PCs),]
    best.PCs <- best.PCs[1]
}
all_PCs_r_best <- all_PCs_r[[best.index]]

fwrite(res, file=paste0(outFolder,"results/",cluster,".",treat,".eGenes-per-GEPCs.txt"), sep='\t', quote=F, row.names=F)

# save the best results:
fwrite(all_PCs_r_best, paste0(outFolder,"results/",cluster,".",treat,".best_", best.PCs, ".GEPCs.txt"), sep='\t', quote=F, row.names=F)

# subset to significant only:
pc_signif_pairs <- all_PCs_r_best[all_PCs_r_best$qval<0.1,]
pairs <- pc_signif_pairs[,c("phenotype_id","variant_id"),] #geneid and snpid
fwrite(pairs, file=paste0(outFolder,"results/",cluster,".",treat,".",best.index,"_significant_topeeQTL_pairs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
##save SNP IDs
snp_region <- transform(pc_signif_pairs,chr=paste0("chr",sapply(strsplit(variant_id,":"),function(y)y[1])),pos=sapply(strsplit(variant_id,":"),function(y)y[2]))
snp_region_o <- snp_region[,c("chr","pos")] 
fwrite(snp_region_o, file=paste0(outFolder,"results/",cluster,".",treat,".",best.index,"_significant_topeeQTL_snps_region.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
snps <- pc_signif_pairs[,c("variant_id")] 
fwrite(snps, file=paste0(outFolder,"results/",cluster,".",treat,".",best.index,"_significant_topeeQTL_snps.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

all_PCs_r_best$pval_beta[all_PCs_r_best$pval_beta<1e-20] <- 1e-20
pc_results_bp <- all_PCs_r_best %>% select(phenotype_id, variant_id, pval_beta) %>% filter(!is.na(pval_beta)) %>%
            arrange(pval_beta) %>%
            mutate(r=rank(pval_beta, ties.method = "random"),
                   pexp=r/length(pval_beta),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)))

pc_results_bp$group <- "Corrected p-value"
pc_results_bp <- pc_results_bp %>% dplyr::rename(pval = pval_beta) 

png(width = 10, height = 10, file=paste0(outFolder,"figures/",cluster,".",treat,".",best.index,"_eGene_qqplotpermuted_pvalue_clipped.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
    p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(cluster," ",treat," ",best.index," eGene QQ Plot with SNP:eGenes = ",max(res$eGenes))) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()
    return(res[res$eGenes==max(res$eGenes),])
})

names(bestPCtable) <- names(counts_ls)

#######################################################3
##########################################################
#making QC tables
bestPCtableu <- ldply(bestPCtable, data.frame)
best_df <- ldply(lapply(names(counts_ls),function(c){
    cat("running", c, "\n")
    best.PCs <- subset(bestPCtableu, .id==c)$PCs
    #pheno <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/phenotypes.",c,".",treat,".residuals_voom.sort.bed.gz"))
    pheno <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals_ctrlonly/phenotypes.",c,".",treat,".residuals_voom.sort.bed.gz"))
    pc_signif_pairs <- fread(paste0(outFolder,"results/",c,".",treat,".best_", best.PCs, ".GEPCs.txt"))
    df <- data.frame(cluster=c,PCs=best.PCs,numInd=length(colnames(pheno))-4,testedgenes=dim(pc_signif_pairs)[1],eGenes_10=dim(pc_signif_pairs[pc_signif_pairs$qval<0.1,])[1],eGenes_5=dim(pc_signif_pairs[pc_signif_pairs$qval<0.05,])[1])
    return(df)
}),data.frame)
fwrite(best_df, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"bestPCs_table.txt"))

