####
####
library(tidyverse)
library(Matrix)
## library(DESeq2)
## library(biobroom)
library(data.table)

##
## library(cowplot)
## library(RColorBrewer)
## library(scales)
## library(viridis)
## library(circlize)
## library(ComplexHeatmap)
## library(openxlsx)
## library(ggrastr)

###

rm(list=ls())

###
### individual motif enrichment test   

outdir <- "./3_regulatory.outs/enrich_motif/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings = F, recursive = T)


###
###
args <- commandArgs(trailingOnly=T)
if ( length(args) > 0){
    ##
    irow <- as.integer(args[1])
}else{
    ##
    irow <- 1
}

cat(irow, "\n")

###############################################
### test which motifs are enriched in DEGs 
###############################################

 
indir <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"
###
### input-1, differential gene results
fn_rna <- paste(indir, "Ali_RNA_psycho.diff.DESeq.txt.gz", sep="")
res_rna <- fread(fn_rna, header=T, data.table=F) 

## fn <- paste(outdir2, "Ali_RNA.DEGs.txt", sep="") 
## DEGs <- read.table(fn, header=F)$V1

###
### input-2, differential TF motifs results   
fn_diff_motif <- paste(indir, "2_psycho_plotData.comb.txt.gz", sep="")
res_motif <- fread(fn_diff_motif, header=T, data.table=F)


###
### input-3, summary cluster from Ali
fn_cl <- "./3_regulatory.outs/list_cluster_var_nDEGs_nTF_new.txt"
df_cl <- read.table(fn_cl, header=T, sep="\t")
df_cl <- df_cl%>%filter(nDEG>0)  ## at least 1 DEGs 


###
### input-4, motif annotation in genes
## fn_annot <- "../1.2_ArchR_process/4_motif_output/motif_anno/2_comb_motif.annot.txt.gz" 
## anno_motif <- fread(fn_annot, header=T, data.table=F)
## gene2 <- anno_motif$gene
dir_annot <- "../1.2_ArchR_process/4b_motif_output/motif_anno/"
fn_ls <- list.files(dir_annot, "2_comb_cluster")
anno_df <- data.frame(cluster=gsub("2_comb_cluster_|.motif.*", "", fn_ls), fn_name=fn_ls)


###
### input-5, motif list
fn_motif <- "./3_regulatory.outs/motif.list.txt"
motif_list <- unique(read.table(fn_motif)$V1)



#####################################
### enrich for each row of df_cl
#####################################


var0 <- df_cl$psycho_variable[irow]
cl_rna <- gsub("^R", "C", df_cl$rna_cluster[irow])
cl_atac <- gsub("^A", "C", df_cl$atac_cluster[irow])

cat("rna:", cl_rna, "atac:", cl_atac, "variable", var0, "\n")  


###
### annotation 
fn0 <- anno_df%>%filter(cluster==cl_atac)%>%pull(fn_name)
anno_fn <- paste(dir_annot, fn0, sep="")
anno_motif <- fread(anno_fn, header=T, sep="\t")
gene2 <- unique(anno_motif$gene)


res0_rna <- res_rna%>%dplyr::filter(cluster==cl_rna, var==var0)
## res0_motif <- res_motif%>%dplyr::filter(Cluster==cl_atac, psycho_variable == var0)

gene_test <- unique(res0_rna$identifier)
##gene_test <- unique(res0_rna$identifier) ## keep the same to overall test 

DEG <- res0_rna%>%dplyr::filter(padj<0.1)%>%pull(identifier)%>%unique()
DEG <- intersect(DEG, gene_test)
notDEG <- setdiff(gene_test, DEG)

df_all <- NULL
for ( ii in motif_list){
###
   anno2 <- anno_motif%>%dplyr::filter(motif_name==ii)
   gene_motif <- intersect(unique(anno2$gene), gene_test)
   ##gene_motif <- unique(anno2$gene)
   ## 
   interest.in <- length(intersect(DEG, gene_motif))
   interest.not <- length(setdiff(DEG, gene_motif))
    
   not.interest.in <- length(intersect(notDEG, gene_motif))
   not.interest.not <- length(setdiff(notDEG, gene_motif))

   if ( length(gene_motif) > 0){    
       ###
       dmat <- matrix(c(interest.in, interest.not, not.interest.in, not.interest.not), 2, 2)
       enrich <- fisher.test(dmat)  ## two side test
       enrich2 <- fisher.test(dmat, alternative="greater") ## one side test
       
       df0 <- data.frame("cluster_rna"=cl_rna, "psycho_variable"=var0, "cluster_atac"=cl_atac, "motif_name"=ii, 
          "interest.in" = interest.in, "interest.not" = interest.not,
          "not.interest.in" = not.interest.in, "not.interest.not" = not.interest.not,
           "odds" = enrich$estimate, "pval" = enrich2$p.value,
           "lower" = enrich$conf.int[1], "upper" = enrich$conf.int[2])

       cat(ii, "\n")
       df_all <- rbind(df_all, df0)
    } ###

} ###

df_all <- df_all%>%
    mutate(log_odds = log(odds), log_lower = log(lower), log_upper = log(upper),
           se = abs(log_odds - log_lower)/1.96)

##
opfn <- paste(outdir, "rna", cl_rna, "_", var0, "_atac", cl_atac, "_enrich.motif.txt", sep="")
write.table(df_all, file = opfn, quote = F, sep = "\t", row.names=F, col.names = T)

###
### END




