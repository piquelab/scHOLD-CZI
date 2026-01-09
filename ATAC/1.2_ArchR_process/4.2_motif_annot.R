##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(SeuratDisk) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(SeuratData) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(Signac) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
###library(SeuratWrappers) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(SeuratObject) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(ArchR) ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(chromVAR)
library(JASPAR2020)
library(JASPAR2022)
library(chromVARmotifs)  ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(motifmatchr)
library(GenomicRanges)

##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)




###
###
rm(list=ls())

outdir <- "./4b_motif_output/motif_anno/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

###
outdir2 <- "./4b_motif_output/motif_anno/motif_list/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)


### Tutorial for files 
### ./4b_motif_output/motif_anno/motif_list/2_motif_***.annotation.txt for each motif include peaks and genes closest to peaks
### ./4b_motif_output/motif_anno/2_comb_cluster_*_motif.annot.txt.gz, combine all the motifs and only include genes infor for each cluster
### ./4b_motif_output/motif_anno/1_motif.list.txt, including 692 motifs, just one column
### ./4b_motif_output/motif_anno/1_motif.infor.txt, motif detailed information for 692 motifs
### ./4b_motif_output/motif_anno/1_peak.annotation.rds, peak annotation to close genes 

###
## extract motif information

## mat <- read_rds("./3_reCallPeaks_output/Peak_matrix/HOLD1-ATAC-CTRL_peakMatrix.rds")
### annotation motif-genes




###
### motif match infor
fn <- "./4b_motif_output/Motif_jaspar2022.match.mat.rds"
motif_match <- read_rds(fn)


###
### peak annotation 
df_row <- rowRanges(motif_match)
anno <- mcols(df_row)
rownames(anno) <- NULL

##  
anno <- anno%>%as.data.frame()%>%
    dplyr::select(gene=nearestGene, distToGeneStart, distToTSS, peakType)

peak2 <- paste(seqnames(df_row), start(df_row), end(df_row), sep="_")
anno$peak <- peak2
anno$chr <- gsub("^chr", "", as.character(seqnames(df_row)))


###
### just run in the first time
opfn <- paste(outdir, "1_peak.annotation.rds", sep="")
write_rds(anno, file=opfn)

## fn <- paste(outdir, "1_peak.annotation.rds", sep="")
## anno <- read_rds(fn)


###
### motif annotation
fn <- "./4b_motif_output/Motif_jaspar2022.motifinfor.ArchR.txt"
motif_df <- read.table(fn, header=T, sep="\t")
motif_df <- motif_df[,c(1:3, 5:6)]


## x <- motif_df%>%drop_na(symbol)


## just run in the first time 
opfn <- paste(outdir, "1_motif.infor.txt", sep="")
write.table(motif_df, file=opfn, row.names=F, quote=F, sep="\t")

opfn <- paste(outdir, "1_motif.list.txt", sep="")
write.table(unique(motif_df$motif_name), file=opfn, row.names=F, quote=F, col.names=F)


## anno <- anno%>%dplyr::filter(chr%in%as.character(1:22))

###
mat2 <- assay(motif_match)
rownames(mat2) <- peak2

cat( "The order of motif", identical(colnames(mat2), motif_df$motif_name), "\n")


 
motif_anno <- NULL
for ( i in 1:nrow(motif_df)){
    ##
    is_peak <- mat2[,i]
    peakSel <- rownames(mat2)[is_peak]

    cat(i, motif_df$motif_name[i], sum(is_peak), "\n")
    
    if ( sum(is_peak) > 0 ){
    ## 100,000        
    anno2 <- anno%>%dplyr::filter(peak%in%peakSel, abs(distToGeneStart)<1e+05)
    
    anno2$motif_name <- motif_df$motif_name[i]
    anno2$motif_name0 <- motif_df$name[i]
    anno2$motif_ID <- motif_df$ID[i]
    anno2$motif_symbol <- motif_df$symbol[i]
    ## anno2$motif_family <- motif_df$family[i]
    ii <- motif_df$motif_name[i]
    opfn <- paste(outdir2, "2_motif_",ii, ".annotation.txt", sep="")
    write.table(anno2, file=opfn, row.names=F, quote=F, sep="\t")    
    }    

}



############################################################################
### combine motif and generate the annotation of TF-genes with all peaks 
#############################################################################

## outdir <- "./4b_motif_output/motif_anno/"
## outdir2 <- "./4b_motif_output/motif_anno/motif_list/"
## fn <- "./4b_motif_output/motif_anno/1_motif.list.txt"
## motif_list <- read.table(fn, header=F)$V1
 
## anno_comb <- map_dfr(motif_list, function(ii){
##     ## 
##     infn_motif <- paste(outdir2, "2_motif_", ii, ".annotation.txt", sep="")
##     anno2 <- read.table(infn_motif, header=T)

##     cat(ii, nrow(anno2), "\n")    
##     ###
##     df0 <- data.frame(gene=unique(anno2$gene))
##     df0$motif_name <- ii
##     df0$motif_name0 <- anno2$motif_name0[1]
##     df0$motif_ID <- anno2$motif_ID[1]
##     df0
## })

## anno_comb <- anno_comb%>%drop_na(gene)
## opfn2 <- paste(outdir, "2_comb_motif.annot.txt.gz", sep="")
## fwrite(anno_comb, file=opfn2, sep="\t", quote=F)


 
#####################################################################################################
### cell type peaks with gene annotation, motif-peak-gene annotation for each cell type peaks (1%)  
######################################################################################################

###
outdir <- "./4b_motif_output/motif_anno/"
outdir2 <- "./4b_motif_output/motif_anno/motif_list/"
fn <- "./4b_motif_output/motif_anno/1_motif.list.txt"
motif_list <- read.table(fn, header=F)$V1


fn_peak <- "./3_reCallPeaks_output/option_nFeature15K/Peak_matrix_cluster/1_celltype_peak.txt.gz"
peak_df <- fread(fn_peak, sep="\t", data.table = F)
 
clusters <- sort(unique(peak_df$Cluster))


###
### 
for ( cl in clusters){
   ###
   ### 
   peak2 <- peak_df%>%filter(Cluster==cl)%>%pull(peak_name)%>%unique()
   ###
   time0 <- Sys.time() 
   anno_comb <- map_dfr(motif_list, function(ii){
     ###
     infn_motif <- paste(outdir2, "2_motif_", ii, ".annotation.txt", sep="")
     anno2 <- read.table(infn_motif, header=T)
     anno2 <- anno2%>%filter(peak%in%peak2)  
       
     ###
     df0 <- data.frame(gene=unique(anno2$gene))
     df0$cluster <- cl  
     df0$motif_name <- ii
     df0$motif_name0 <- anno2$motif_name0[1]
     df0$motif_ID <- anno2$motif_ID[1]
     df0
   })
   time1 <- Sys.time()
   diff <- difftime(time1, time0, units="mins")

   cat(cl, length(motif_list), diff, "\n") 
    
   ###
   ###
   opfn <- paste(outdir, "2_comb_cluster_", cl, "_motif.annot.txt.gz", sep="")
   fwrite(anno_comb, file = opfn, sep = "\t", quote=F)    
} 


    
## anno_comb <- anno_comb%>%drop_na(gene)
## opfn2 <- paste(outdir, "2_comb_motif.celltype.annot.txt.gz", sep="")
## fwrite(anno_comb, file=opfn2, sep="\t")

###
### summary 
fn_peak <- "./3_reCallPeaks_output/option_nFeature15K/Peak_matrix_cluster/1_celltype_peak.txt.gz"
peak_df <- fread(fn_peak, sep="\t", data.table = F)

df_summ <- peak_df%>%dplyr::group_by(Cluster)%>%
    dplyr::summarize(npeak=n(), .groups="drop")%>%dplyr::ungroup()

df_summ <- df_summ%>%mutate(cl_val=as.numeric(gsub("C", "", Cluster)))%>%arrange(cl_val)

##
### summarize number of gene  by motifs
ngene_df <-map_dfr(1:nrow(df_summ), function(i){
   ###
   cl <- df_summ$Cluster[i]
   fn <- paste(outdir, "2_comb_cluster_", cl, "_motif.annot.txt.gz", sep="")
   annot <- fread(fn, sep="\t", data.table=FALSE)
   cat(cl, nrow(annot), "\n")

   ### 
   summ0 <- annot%>%dplyr::group_by(motif_name)%>%dplyr::summarize(ngene=length(unique(gene)), .groups="drop")%>%
       ungroup()
   ###
   summ2 <- data.frame(Cluster=cl, ngene_motif_min=min(summ0$ngene), ngene_motif_max=max(summ0$ngene),
                       ngene_motif_median=median(summ0$ngene))
   summ2
})    

###
### save
df_summ <- df_summ%>%left_join(ngene_df, by="Cluster")

opfn <- paste(outdir, "2_summary.tsv", sep="")
write_tsv(df_summ, file=opfn)

opfn2 <- paste(outdir, "2_summary.xlsx", sep="")
write.xlsx(df_summ, file=opfn2)



##################################################
### hist plots for number of gene by motif 
###################################################

outdir <- "./4b_motif_output/motif_anno/"

fn_ls <- list.files(outdir, "_motif.annot.txt.gz")
summ_all <-map_dfr(fn_ls, function(ii){
   ###
   fn <- paste(outdir, ii, sep="")
   annot <- fread(fn, sep="\t", data.table=FALSE)
    
   cl <- gsub(".*cluster_|_motif.*", "", ii)
   cat(cl, nrow(annot), "\n")    
   ### 
   summ0 <- annot%>%dplyr::group_by(cluster, motif_name)%>%
       dplyr::summarize(ngene=length(unique(gene)), .groups="drop")%>%
       ungroup()
   ### 
})    

summ_all <- summ_all%>%
    dplyr::mutate(cl_val=as.numeric(gsub("^C", "", cluster)),
                  cluster_sort=fct_reorder(cluster, cl_val))

                                     
###
### hist plots
p0 <- ggplot(summ_all, aes(x=ngene))+
    geom_histogram(fill="white", color="grey50")+
    xlab("#genes")+ylab("#motifs")+
    facet_wrap(~cluster_sort, nrow=3, scales="free")+
    theme_bw()+
    theme(axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          strip.text=element_text(size=12))

figfn <- paste(outdir, "Figure2.0_hist.pdf", sep="")
ggsave(figfn, p0, width=9, height=6)


###
### END






###
## fn <- "../2_Differential_analysis/0_pseudobulk.outs/Cluster_res0.07_combat_th20/C0_combat_seq.adjusted.rds"
## x <- read_rds(fn)


###
## fn <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/sigDEGs/ALL.0.2.11.C0.sigDEGs_age-RNA-CTRL.SES_PCs_sex_age_and_treats_generem.txt"
## x <- read.table(fn, header=T)
  


 


###
### TF motif activity matrix 
## motif_mat <- getMatrixFromProject(proj, useMatrix="Motif_JASPAR2022_Matrix")
## opfn <- paste(outdir, "Motif_jaspar2022.activity.mat.rds", sep="")
## write_rds(motif_mat, file=opfn)


              
###
### TF motif match matrix 
## motif_match <- getMatches(proj, name="Motif_JASPAR2022") ### 
## opfn <- paste(outdir, "Motif_jaspar2022.match.mat.rds", sep="")
## write_rds(motif_match, file=opfn)


## df_col <- colData(mat)
## df_row <- rowRanges(mat)
## mat2 <- assay(mat)
## meta <- metadata(mat)

###
### modify meta data and not need 
## proj <- loadArchRProject(path=dir_ArchR)
## x <- getCellColData(proj)
## x2 <- read_rds("./2_Integrate_output/2_atac.meta.rds")
## x2$Cluster2 <- x$Cluster2

## proj$predictedCell_Un <- x2$predictedCell_Un
## proj$predictedGroup_Un <- x2$predictedGroup_Un
## proj$predictedScore_Un <- x2$predictedScore_Un
## ###
## saveArchRProject(proj)
