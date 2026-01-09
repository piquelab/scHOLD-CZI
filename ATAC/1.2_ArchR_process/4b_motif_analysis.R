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
library(plyranges)
 
##
## library(cowplot)
## library(ggdendro)
## library(RColorBrewer)
## library(ComplexHeatmap)
## library(viridis)
## library(scales)
## library(circlize)
## library(ggrepel)
## library(ggrastr)
## library(openxlsx)


###
###
rm(list=ls())


###
### get motif annotation without X chromosome

outdir <- "./4b_motif_output/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


###
### calculate motif activity

addArchRThreads(threads=16)
addArchRGenome("hg38")


#############################
### Calculate motif 
#############################

###
### 0-read data

dir_ArchR <- "./ArchR_all_output/" 
proj <- loadArchRProject(path=dir_ArchR)


###
### 1-Motif match matrix 
 
### scan motifs from cisdb and jaspar2022
## proj <- addMotifAnnotations(ArchRProj=proj, motifSet="cisbp", annoName="Motif_cisbp")
## proj <- addMotifAnnotations(ArchRProj=proj, motifSet="jaspar2022", annoName="Motif_JASPAR2022")
proj <- addMotifAnnotations(ArchRProj=proj, motifSet="jaspar2022", annoName="Motif_JASPAR2022_excludeX")


###
### 2-calculate motif activity 
 
bg <- getBgdPeaks(proj) 
## bg2 <- subset(bg, subset=seqnames(rowRanges(bg))!="chrX")

proj2 <- addDeviationsMatrix(
    ArchRProj = proj,
    bgdPeaks = bg,
    peakAnnotation = "Motif_JASPAR2022_excludeX",
    matrixName = "Motif_JASPAR2022_excludeX_Matrix")

###
### save 
saveArchRProject(proj2, outputDirectory=dir_ArchR)





#########################################
### get motif-related information
#########################################


outdir <- "./4b_motif_output/"

##
dir_ArchR <- "./ArchR_all_output/"     
proj <- loadArchRProject(path=dir_ArchR)


###
### motif annotation 
infn <- paste(dir_ArchR, "Annotations/Motif_JASPAR2022_excludeX-In-Peaks-Summary.rds", sep="")
motif <- read_rds(infn)

motif_anno <- as.data.frame(motif$motifSummary)
motif_anno <- motif_anno%>%rownames_to_column(var="motif_name")
 
### save output include motif ID information 
opfn <- paste(outdir, "Motif_jaspar2022.motifinfor.ArchR.txt", sep="")
write.table(motif_anno, file=opfn, sep="\t", quote=F, row.names=F, col.names=T)



###
### motif infor from ArchR
## motif_DF <- getVarDeviations(proj, name="Motif_JASPAR2022_Matrix", plot=F)
## opfn <- paste(outdir, "Motif_jaspar2022.motifinfor.S4.rds", sep="")
## write_rds(motif_DF, file = opfn)

 
###
### TF motif activity matrix
### this contains deviations and z-score. We use z-score for downstream analysis, the same to Signac package 

motif_mat <- getMatrixFromProject(proj, useMatrix="Motif_JASPAR2022_excludeX_Matrix") 
opfn <- paste(outdir, "Motif_jaspar2022.activity.mat.rds", sep="")
write_rds(motif_mat, file=opfn)

              
###
### TF motif match matrix 
motif_match <- getMatches(proj, name="Motif_JASPAR2022_excludeX") ### 
opfn <- paste(outdir, "Motif_jaspar2022.match.mat.rds", sep="")
write_rds(motif_match, file=opfn)

###
### END
