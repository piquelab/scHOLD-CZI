##
###
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
##
## library(cowplot)
## library(RColorBrewer)
## library(ComplexHeatmap)
## library(viridis)
## library(circlize)
## library(ggrepel)
## library(ggrastr)
## library(openxlsx)

## source("./script/RNAIntegration.R")
rm(list=ls())

args=commandArgs(trailingOnly=T)
if ( length(args)>0){
   ###
   i <- as.integer(args[1])
}else{
   i <- 1
}   
    
###
### use new folder to process data
option <- "harmony_default"
outdir <- paste("./2_Integrate_output/", option, "/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


addArchRThreads(threads=10)
addArchRGenome("hg38")

dir_ArchR <- "./ArchR_all_output/"
infn_RNA <- "./2_Integrate_output/scHOLD_RNA_all.outs/scHOLD_rna_LogNorm.seurat.rds"

opt_harmony <- "Harmony"
opt_umap <- "UMAP"

###
### read data 
proj <- loadArchRProject(path=dir_ArchR)
scRNA <- read_rds(infn_RNA)

### Imputation using MAGIC
## proj <- addImputeWeights(proj, reducedDims=opt_harmony)

## meta <- scRNA@meta.data
## res_atac <- c(0.05, 0.1, 0.15)
## proj2 <- addClusters(input=proj, reducedDims="Harmony", resolution=, force=TRUE)
## proj2 <- proj


###
res_scRNA <- c(0.1, 0.15, 0.2, 0.3, 0.4)


ii <- res_scRNA[i]
ii2 <- paste("RNA_snn_res.", ii, sep="")

cat(ii2, "\n")

proj2 <- addGeneIntegrationMatrix(ArchRProj=proj,
    useMatrix="GeneScoreMatrix",
    matrixName="GeneIntegrationMatrix",
    reducedDims=opt_harmony,
    seRNA=scRNA,
    addToArrow=FALSE,
    plotUMAP=FALSE,
    groupRNA=ii2,
    nameCell="predictedCell_Un",
    nameGroup="predictedGroup_Un",
    nameScore="predictedScore_Un")
 
###
### extract predicted labels 
x <- as.data.frame(getCellColData(proj2))
x2 <- x%>%dplyr::select(predicted.id=predictedGroup_Un, predicted.score=predictedScore_Un)
umap <- getEmbedding(proj2, embedding=opt_umap)
names(umap) <- c("UMAP_1", "UMAP_2")
    
cat("RNA", ii, identical(rownames(x2), rownames(umap)), "\n")

if ( ! identical(rownames(x2), rownames(umap))) umap <- umap[rownames(x2),]


###
### save output
df0 <- cbind(x2, umap)    
opfn <- paste(outdir, "1_RNA.res", ii, ".meta.rds", sep="")
write_rds(df0, file=opfn)    

###
### END



### 
### change cluster with 0.07 resolution and name cluster id by seurat's way
### We don't need run this step and keep current cluster id. For each following analysis, we change cluster id if by this way.
### For Call peaks, we change cluster id in the new cl=oloumn named 'cluster2' in ArchR project.   
## x <- as.data.frame(getCellColData(proj))
## for ( ii in c(0.1, 0.15, 0.2, 0.22, 0.25, 0.3, 0.4, 0.5)){
##     ##
##     fn <- paste(outdir, "1_RNA.res", ii, ".meta.rds", sep="")
##     df0 <- read_rds(fn)
    
##     cat(ii, identical(rownames(x), rownames(df0)), "\n")
    
##     if ( identical(rownames(x), rownames(df0)) ) df0$Clusters <- x$Cluster2

##     ### save
##     opfn <- fn
##     write_rds(df0, file=opfn)
## }    


