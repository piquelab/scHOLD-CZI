##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
##library(SeuratDisk) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
##library(SeuratData) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(Signac) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
###library(SeuratWrappers) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
##library(SeuratObject) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(ArchR) ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)


rm(list=ls())


###
### passing arguments
args <- commandArgs(trailingOnly=T)
if ( length(args)>0){
    ###
    batch <- args[1]
    option <- args[2]
    ncpu <- as.integer(args[3])
}else{
    batch <- "HOLD1-ATAC-CTRL"
    option <- "option_nFeature15K"
    ncpu <- 5
}        


###
### set arguments 

outdir <- paste("./3_reCallPeaks_output/", option, "/Peak_matrix_lib/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


addArchRThreads(threads=ncpu)
addArchRGenome("hg38")

dir_ArchR <- "./ArchR_all_output/"


###
### read data 
proj <- loadArchRProject(path=dir_ArchR)


## x <- sort(unique(proj$Cluster2))
## x_df <- data.frame(cl=x, cl_val=as.integer(as.character(gsub("C", "", x))))
## x_df <- x_df%>%dplyr::arrange(cl_val)
## write.table(x_df$cl, "cluster.txt", quote=F, row.names=F, col.names=F)

## x <- getCellColData(proj)
## write.table(sort(unique(x$EXP)), "HOLD.lib.txt", quote=F, row.names=F, col.names=F)
idx_exp <- BiocGenerics::which(proj$EXP==batch)
cellSel <- proj$cellNames[idx_exp]
proj2 <- proj[cellSel,]

cat(batch, length(cellSel), "\n")


###
### get peakMatrix and save 
peak <- getMatrixFromProject(proj2, useMatrix="PeakMatrix")
###
opfn <- paste(outdir, batch, "_peakMatrix.rds", sep="")
write_rds(peak, file=opfn)

###
### END


