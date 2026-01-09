##
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(SeuratDisk) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(SeuratData) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(Signac) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
##library(SeuratWrappers) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(SeuratObject) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
library(ArchR) ###, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(BSgenome.Hsapiens.UCSC.hg38)

##
## library(cowplot)
## library(RColorBrewer)
## library(ComplexHeatmap)
## library(viridis)
## library(circlize)
## library(ggrepel)
## library(ggrastr)
## library(openxlsx)



###
###
rm(list=ls())

option <- "option_nFeature15K"

outdir <- paste("./3_reCallPeaks_output/", option, "/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)
 

###
### use 0.12 resolution

addArchRThreads(threads=16)
addArchRGenome("hg38")

dir_ArchR <- "./ArchR_all_output/"


###
### read data 
proj <- loadArchRProject(path=dir_ArchR)

## x <- as.data.frame(getCellColData(proj))

df_cl <- as.data.frame(table(proj$Clusters))
names(df_cl) <- c("Cluster", "ncell")
df_cl <- df_cl%>%dplyr::arrange(dplyr::desc(ncell))
new_cl <- paste("C", 0:(nrow(df_cl)-1), sep="")
names(new_cl) <- df_cl$Cluster 

###
proj$Cluster2 <- new_cl[proj$Clusters]



###
### option-2, without X chromosome

###
### call peaks
proj2 <- addGroupCoverages(proj, groupBy="Cluster2")
path_Macs2 <- "/wsu/home/groups/piquelab/apps/el7/anaconda3python/envs/macs2/bin/macs2"
proj2 <- addReproduciblePeakSet(
    ArchRProj=proj2,
    groupBy="Cluster2",
    excludeChr=c("chrM", "chrX", "chrY"),
    pathToMacs2=path_Macs2)
 
###
## peakSel <- getPeakSet(proj2)
proj2 <- addPeakMatrix(proj2)


###
### save
dir_ArchR <- "./ArchR_all_output/"
saveArchRProject(proj2, outputDirectory=dir_ArchR)


###
### extract information 
## x <- getCellColData(proj)

## df <- data.frame(EXP=sort(unique(x$EXP)))%>%
##     mutate(EXP_val=as.numeric(gsub("^HOLD|-ATAC-.*", "", EXP)))%>%
##     arrange(EXP_val)
## write.table(df$EXP, "HOLD.lib.txt", quote=F, row.names=F, col.names=F)

## write.table(sort(unique(x$Cluster2)), "cluster.txt", quote=F, row.names=F, col.names=F)

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
