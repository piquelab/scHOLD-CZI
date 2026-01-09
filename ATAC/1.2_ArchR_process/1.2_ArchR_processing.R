###
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
## library(SeuratWrappers)
library(SeuratObject)
library(ArchR)


##
## library(cowplot)
## library(RColorBrewer)
## library(ComplexHeatmap)
## library(viridis)
## library(circlize)
## library(ggrepel)
## library(ggrastr)
## library(openxlsx)


rm(list=ls())


dir_ArchR <- "./ArchR_all_output/"
outdir2 <- "./1_results_all/Batch_correct_Harmony/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

###
addArchRThreads(threads=20)
addArchRGenome("hg38")


################################################
### step-1 create ArchR project 
################################################

####
#### I try parallely submmit the jobs 1.0_createArrowFiles.R and but it failed.
###  So I recommend using loop with parallel in the following way 
#### 

####
### (1) create arrow files 

### cellranger path
basefolder <- "/rs/rs_grp_schold/CZI/ATAC/counts_cellranger_atac/"    

###
### fastdemux clean data
demuxfn <- "./0_demux.outs/2_demux.SNG.correctBatch.rds"
demux <- read_rds(demuxfn) 

####
conditions <- read.table("expNameFiles.txt")$V1  
for ( sample_id in conditions){

###
###    
dir0 <- paste(basefolder, sample_id, "/", sep="")
inputFile <- paste(dir0, "fragments.tsv.gz", sep="")

###
### get validBarcodes
## currently we fill filter after create arrow files not now.
barSel <- demux%>%dplyr::filter(EXP==sample_id)%>%pull(BARCODE)
## barfn <- paste(dir0,  "filtered_peak_bc_matrix/barcodes.tsv", sep="")
## barSel <- read.table(barfn)$V1

###
### create arrow files
time0 <- Sys.time()    
ArrowFiles <- createArrowFiles(inputFiles=inputFile, sampleNames=sample_id,
      excludeChr=c("chrM", "chrY"),
      nChunk = 5,
      validBarcodes=barSel)

time1 <- Sys.time()
elapsed <- difftime(time1, time0, units="mins")

cat("create Arrow file", elapsed, "mins\n")

}


###
### (2), Create ArchRProject

dir_ArchR <- "./ArchR_all_output/"

basefolder <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/"
dir_arrow <- paste(basefolder, "analyses_correct_cellRanger_2025_11_14/1.2_ArchR_process/", sep="")

arrowFiles <- list.files(dir_arrow, "\\.arrow$")  ## 31 libraries 
arrowFiles <- arrowFiles[!grepl("DEX", arrowFiles)]
arrowFiles <- arrowFiles[!grepl("HOLD7", arrowFiles)]

arrowFiles_fullpath <- paste(dir_arrow, arrowFiles, sep="")  ## 26 libraries except DEX and HOLD7 

### 522,094 cells
proj <- ArchRProject(ArrowFiles=arrowFiles_fullpath, outputDirectory=dir_ArchR, copyArrows=TRUE) 

## x <- as.data.frame(getCellColData(proj))
## opfn <- "./1_results_all/Batch_correct_Harmony/0_all.meta.rds"
## write_rds(x, file=opfn)

###
### 3, add meta data and filtering by demuxlet
demuxfn <- "./0_demux.outs/2_demux.SNG.correctBatch.rds"
demux <- read_rds(demuxfn)
demux <- demux%>%mutate(NEW_BARCODE2=paste(EXP, BARCODE, sep="#"), treat=gsub(".*-", "", EXP))


###
### have been subset in this create arrow files 
## idx <- BiocGenerics::which(proj$cellNames%in%demux$NEW_BARCODE2)  ### 522,094 cells
## cellSel <- proj$cellNames[idx]
## proj2 <- subsetArchRProject(ArchRProj=proj, cells=cellSel,
##      outputDirectory=dir_ArchR,
##      dropCells=TRUE, force = TRUE)
     
## proj2 <- proj[cellSel, ]  ## 540,328 cells

proj2 <- proj

###
### add new column to cellColData
cellSel <- proj2$cellNames
demux2 <- demux%>%dplyr::filter(NEW_BARCODE2%in%cellSel)

sampleID <- demux2$BestSample
names(sampleID) <- demux2$NEW_BARCODE2
sampleID2 <- sampleID[cellSel]
proj2$sampleID <- sampleID2
 
proj2$EXP  <- proj2$Sample ##gsub("#.*", "", proj2$cellNames)
proj2$treat <- gsub(".*-", "", proj2$EXP)

###
### save ArchR project
proj2 <- saveArchRProject(ArchRProj=proj2, outputDirectory=dir_ArchR)

###
### save meta data
x2 <- as.data.frame(getCellColData(proj2))
opfn <- "./1_results_all/Batch_correct_Harmony/1_filtered.meta.rds"
write_rds(x2, file=opfn)
 



################################################################
### step 2, dimensionality reduction and clustering 
################################################################

###
## outdir2 <- "./1_results_all/Batch_correct_Harmony/"
dir_ArchR <- "./ArchR_all_output/"

###
proj <- loadArchRProject(dir_ArchR)

###
### (1) Dimensionality reduction, harmony batch correction and  UMAP. try varFeatures=16000, 18000 and 20000 not work

### here We try different resolutions 
## varFeatures=15,000, IterativeLSI, Harmony, UMAP, imputation currectly using harmony, Cluster and Cluster2 (seurat cluster), default
## sampleCellsFinal=300K, IterativeLSI2, Harmony2, UMAP2, Cluster_res0.12
## sampleCellsFinal=350K, IterativeLSI3, Harmony3, UMAP3, Cluster_res0.12_sampleCell350K
## sampleCellsFinal=380K

### option 
## opt_num <- ""
## opt <- "default"
opt_LSI <- "IterativeLSI" ###paste("IterativeLSI", opt_num, sep="")
opt_harmony <- "Harmony"  ### paste("Harmony", opt_num, sep="")
opt_umap <- "UMAP"        ### paste("UMAP", opt_num, sep="")
opt_cl <- "Cluster_res0.12_default"              ### paste("Cluster_res0.12_", opt, sep="")

##### LSI
## proj2 <- addIterativeLSI(ArchRProj=proj, useMatrix="TileMatrix", name=opt_LSI,
##       sampleCellsFinal=380000, projectCellsPre=TRUE, force=TRUE)
proj2 <- addIterativeLSI(ArchRProj=proj, useMatrix="TileMatrix", name=opt_LSI,
      varFeatures=15000, force=TRUE)

###        
#### Harmony
proj2 <- addHarmony(ArchRProj=proj2, reducedDims=opt_LSI, name=opt_harmony, groupBy="EXP", force=TRUE)
#### UMAP
proj2 <- addUMAP(ArchRProj=proj2, reducedDims=opt_harmony, name=opt_umap, force=TRUE)
### cluster 
proj2 <- addClusters(input=proj2, reducedDims=opt_harmony, resolution=0.12, name=opt_cl,
         force=TRUE) 

### Imputation using MAGIC
proj2 <- addImputeWeights(proj2, reducedDims=opt_harmony)
proj2 <- saveArchRProject(ArchRProj=proj2)



###
### (2) loop for different resolution

option <- "Batch_correct_Harmony"   ###"option_estLSI_ncell380K"
outdir2 <- paste("./1_results_all/", option, "/", sep="")
dir.create(outdir2, showWarnings=F, recursive=T)

###
dir_ArchR <- "./ArchR_all_output/"
proj2 <- loadArchRProject(dir_ArchR)

cellSel <- proj2$cellNames
umap <- getEmbedding(proj2, embedding="UMAP")


res_ls <- c(0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.15)
df_cl <- map_dfc(res_ls, function(ii){
###
proj2 <- addClusters(input=proj2, reducedDims="Harmony", resolution=ii, name="Cluster_TMP", force=TRUE)
x <- as.data.frame(getCellColData(proj2, select=c("Cluster_TMP"))) 
names(x) <- paste("ATAC_cluster_res.", ii, sep="")
 
cat("resolution", ii, identical(cellSel, rownames(x)), "\n")

if ( ! identical(cellSel, rownames(x))) x <- x[cellSel,]

x
})

cat(identical(rownames(df_cl), rownames(umap)), "\n")

df2 <- cbind(df_cl, umap)
opfn <- paste(outdir2, "1.0_cluster.rds", sep="")
write_rds(df2, file=opfn)



####
#### (3), run cluster analysis with the defined resolution when we decided which resolution will be used in the final

###
###
dir_ArchR <- "./ArchR_all_output/"
proj <- loadArchRProject(path=dir_ArchR) 

### cluster
proj2 <- addClusters(input=proj, reducedDims="Harmony", resolution=0.12, force=TRUE) 
proj2 <- saveArchRProject(ArchRProj=proj2)






####
### add 0.07 resolution in 1.0_cluster.rds

## dir_ArchR <- "./ArchR_all_output/"
## proj <- loadArchRProject(path=dir_ArchR)
## x <- as.data.frame(getCellColData(proj))
 
## fn <- paste(outdir2, "1.0_cluster.rds", sep="")
## df2 <- read_rds(fn)

## x <- getCellColData(proj)
## identical(rownames(x), rownames(df2))
## df2$ATAC_cluster_res.0.07 <- x$Clusters

## colSort <- sort(colnames(df2))

## df2 <- df2[, colSort]

## opfn <- paste(outdir2, "1.0_cluster.rds", sep="")
## write_rds(df2, file=opfn)

###
### END 



