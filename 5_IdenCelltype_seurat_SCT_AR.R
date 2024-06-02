##
#rm(list=ls())

library("rhdf5")
library("corpcor")
library(Matrix)
library(MASS)
library(scales)
library(tidyverse)
library(parallel)
library(data.table)
library(future)
library(purrr)
library(furrr)
library(Rcpp)
library("BiocParallel")
##
library(Seurat)
library(SeuratDisk) #not installing
library(harmony)
library(annotables)
library(biobroom)
library(org.Hs.eg.db)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(viridis)
theme_set(theme_grey())

library(patchwork)
library(readr)


#####################################################################
### 07/28/2023, Ali R                                          ###### 
###  SCAIP7-18 Multiple Reference Mapping Cell type annotation ######
###  Transfer of annotation and seurat assays                  ######
###  Plots and figure of annotated results                     ######
#####################################################################


###  
args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",0.3,"ALL","fastdemux") #for testing
base <- args[1]
resset <- as.numeric(args[2])
project <- args[3]
method <- args[4]

cat("resolution=",resset,"\nproject=",project,"\n","method=",method,"\n")

basefolder=gsub("analysis/","counts_cellranger_hg38/",base)
outdir=paste0(base,"5b_IdenCelltype_",method,"/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outdir,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy="multicore", workers=10)
options(future.globals.maxSize=10*20124^3)
plan()


##################### 
### 1 query data  ### 
##################### 

#### run sc transform on the count data prior to normalization

#load the renamed seurat object, prior to norm, run sct on it.
opfn_i <- file.info(dir(paste0(base,"2.1_mergeCellRangerAnd",method,"/"), full.names=T, pattern=paste0(project,".seuratObj-merge.md.post-merge-demux.")))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)


#################### previous cell typing ##############

# store mitochondrial percentage in object meta data
sc <- PercentageFeatureSet(sc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
sc <- SCTransform(sc, vars.to.regress = "percent.mt", verbose = FALSE)
# by defualt it runs NormalizeData(), ScaleData(), and FindVariableFeatures().
# nfeatures = 2000 in FindVariableFeatures. 



############################
### load reference data  ### 
############################


### annotation using seurat pbmc reference
ref <- LoadH5Seurat("/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/pbmc_multimodal.h5seurat")


# find anchors with reference.reduction spca (superviside PCA)
anchors <- FindTransferAnchors(reference=ref, query=sc,
           normalization.method="SCT", reference.reduction="spca", dims=1:50)
#found 24 anchors

### write the anchors as an object
opfn <- paste0(outdir,project,".seuratObj-anchors-normSCT-supervPCA-50dim-Anchors",Sys.Date(),".rds")
write_rds(anchors, opfn)


## pred <- TransferData(anchorset=anchors, refdata=ref$celltype.l1, dims=1:30)
## pred <- TransferData(anchorset=anchors, refdata=ref$celltype.l2, dims=1:30)
## opfn <- "./5_IdenCelltype_output/"

sc <- MapQuery(
  anchorset = anchors,
  query = sc,
  reference = ref,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)



## write the sc object, after annotating
opfn <- paste0(outdir,project,".seuratObj-sc-annotated-multimodal-ref-mapping-",Sys.Date(),".rds")
write_rds(sc, opfn)

####### Plot
fname=paste0(figuredir,project,".Figure1.1_umap_cellTypes_by_ref_",Sys.Date(),".png");
png(fname,width=7000,height=5000, res=240)
p1 = DimPlot(sc, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 10, repel = TRUE) + NoLegend()
p2 = DimPlot(sc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 10,repel = TRUE) + NoLegend()
p = p1 + p2
print(p)
dev.off()


fname=paste0(figuredir,project,".Figure1.1_umap_cellTypes_by_ref_",Sys.Date(),".pdf");
pdf(fname,width = 20, height = 15)
p1 = DimPlot(sc, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 10, repel = TRUE) + NoLegend()
p2 = DimPlot(sc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 10,repel = TRUE) + NoLegend()
p = p1 + p2
print(p)
dev.off()

# celltype.l1: cell type annotation with granularity level 1 -> 8 categories

# celltype.l2: cell type annotation with granularity level 2 -> 30 clusters

# make initial umap group by 
fname=paste0(figuredir,project,".Figure1.2_umap_cellTypes_by_ref_granLevel1-8cell-types",Sys.Date(),".png");
png(fname,width=3000,height=3000, res=240)
fig1 <- DimPlot(sc, reduction = "ref.umap", label=T, group.by="predicted.celltype.l1", label.size=10, raster=F, pt.size = 0.25)+ #, )+ #, , pt.size = 0.5, cols=col0)+
        theme_bw()+
        theme(legend.position="right",
              axis.title=element_text(size=25),
              axis.text=element_text(size=25),
              legend.text=element_text(size=25), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(), 
              plot.title=element_text(size=35),  
              text=element_text(size=35),         
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()

# make initial umap group by 
fname=paste0(outdir,project,".Figure1.3_umap_cellTypes_by_ref_granLevel1-30cell-types",Sys.Date(),".png");
png(fname,width=3000,height=3000, res=240)
fig1 <- DimPlot(sc, reduction = "ref.umap", label=T, group.by="predicted.celltype.l2", label.size=7.5, raster=F, pt.size = 0.25, repel=T)+ #, )+ #, , pt.size = 0.5, cols=col0)+
        theme_bw()+
        theme(legend.position="none",
              axis.title=element_text(size=25),
              axis.text=element_text(size=25),
              legend.text=element_text(size=25), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(), 
              plot.title=element_text(size=35),  
              text=element_text(size=35),         
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()

### prediction score

fname=paste0(outdir,project,".Figure2.2_predictionscore_level2",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig1 <- FeaturePlot(sc, features = c("CD4 Naive", "CD14 Mono", "B memory", "CD8 TEM"),  reduction = "ref.umap", 
                    cols = c("lightgrey", "darkred"), ncol = 4) & 
                    theme(plot.title = element_text(size = 13))
print(fig1)
dev.off()


fname=paste0(outdir,project,".Figure2.1_predictionscore_level1",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig1 <- FeaturePlot(sc, features = c("CD4 T", "Mono", "B", "CD8 T"),  reduction = "ref.umap", 
                    cols = c("lightgrey", "darkred"), ncol = 4) & 
                    theme(plot.title = element_text(size = 13))
print(fig1)
dev.off()



level1 <- sc@meta.data %>% group_by(predicted.celltype.l1) %>% summarise(ncell=n())
fname <- paste0(outdir, project,".Tabel1_cell_type_count_level_1", Sys.Date(), ".csv")
write.csv(level1, fname, quote=F, row.names=F)


level2 <- sc@meta.data %>% group_by(predicted.celltype.l2) %>% summarise(ncell=n())
fname <- paste0(outdir, project,".Tabel2_cell_type_count_level_2_30celltypes", Sys.Date(), ".csv")
write.csv(level2, fname, quote=F, row.names=F)




#### are there any cells in the query that are missed when we look at the reference
#merge reference and query
ref$id <- 'reference'
sc$id <- 'query'
refquery <- merge(ref, sc)
refquery[["spca"]] <- merge(ref[["spca"]], sc[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)

refquery$dataset <- "NA"
refquery$dataset[refquery$id=="query"] <- "CZI1"
refquery$dataset[refquery$id=="reference"] <- "PBMC Ref"

table(refquery$dataset)
table(refquery$id)


fname=paste0(outdir,project,".Figure3.1_query_vs_ref_leftover_cells",Sys.Date(),".png");
png(fname,width=5000,height=4000, res=240)
fig1 <- DimPlot(refquery, group.by = 'dataset', shuffle = T, raster=F)+
        theme_bw()+
        theme(legend.position="right",
              axis.title=element_text(size=25),
              axis.text=element_text(size=25),
              legend.text=element_text(size=35), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(), 
              plot.title=element_text(size=45),  
              text=element_text(size=35),         
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()



# write the merged ref and query: 
opfn <- paste0(outdir,project,".seuratObj-combined-refCITE-queryCZI1-",Sys.Date(),".rds")
write_rds(refquery, opfn)



#######################################################################
####### assign the cell annotations to the original clusters ##########
#######################################################################

# load the annotated seurat object 
outdir <- "./5b_IdenCelltype_output_renamed/"

## write the sc object, after annotating
opfn_i <- file.info(dir(outdir, full.names=T, pattern=paste0(project,".seuratObj-sc-annotated-multimodal-ref-mapping")))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)


### load the renamed seurat object
## read the sc object
outFolder <- "./2.1_mergeCellRangerAndDemuxlet_renamed/"

opfn_i <- file.info(dir(outFolder, full.names=T, pattern=paste0(project,".seuratObj-post-clustering-res",resset)))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sccr <- read_rds(opfn)


# transfer the annotations from the sc multimodel reference mapped seurat object to the celllranger seurat object with clusters

length(intersect(sc@meta.data$NEW_BARCODE, sccr@meta.data$NEW_BARCODE))
#160083
length(setdiff(sc@meta.data$NEW_BARCODE, sccr@meta.data$NEW_BARCODE))
#2110 # this is the difference in the numb. cells of the two objects. Most likley due to the SCT filtering, less of the cells were removed in MRM object


sc2 <- sccr

columns_to_take <- c("predicted.celltype.l1.score", "predicted.celltype.l1", "predicted.celltype.l2.score", "predicted.celltype.l2")

# Merge the selected columns from dataframe1 into dataframe2 based on the barcode column
sc2@meta.data <- merge(sc2@meta.data, sc@meta.data[, c("NEW_BARCODE", columns_to_take)], by = "NEW_BARCODE", all.x = TRUE)

rownames(sc2@meta.data) <- sc2@meta.data$NEW_BARCODE

## make UMAP for the transfered annotation 

# both levels together
fname=paste0(figuredir,project,".Figure4.1_umap_cellTypes_by_ref_",Sys.Date(),".png");
png(fname,width=7000,height=5000, res=240)
p1 = DimPlot(sc2, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 10, repel = TRUE) + NoLegend()
p2 = DimPlot(sc2, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 10,repel = TRUE) + NoLegend()
p = p1 + p2
print(p)
dev.off()


# make initial umap group by 
fname=paste0(figuredir,project,".Figure4.2_umap_cellTypes_by_ref_granLevel1-8cell-types",Sys.Date(),".png");
png(fname,width=3000,height=3000, res=240)
fig1 <- DimPlot(sc2, reduction = "umap", label=T, group.by="predicted.celltype.l1", label.size=10, raster=F, pt.size = 0.25)+ #, )+ #, , pt.size = 0.5, cols=col0)+
        theme_bw()+
        theme(legend.position="right",
              axis.title=element_text(size=25),
              axis.text=element_text(size=25),
              legend.text=element_text(size=25), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(), 
              plot.title=element_text(size=35),  
              text=element_text(size=35),         
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()

