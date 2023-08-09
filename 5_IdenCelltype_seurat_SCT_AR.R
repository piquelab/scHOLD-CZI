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
library(SeuratDisk)
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


setwd("/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/")

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

###  
outdir <- "./5b_IdenCelltype_output_renamed/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F) 


future::plan(strategy="multicore", workers=10)
options(future.globals.maxSize=10*20124^3)
plan()



##################### 
### 1 query data  ### 
##################### 

#### run sc transform on the count data prior to normalization
#sc <- read_rds("/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/2_mergeCellRangerAndDemuxlet/seuratObj-merge.md.post-merge-demux.2023-06-23.rds")
#head(sc@meta.data)

#load the renamed seurat object, prior to norm, run sct on it.
sc <- read_rds("/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/2.1_mergeCellRangerAndDemuxlet_renamed/seuratObj-merge.md.post-merge-demux.2023-07-25.rds")
head(sc@meta.data)

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
opfn <- paste0(outdir,"seuratObj-anchors-normSCT-supervPCA-50dim-24kAnchors",Sys.Date(),".rds")
write_rds(anchors, opfn)


## write the sc object
opfn <- paste0(outdir,"seuratObj-sc-transformed-mtregressed-3k-features-5kcell-",Sys.Date(),".rds")
write_rds(sc, opfn)


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
opfn <- paste0(outdir,"seuratObj-sc-annotated-multimodal-ref-mapping-",Sys.Date(),".rds")
write_rds(sc, opfn)



####### Plot





fname=paste0(outdir,"Figure1.1_umap_cellTypes_by_ref_",Sys.Date(),".png");
png(fname,width=7000,height=5000, res=240)
p1 = DimPlot(sc, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 10, repel = TRUE) + NoLegend()
p2 = DimPlot(sc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 10,repel = TRUE) + NoLegend()
p = p1 + p2
print(p)
dev.off()


fname=paste0(outdir,"Figure1.1_umap_cellTypes_by_ref_",Sys.Date(),".pdf");
pdf(fname,width = 20, height = 15)
p1 = DimPlot(sc, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 10, repel = TRUE) + NoLegend()
p2 = DimPlot(sc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 10,repel = TRUE) + NoLegend()
p = p1 + p2
print(p)
dev.off()

# celltype.l1: cell type annotation with granularity level 1 -> 8 categories

# celltype.l2: cell type annotation with granularity level 2 -> 30 clusters



# make initial umap group by 
fname=paste0(outdir,"Figure1.2_umap_cellTypes_by_ref_granLevel1-8cell-types",Sys.Date(),".png");
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
fname=paste0(outdir,"Figure1.3_umap_cellTypes_by_ref_granLevel1-30cell-types",Sys.Date(),".png");
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

fname=paste0(outdir,"Figure2.2_predictionscore_level2",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig1 <- FeaturePlot(sc, features = c("CD4 Naive", "CD14 Mono", "B memory", "CD8 TEM"),  reduction = "ref.umap", 
                    cols = c("lightgrey", "darkred"), ncol = 4) & 
                    theme(plot.title = element_text(size = 13))
print(fig1)
dev.off()


fname=paste0(outdir,"Figure2.1_predictionscore_level1",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig1 <- FeaturePlot(sc, features = c("CD4 T", "Mono", "B", "CD8 T"),  reduction = "ref.umap", 
                    cols = c("lightgrey", "darkred"), ncol = 4) & 
                    theme(plot.title = element_text(size = 13))
print(fig1)
dev.off()



level1 <- sc@meta.data %>% group_by(predicted.celltype.l1) %>% summarise(ncell=n())
fname <- paste0(outdir, "Tabel1_cell_type_count_level_1", Sys.Date(), ".csv")
write.csv(level1, fname, quote=F, row.names=F)


level2 <- sc@meta.data %>% group_by(predicted.celltype.l2) %>% summarise(ncell=n())
fname <- paste0(outdir, "Tabel2_cell_type_count_level_2_30celltypes", Sys.Date(), ".csv")
write.csv(level2, fname, quote=F, row.names=F)




#### are there any cells in the query that are missed when we look at the reference
#merge reference and query
ref$id <- 'reference'
sc$id <- 'query'
refquery <- merge(ref, sc)
refquery[["spca"]] <- merge(ref[["spca"]], sc[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)

refquery$dataset <- "NA"
refquery$dataset[refquery$id=="query"] <- "SCAIP7-18"
refquery$dataset[refquery$id=="reference"] <- "PBMC Ref"

table(refquery$dataset)
table(refquery$id)


fname=paste0(outdir,"Figure3.1_query_vs_ref_leftover_cells",Sys.Date(),".png");
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
opfn <- paste0(outdir,"seuratObj-combined-refCITE-querySCAIP7-18-",Sys.Date(),".rds")
write_rds(refquery, opfn)



#######################################################################
####### assign the cell annotations to the original clusters ##########
#######################################################################

# load the annotated seurat object 
outdir <- "./5b_IdenCelltype_output_renamed/"

## write the sc object, after annotating

opfn <- paste0(outdir,"seuratObj-sc-annotated-multimodal-ref-mapping-2023-07-28.rds")
sc <- read_rds(opfn)


### load the renamed seurat object
## read the sc object
outFolder <- "./2.1_mergeCellRangerAndDemuxlet_renamed/"

opfn <- paste0(outFolder,"seuratObj-post-clustering-res0.5.2023-07-27.rds")
sccr <- read_rds(opfn)


# transfer the annotations from the sc multimodel reference mapped seurat object to the celllranger seurat object with clusters

length(intersect(sc@meta.data$NEW_BARCODE, sccr@meta.data$NEW_BARCODE))
#567386
length(setdiff(sc@meta.data$NEW_BARCODE, sccr@meta.data$NEW_BARCODE))
#3041 # this is the difference in the numb. cells of the two objects. Most likley due to the SCT filtering, less of the cells were removed in MRM object


sc2 <- sccr

columns_to_take <- c("predicted.celltype.l1.score", "predicted.celltype.l1", "predicted.celltype.l2.score", "predicted.celltype.l2")

# Merge the selected columns from dataframe1 into dataframe2 based on the barcode column
sc2@meta.data <- merge(sc2@meta.data, sc@meta.data[, c("NEW_BARCODE", columns_to_take)], by = "NEW_BARCODE", all.x = TRUE)

rownames(sc2@meta.data) <- sc2@meta.data$NEW_BARCODE

## make UMAP for the transfered annotation 

# both levels together
fname=paste0(outdir,"Figure4.1_umap_cellTypes_by_ref_",Sys.Date(),".png");
png(fname,width=7000,height=5000, res=240)
p1 = DimPlot(sc2, reduction = "umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 10, repel = TRUE) + NoLegend()
p2 = DimPlot(sc2, reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 10,repel = TRUE) + NoLegend()
p = p1 + p2
print(p)
dev.off()


# make initial umap group by 
fname=paste0(outdir,"Figure4.2_umap_cellTypes_by_ref_granLevel1-8cell-types",Sys.Date(),".png");
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


library(viridis)

#distinct_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#                     "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
#                     "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
#                     "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
#                     "#8c6d31", "#393b79", "#637939", "#f17cb0", "#843c39",
#                     "#6b6ecf", "#b5cf6b", "#d6616b", "#bd9e39", "#e7969c")

#distinct_colors <- sample(distinct_colors) # which resultedd in this following color set in this order:
#distinct_colors <- c("#8c6d31", "#2ca02c", "#b5cf6b", "#9edae5", "#393b79", "#c49c94", "#9467bd", "#ffbb78", "#17becf", "#ff7f0e", "#e377c2", "#dbdb8d", "#aec7e8", "#98df8a", "#c5b0d5", "#7f7f7f", "#ff9896", "#8c564b", "#bd9e39", "#bcbd22", "#d62728", "#f7b6d2", "#6b6ecf", "#637939", "#1f77b4", "#e7969c", "#843c39", "#f17cb0", "#c7c7c7", "#d6616b")

distinct_colors <- c("#ffbb78", "#c5b0d5", "#e377c2", "#bd9e39", "#9467bd", "#7f7f7f", "#8c564b", "#393b79", "#ff7f0e", "#2ca02c", "#17becf", "#c49c94", "#98df8a", "#6b6ecf", "#bcbd22", "#b5cf6b", "#dbdb8d", "#843c39", "#d62728", "#f17cb0", "#aec7e8", "#9edae5", "#637939", "#e7969c", "#f7b6d2", "#1f77b4", "#c7c7c7", "#8c6d31", "#ff9896", "#d6616b")

sc2@meta.data$predicted.celltype.l2 <- fct_infreq(sc2@meta.data$predicted.celltype.l2)

#quoted_elements <- sprintf("\"%s\"", distinct_colors)
#result <- paste(quoted_elements, collapse = ", ")
#cat(result)

# make initial umap group by 
fname=paste0(outdir,"Figure4.3.2_umap_cellTypes_by_ref_granLevel2-30cell-types",Sys.Date(),".png");
png(fname,width=3500,height=3000, res=240)
fig1 <- DimPlot(sc2, reduction = "umap", label=T, group.by="predicted.celltype.l2", label.size=7.5, raster=F, pt.size = 0.25, repel=T)+ #, )+ #, , pt.size = 0.5, cols=col0)+
        theme_bw()+
        #scale_color_viridis_d(option = "H") +
        #scale_color_brewer(palette = "Set1") +
        scale_color_manual(values = distinct_colors) +
        theme(legend.position="bottom",
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




# prediction score for the new cluster

### prediction score

## transfer prediction score assay data from the annotated seurat object
annot_L1 <- GetAssayData(object =  sc[['prediction.score.celltype.l1']], slot = 'data')
cells_to_take <- colnames(sc2)
annot_L1 <- annot_L1[, cells_to_take]
sc2[["prediction.score.celltype.l1"]] <- CreateAssayObject(data = annot_L1 )

annot_L2 <- GetAssayData(object =  sc[['prediction.score.celltype.l2']], slot = 'data')
cells_to_take <- colnames(sc2)
annot_L2 <- annot_L2[, cells_to_take]
sc2[["prediction.score.celltype.l2"]] <- CreateAssayObject(data = annot_L2 )


###### save the seurat object with transfered assaay data: 
opfn <- paste0(outdir,"seuratObj-MRM-annotations-assays-trasfered-to-SCAIP7-18-object",Sys.Date(),".rds")
write_rds(sc2, opfn)



fname=paste0(outdir,"Figure5.1_predictionscore_level2_",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig1 <- FeaturePlot(sc2, features = c("CD4 Naive", "CD14 Mono", "B memory", "CD8 TEM"),  reduction = "umap", 
                    cols = c("lightgrey", "darkred"), ncol = 4) & 
                    theme(plot.title = element_text(size = 13))
print(fig1)
dev.off()


fname=paste0(outdir,"Figure5.2_predictionscore_level1_",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig1 <- FeaturePlot(sc2, features = c("CD4 T", "Mono", "B", "CD8 T"),  reduction = "umap", 
                    cols = c("lightgrey", "darkred"), ncol = 4) & 
                    theme(plot.title = element_text(size = 13))
print(fig1)
dev.off()



level1 <- sc2@meta.data %>% group_by(predicted.celltype.l1) %>% summarise(ncell=n())
fname <- paste0(outdir, "Tabel1_cell_type_count_level_1_reassigned_", Sys.Date(), ".csv")
write.csv(level1, fname, quote=F, row.names=F)


level2 <- sc2@meta.data %>% group_by(predicted.celltype.l2) %>% summarise(ncell=n())
fname <- paste0(outdir, "Tabel2_cell_type_count_level_2_30celltypes_reassigned_", Sys.Date(), ".csv")
write.csv(level2, fname, quote=F, row.names=F)



#Default(RNA.object) <- "integrated.adt"


# prediction score as violin plots for the cluster

fname=paste0(outdir,"Figure6.1_VlnPlot_predictionscore_level1_",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig0 <- VlnPlot(sc2, features = "predicted.celltype.l1.score", ncol = 1, pt.size = 0, group.by='predicted.celltype.l1')
print(fig0)
dev.off()

# assign the same colors as the UMAP
distinct_colors <- c("#ffbb78", "#c5b0d5", "#e377c2", "#bd9e39", "#9467bd", "#7f7f7f", "#8c564b", "#393b79", "#ff7f0e", "#2ca02c", "#17becf", "#c49c94", "#98df8a", "#6b6ecf", "#bcbd22", "#b5cf6b", "#dbdb8d", "#843c39", "#d62728", "#f17cb0", "#aec7e8", "#9edae5", "#637939", "#e7969c", "#f7b6d2", "#1f77b4", "#c7c7c7", "#8c6d31", "#ff9896", "#d6616b")

fname=paste0(outdir,"Figure6.3_VlnPlot_predictionscore_level2_",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig0 <- VlnPlot(sc2, features = "predicted.celltype.l2.score", ncol = 1, pt.size = 0, cols = distinct_colors, group.by='predicted.celltype.l2')
print(fig0)
dev.off()

# get the prediction scores of the 22 numerical clusters

fname=paste0(outdir,"Figure6.4_VlnPlot_predictionscore_seurat_clusters_",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig0 <- VlnPlot(sc2, features = "predicted.celltype.l2.score", ncol = 1, pt.size = 0, group.by='seurat_clusters')
print(fig0)
dev.off()


# get the prediction scores of the 22 numerical clusters

fname=paste0(outdir,"Figure6.5_VlnPlot_predictionscore_seurat_clusters_",Sys.Date(),".png");
png(fname,width=6000,height=1500, res=240)
fig0 <- VlnPlot(sc2, features = "predicted.celltype.l1.score", ncol = 1, pt.size = 0, group.by='seurat_clusters')
print(fig0)
dev.off()



####################################
############# heat map #############
####################################

######### Level 1: 
# heatmap plus the UMAPs
p1 <- DimPlot(sc2, reduction="umap", group.by="seurat_clusters",
   repel=T, pt.size=0.2, label=T, raster=T, label.size = 7)+
   ggtitle("Seurat clusters")+ 
   theme_bw()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=20),
         axis.title=element_text(size=20),
         axis.text=element_text(size=15))


###
p2 <- DimPlot(sc2, reduction="umap", group.by="predicted.celltype.l1",
   repel=T, pt.size=0.2, label=T,  raster=T, label.size = 7)+
   ggtitle("Multimodal Reference Mapping Annotation")+ 
   theme_bw()+ 
   theme(plot.title=element_text(hjust=0.5, size=20),
         legend.position="none",
         axis.title=element_text(size=20),
         axis.text=element_text(size=15))
         ## legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         ## legend.key.size=grid::unit(1,"lines"))


x <- as.data.frame(sc2@meta.data)

x$predicted.celltype.l1 <- fct_infreq(x$predicted.celltype.l1)

df1 <- x%>%
   group_by(predicted.celltype.l1, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))%>%ungroup()

CL_value <- c("0"=1, "1"=2, "3"=3, "6"=4, "7"=5, "8"=6, "9"=7, "10"=8, "12"=9, "13"=10, "14"=11,
              "2"=12, "4"=13, "5"=14, "11"=15)

df1 <- df1%>%mutate(CL_val=as.numeric(CL_value[as.character(seurat_clusters)]),
                    Cluster2=fct_reorder(seurat_clusters, CL_val))

p3 <- ggplot(df1, aes(x=Cluster2, y=predicted.celltype.l1, fill=Perc))+
   geom_tile()+
   scale_fill_gradient(
      low="#ffffcc", high="#e31a1c", limits=c(0,1), n.breaks=5,
      guide=guide_legend(keywidth=unit(0.5, "cm"),
                         keyheight=unit(1,"cm")))+
   xlab("Cluster")+
   ylab("Annotation Level1")+
   ggtitle("Fraction of cells")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         plot.title=element_text(hjust=0.5, size=20),
         axis.text.x=element_text(hjust=0.5, size=15,margin = margin(r = 15)),
         axis.text.y=element_text(size=20),
         axis.title.x=element_text(size=20),
         axis.title.y=element_text(size=20),
         legend.text=element_text(size=20))

p1 <- as_grob(p1)
p2 <- as_grob(p2)
p3 <- as_grob(p3)
###

figfn <- paste0(outdir, "test-Figure7.3_comb_annot_Level1_",Sys.Date(),".png")
png(figfn, width=2800, height=750, res=120)
plot_grid(p1, p2, p3, nrow=1, ncol=3, labels="AUTO", label_fontface="plain", label_x=0.1,
          align="h", axis="tb", rel_widths=c(1, 1, 1.3))
dev.off()


########## solo heatmap
p3 <- ggplot(df1, aes(x=Cluster2, y=predicted.celltype.l1, fill=Perc))+
   geom_tile()+
   scale_fill_gradient(
      low="#ffffcc", high="#e31a1c", limits=c(0,1), n.breaks=5,
      guide=guide_legend(keywidth=unit(0.5, "cm"),
                         keyheight=unit(1,"cm")))+
   xlab("Cluster")+
   ylab("Annotation Level1")+
   ggtitle("Fraction of cells")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         plot.title=element_text(hjust=0.5, size=25),
         axis.text.x=element_text(hjust=0.5, size=20),
         axis.text.y=element_text(size=20),
         axis.title.x=element_text(size=20),
         axis.title.y=element_text(size=20),
         legend.text=element_text(size=20) 
)


figfn <- paste0(outdir, "Figure7.4_heatmap_comb_annot_Level1_",Sys.Date(),".png")
png(figfn, width=1500, height=700, res=120)
print(p3)
dev.off()




######### Level2 
# heatmap plus the UMAPs
p1 <- DimPlot(sc2, reduction="umap", group.by="seurat_clusters",
   repel=T, pt.size=0.2, label=T, raster=T, label.size = 7)+
   ggtitle("Seurat clusters")+ 
   theme_bw()+
   theme(legend.position="none",
         plot.title=element_text(hjust=0.5, size=20),
         axis.title=element_text(size=20),
         axis.text=element_text(size=15))


###
p2 <- DimPlot(sc2, reduction="umap", group.by="predicted.celltype.l2",
   repel=T, pt.size=0.2, label=T,  raster=T, label.size = 5)+
   ggtitle("Multimodal Reference Mapping Annotation")+ 
   theme_bw()+ 
   theme(plot.title=element_text(hjust=0.5, size=20),
         legend.position="none",
         axis.title=element_text(size=20),
         axis.text=element_text(size=20))
         ## legend.title=element_blank(),
         ## legend.text=element_text(size=8),
         ## legend.key.size=grid::unit(1,"lines"))


x <- as.data.frame(sc2@meta.data)

x$predicted.celltype.l2 <- fct_infreq(x$predicted.celltype.l2)

df1 <- x%>%
   group_by(predicted.celltype.l2, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))%>%ungroup()

CL_value <- c("0"=1, "1"=2, "3"=3, "6"=4, "7"=5, "8"=6, "9"=7, "10"=8, "12"=9, "13"=10, "14"=11,
              "2"=12, "4"=13, "5"=14, "11"=15)

df1 <- df1%>%mutate(CL_val=as.numeric(CL_value[as.character(seurat_clusters)]),
                    Cluster2=fct_reorder(seurat_clusters, CL_val))

p3 <- ggplot(df1, aes(x=Cluster2, y=predicted.celltype.l2, fill=Perc))+
   geom_tile()+
   scale_fill_gradient(
      low="#ffffcc", high="#e31a1c", limits=c(0,1), n.breaks=5,
      guide=guide_legend(keywidth=unit(0.5, "cm"),
                         keyheight=unit(1,"cm")))+
   xlab("Cluster")+
   ylab("Annotation Level2")+
   ggtitle("Fraction of cells")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         plot.title=element_text(hjust=0.5, size=20),
         axis.text.x=element_text(hjust=0.5, size=15),
         axis.text.y=element_text(size=15),
         axis.title.x=element_text(size=20),
         axis.title.y=element_text(size=20),
         legend.text=element_text(size=20))

p1 <- as_grob(p1)
p2 <- as_grob(p2)
p3 <- as_grob(p3)
###

figfn <- paste0(outdir, "Figure7.5_comb_annot_Level2_",Sys.Date(),".png")
png(figfn, width=2800, height=850, res=120)
plot_grid(p1, p2, p3, nrow=1, ncol=3, labels="AUTO", label_fontface="plain", label_x=0.1,
          align="h", axis="tb", rel_widths=c(1, 1, 1.3))
dev.off()



########## solo heatmap
p3 <- ggplot(df1, aes(x=Cluster2, y=predicted.celltype.l2, fill=Perc))+
   geom_tile()+
   scale_fill_gradient(
      low="#ffffcc", high="#e31a1c", limits=c(0,1), n.breaks=5,
      guide=guide_legend(keywidth=unit(0.5, "cm"),
                         keyheight=unit(1,"cm")))+
   xlab("Cluster")+
   ylab("Annotation Level2")+
   ggtitle("Fraction of cells")+ 
   theme_bw()+
   theme(legend.title=element_blank(),
         plot.title=element_text(hjust=0.5, size=25),
         axis.text.x=element_text(hjust=0.5, size=20),
         axis.text.y=element_text(size=20),
         axis.title.x=element_text(size=20),
         axis.title.y=element_text(size=20),
         legend.text=element_text(size=20) 
)


figfn <- paste0(outdir, "Figure7.6_heatmap_comb_annot_Level2_",Sys.Date(),".png")
png(figfn, width=1500, height=1500, res=120)
print(p3)
dev.off()



################################################
################# Dot Plot #####################
################################################


## level 1
x <- as.data.frame(sc2@meta.data)

x$predicted.celltype.l1 <- fct_infreq(x$predicted.celltype.l1)
x$seurat_clusters <- fct_infreq(x$seurat_clusters)

df1 <- x%>%
   group_by(predicted.celltype.l1, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))%>%ungroup()

CL_value <- c("0"=1, "1"=2, "3"=3, "6"=4, "7"=5, "8"=6, "9"=7, "10"=8, "12"=9, "13"=10, "14"=11,
              "2"=12, "4"=13, "5"=14, "11"=15)

df1 <- df1%>%mutate(CL_val=as.numeric(CL_value[as.character(seurat_clusters)]),
                   Cluster2=fct_reorder(seurat_clusters, CL_val))



p <- ggplot(df1, aes(x=Cluster2, y=predicted.celltype.l1))+
     geom_point(aes(fill=Perc), size=4.5,  shape=21, stroke=0)+
     scale_fill_gradient(name="Average expression",
          low="white", high="#de2d26", na.value=NA,
          breaks=waiver(), n.breaks=5, limits=c(0,1),
          guide=guide_colourbar(barwidth=grid::unit(0.8,"lines"),
          barheight=grid::unit(5,"lines")
          )
          )+
     xlab("Cluster")+
     ylab("Annotation Level1")+
     ggtitle("Fraction of cells")+ 
##    scale_size_binned(range=c(0.1, 7), limits=c(0.6,1), guide="none")+
     #scale_color_manual(values=c("gr2"="black"), na.value=NA,  guide="none")+    
     ##scale_size_manual(values=c("a1"=0.1, "a2"=4), guide="none")+ 
     theme(#axis.title=element_blank(),
          ##axis.text.y=element_text(size=9),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          panel.background=element_blank(),
          panel.border=element_rect(color="black", fill=NA, size=1.5),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          ##axis.line=element_line(color="black", size=2),
          legend.title=element_blank())
          #legend.text=element_text(size=8))

##
figfn <- paste0(outdir, "Figure8.1_dotplot_level1.png",Sys.Date(),".png")
png(figfn, width=1500, height=1000, res=240)
print(p)
dev.off()


# Level 2

x <- as.data.frame(sc2@meta.data)

x$predicted.celltype.l2 <- fct_infreq(x$predicted.celltype.l2)
x$seurat_clusters <- fct_infreq(x$seurat_clusters)

df1 <- x%>%
   group_by(predicted.celltype.l2, seurat_clusters)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters)%>%
   mutate(Perc=Freq/sum(Freq))%>%ungroup()

CL_value <- c("0"=1, "1"=2, "3"=3, "6"=4, "7"=5, "8"=6, "9"=7, "10"=8, "12"=9, "13"=10, "14"=11,
              "2"=12, "4"=13, "5"=14, "11"=15)

df1 <- df1%>%mutate(CL_val=as.numeric(CL_value[as.character(seurat_clusters)]),
                   Cluster2=fct_reorder(seurat_clusters, CL_val))



p <- ggplot(df1, aes(x=Cluster2, y=predicted.celltype.l2))+
     geom_point(aes(fill=Perc), size=4.5,  shape=21, stroke=0)+
     scale_fill_gradient(name="Average expression",
          low="white", high="#de2d26", na.value=NA,
          breaks=waiver(), n.breaks=5, limits=c(0,1),
          guide=guide_colourbar(barwidth=grid::unit(0.8,"lines"),
          barheight=grid::unit(5,"lines")))+
     xlab("Cluster")+
     ylab("Annotation Level2")+
     ggtitle("Fraction of cells")+ 
##    scale_size_binned(range=c(0.1, 7), limits=c(0.6,1), guide="none")+
     #scale_color_manual(values=c("gr2"="black"), na.value=NA,  guide="none")+    
     ##scale_size_manual(values=c("a1"=0.1, "a2"=4), guide="none")+ 
     theme(#axis.title=element_blank(),
          ##axis.text.y=element_text(size=9),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          panel.background=element_blank(),
          panel.border=element_rect(color="black", fill=NA, size=1.5),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          ##axis.line=element_line(color="black", size=2),
          legend.title=element_blank())
          #legend.text=element_text(size=8))

##
figfn <- paste0(outdir, "Figure8.2_dotplot_level2.png",Sys.Date(),".png")
png(figfn, width=2000, height=1500, res=240)
print(p)
dev.off()



#############################################################
################# violin per cluster ########################
#############################################################


######### Level 1 annotations 

clustList <- c(0:22)
plotData <- FetchData(sc2,vars= c("seurat_clusters", "predicted.celltype.l1", "predicted.celltype.l1.score"))
plotData$predicted.celltype.l1 <- fct_infreq(plotData$predicted.celltype.l1)

hue_pal()(8)

mycolors= c("CD4 T"="#F8766D", "CD8 T"="#CD9600", "NK"="#7CAE00", "B"="#00BE67", "Mono"="#00BFC4", "other T"="#00A9FF", "other"="#C77CFF", "DC"="#FF61CC")

### plot 0-14
figs_ls <- lapply(0:11, function(i){
    ##
    clust <- clustList[i]
    #plotDF2 <- plotData[,c(i, 5:7)]
    plotDF2 <- plotData[plotData$seurat_clusters==i,]
    #names(plotDF2)[1] <- "y"
    ##
    p <- ggplot(plotDF2, aes(x=predicted.celltype.l1, y=predicted.celltype.l1.score, fill=predicted.celltype.l1))+
       geom_violin(aes(fill=predicted.celltype.l1),width=1)+
       #geom_boxplot(width=0.2,color="grey", outlier.shape=NA)+
       #geom_jitter(width=0.3, size=0.001, shape=20, color="grey88")+ #color="gray",
       geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.000001, alpha = 0.1, shape=20) +
       ylab("Prediction score (L1)")+
       #xlab(paste0("Annotated cell type in cluster ", i))+
       ggtitle(paste0("Cluster ", i))+ 
       scale_fill_manual(values=mycolors)+ 
       #ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p

})


### output figures
figfn <- paste(outdir, "Vln/Figure1.1_violin_clusters_0-11-points_level1.png", sep="")
png(figfn, width=4500, height=2500, res=240)
plot_grid(plotlist=figs_ls, ncol=3)
dev.off()



## plot 15-22
figs_ls <- lapply(12:22, function(i){
    ##
    clust <- clustList[i]
    #plotDF2 <- plotData[,c(i, 5:7)]
    plotDF2 <- plotData[plotData$seurat_clusters==i,]
    #names(plotDF2)[1] <- "y"
    ##
    p <- ggplot(plotDF2, aes(x=predicted.celltype.l1, y=predicted.celltype.l1.score, fill=predicted.celltype.l1))+
       geom_violin(aes(fill=predicted.celltype.l1),width=1)+
       #geom_boxplot(width=0.2,color="grey", outlier.shape=NA)+
       #geom_jitter(width=0.3, size=0.001, shape=20, color="grey88")+ #color="gray",
       geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.000001, alpha = 0.1, shape=20) +
       ylab("Prediction score (L1)")+
       #xlab(paste0("Annotated cell type in cluster ", i))+
       ggtitle(paste0("Cluster ", i))+ 
       scale_fill_manual(values=mycolors)+ 
       #ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p

})


### output figures
figfn <- paste(outdir, "Vln/Figure1.2_violin_clusters_12-22-points_level1.png", sep="")
png(figfn, width=4500, height=2500, res=240)
plot_grid(plotlist=figs_ls, ncol=3)
dev.off()


######### Level 2 annotations 

clustList <- c(0:22)
plotData <- FetchData(sc2,vars= c("seurat_clusters", "predicted.celltype.l2", "predicted.celltype.l2.score"))
plotData$predicted.celltype.l2 <- fct_infreq(plotData$predicted.celltype.l2)

#hue_pal()(30)

mycolors <-  c("CD4 TCM"="#F8766D", "CD4 Naive"="#EF7F49", "CD8 Naive"="#E58700", "NK"="#D89000", "CD8 TEM"="#C99800", "CD14 Mono"="#B79F00", "B naive"="#A3A500",
             "gdT"="#8AAB00", "B intermediate"="#6BB100", "B memory"="#39B600", "MAIT"="#00BA38", "Treg"="#00BD5F", "CD4 TEM"="#00BF7D", "dnT"="#00C097",
             "CD16 Mono"="#00C0AF", "CD8 TCM"="#00BFC4", "NK_CD56bright"="#00BCD8", "Eryth"="#00B7E9", "cDC2"="#00B0F6", "ILC"="#00A7FF", "CD4 CTL"="#619CFF",
             "ASDC"="#9590FF", "Platelet"="#B983FF", "HSPC"="#D376FF", "pDC"="#E76BF3", "Plasmablast"="#F564E3", "NK Proliferating"="#FD61D1", "CD4 Proliferating"="#FF62BC",
             "cDC1"="#FF67A4", "CD8 Proliferating"="#FF6C91")


### plot 0-14
figs_ls <- lapply(0:11, function(i){
    ##
    clust <- clustList[i]
    #plotDF2 <- plotData[,c(i, 5:7)]
    plotDF2 <- plotData[plotData$seurat_clusters==i,]
    #names(plotDF2)[1] <- "y"
    ##
    p <- ggplot(plotDF2, aes(x=predicted.celltype.l2, y=predicted.celltype.l2.score, fill=predicted.celltype.l2))+
       geom_violin(aes(fill=predicted.celltype.l2),width=1)+
       #geom_boxplot(width=0.2,color="grey", outlier.shape=NA)+
       #geom_jitter(width=0.3, size=0.001, shape=20, color="grey88")+ #color="gray",
       geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.000001, alpha = 0.1, shape=20) +
       ylab("Prediction score (l2)")+
       #xlab(paste0("Annotated cell type in cluster ", i))+
       ggtitle(paste0("Cluster ", i))+ 
       scale_fill_manual(values=mycolors)+ 
       #ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p

})


### output figures
figfn <- paste(outdir, "Vln/Figure2.1_violin_clusters_0-11-points_level2.png", sep="")
png(figfn, width=4500, height=2500, res=240)
plot_grid(plotlist=figs_ls, ncol=3)
dev.off()



## plot 15-22
figs_ls <- lapply(12:22, function(i){
    ##
    clust <- clustList[i]
    #plotDF2 <- plotData[,c(i, 5:7)]
    plotDF2 <- plotData[plotData$seurat_clusters==i,]
    #names(plotDF2)[1] <- "y"
    ##
    p <- ggplot(plotDF2, aes(x=predicted.celltype.l2, y=predicted.celltype.l2.score, fill=predicted.celltype.l2))+
       geom_violin(aes(fill=predicted.celltype.l2),width=1)+
       #geom_boxplot(width=0.2,color="grey", outlier.shape=NA)+
       #geom_jitter(width=0.3, size=0.001, shape=20, color="grey88")+ #color="gray",
       geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.000001, alpha = 0.1, shape=20) +
       ylab("Prediction score (l2)")+
       #xlab(paste0("Annotated cell type in cluster ", i))+
       ggtitle(paste0("Cluster ", i))+ 
       scale_fill_manual(values=mycolors)+ 
       #ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p

})


### output figures
figfn <- paste(outdir, "Vln/Figure2.2_violin_clusters_12-22-points_level2.png", sep="")
png(figfn, width=4500, height=2500, res=240)
plot_grid(plotlist=figs_ls, ncol=3)
dev.off()






### plot 0-14
figs_ls <- lapply(0:7, function(i){
    ##
    clust <- clustList[i]
    #plotDF2 <- plotData[,c(i, 5:7)]
    plotDF2 <- plotData[plotData$seurat_clusters==i,]
    #names(plotDF2)[1] <- "y"
    ##
    p <- ggplot(plotDF2, aes(x=predicted.celltype.l2, y=predicted.celltype.l2.score, fill=predicted.celltype.l2))+
       geom_violin(aes(fill=predicted.celltype.l2),width=1)+
       #geom_boxplot(width=0.2,color="grey", outlier.shape=NA)+
       #geom_jitter(width=0.3, size=0.001, shape=20, color="grey88")+ #color="gray",
       geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.000001, alpha = 0.1, shape=20) +
       ylab("Prediction score (l2)")+
       #xlab(paste0("Annotated cell type in cluster ", i))+
       ggtitle(paste0("Cluster ", i))+ 
       scale_fill_manual(values=mycolors)+ 
       #ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p

})


### output figures
figfn <- paste(outdir, "Vln/Figure2.3_violin_clusters_0-7-points_level2.png", sep="")
png(figfn, width=3500, height=2300, res=240)
plot_grid(plotlist=figs_ls, ncol=2)
dev.off()



## plot 8-13
figs_ls <- lapply(8:15, function(i){
    ##
    clust <- clustList[i]
    #plotDF2 <- plotData[,c(i, 5:7)]
    plotDF2 <- plotData[plotData$seurat_clusters==i,]
    #names(plotDF2)[1] <- "y"
    ##
    p <- ggplot(plotDF2, aes(x=predicted.celltype.l2, y=predicted.celltype.l2.score, fill=predicted.celltype.l2))+
       geom_violin(aes(fill=predicted.celltype.l2),width=1)+
       #geom_boxplot(width=0.2,color="grey", outlier.shape=NA)+
       #geom_jitter(width=0.3, size=0.001, shape=20, color="grey88")+ #color="gray",
       geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.000001, alpha = 0.1, shape=20) +
       ylab("Prediction score (l2)")+
       #xlab(paste0("Annotated cell type in cluster ", i))+
       ggtitle(paste0("Cluster ", i))+ 
       scale_fill_manual(values=mycolors)+ 
       #ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p

})


### output figures
figfn <- paste(outdir, "Vln/Figure2.4_violin_clusters_8-15-points_level2.png", sep="")
png(figfn, width=3500, height=2300, res=240)
plot_grid(plotlist=figs_ls, ncol=2)
dev.off()





## plot 14-22
figs_ls <- lapply(16:22, function(i){
    ##
    clust <- clustList[i]
    #plotDF2 <- plotData[,c(i, 5:7)]
    plotDF2 <- plotData[plotData$seurat_clusters==i,]
    #names(plotDF2)[1] <- "y"
    ##
    p <- ggplot(plotDF2, aes(x=predicted.celltype.l2, y=predicted.celltype.l2.score, fill=predicted.celltype.l2))+
       geom_violin(aes(fill=predicted.celltype.l2),width=1)+
       #geom_boxplot(width=0.2,color="grey", outlier.shape=NA)+
       #geom_jitter(width=0.3, size=0.001, shape=20, color="grey88")+ #color="gray",
       geom_point(position = position_jitter(seed = 1, width = 0.3), size=0.000001, alpha = 0.1, shape=20) +
       ylab("Prediction score (l2)")+
       #xlab(paste0("Annotated cell type in cluster ", i))+
       ggtitle(paste0("Cluster ", i))+ 
       scale_fill_manual(values=mycolors)+ 
       #ggtitle(gene)+
       theme_bw()+
       theme(plot.title=element_text(hjust=0, size=12),
             axis.text.x=element_text(size=9, angle=45, hjust=1),
             axis.text.y=element_text(size=9),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size=10),
             legend.position="none")
     p

})


### output figures
figfn <- paste(outdir, "Vln/Figure2.5_violin_clusters_16-22-points_level2.png", sep="")
png(figfn, width=3500, height=2300, res=240)
plot_grid(plotlist=figs_ls, ncol=2)
dev.off()


