### 
library(Matrix)
library(MASS)
## library(scales) 
library(tidyverse)
library(future)
## library(parallel)
library(data.table)
library(GenomicRanges)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
###library(SeuratWrappers)
##library(cicero)
###library(monocle3)
###library(SummarizedExperiment)

###
library(cowplot)
library(ggtext)
library(RColorBrewer)
library(viridis)
library(ggrastr)
library(scales)

rm(list=ls())


###
### all library RNA-seq data 


outdir <- "./2_Integrate_output/scHOLD_RNA_all.outs/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)

plan("multicore", workers=1)
options(future.globals.maxSize=100*1024^3)



########################################
### step-1, RNA-seq normalized data
#########################################


fn <- "/rs/rs_grp_schold/CZI/RNA/analysis/5b_IdenCelltype_fastdemux/noDEX/ALL.seuratObj-.harmony-sctype-0.1.13.rds"
## fn  <- "/rs/rs_grp_schold/CZI/RNA/analysis/5b_IdenCelltype_fastdemux/ALL.seuratObj-.harmony-sctype-2024-05-15.rds"
sc <- read_rds(fn)

### umap object for plots 
opfn <- paste(outdir, "sc_rna.hmn13.umap.rds", sep="")
write_rds(sc[["umap"]], file=opfn)


###
### get meta data
meta <- sc@meta.data
meta <- meta%>%dplyr::select(-seurat_clusters)

cellName <- rownames(meta)

###
### get counts data 
x <- as(object=sc[["RNA"]], Class = "Assay")
count <- x$counts

cat( "The same cell name between the two data", identical(colnames(count), cellName), "\n")

### generate new meta data
prefix_rna <- "/rs/rs_grp_schold/CZI/RNA/analysis/5b_IdenCelltype_fastdemux/noDEX/"

th_ls <- c(0.15, 0.2, 0.3, 0.4)
df_cl <- map_dfc(th_ls, function(th0){
   ###
   ### 
   cat(th0, "read sc data", "\n")
   t0 <- Sys.time() 
   fn0 <- paste(prefix_rna, "ALL.seuratObj-.harmony-sctype-", th0, ".13.rds", sep="")
   sc0 <- read_rds(fn0)
   t1 <- Sys.time()
   ## 
   elapsed <- difftime(t1, t0, units="mins")
   cat("Read done", elapsed, "mins", "\n") 

   ###
   ### 
   x <- sc0@meta.data
   nn <- paste("RNA_snn_res.", th0, sep="")
   x2 <- x%>%dplyr::select(all_of(nn))

   cat(nn, "the row names", identical(rownames(x2), rownames(meta)), "\n")

   if ( ! identical(rownames(x2), rownames(meta))) x2 <- x2[cellName, ]
   x2
})

###
###

cat("The same row names", identical(rownames(meta), rownames(df_cl)), "\n")
meta_new <- cbind(meta, df_cl)


cat("The same cell names", identical(rownames(meta_new), colnames(count)), "\n")
    
###
sc2 <- CreateSeuratObject(counts=count, project="scHOLD", meta.data=meta_new) 
sc2 <- sc2%>%
    NormalizeData()%>%
    ScaleData()%>%
    FindVariableFeatures()
 
sc2 <- sc2%>%RunPCA(npcs=100)

###
### save 
opfn <- paste(outdir, "scHOLD_rna_LogNorm.seurat.rds", sep="")
write_rds(sc2, file=opfn)



###
### add umap infor
fn <- paste(outdir, "scHOLD_rna_LogNorm.seurat.rds", sep="")
sc <- read_rds(fn)
meta <- sc@meta.data

###
###
fn2 <- paste(outdir, "sc_rna.hmn13.umap.rds", sep="")
umap <- read_rds(fn2)
x <- Embeddings(umap)

cat( "the same row names", identical(rownames(meta), rownames(x)), "\n")

meta2 <- cbind(meta, x)

### save
opfn <- paste(outdir, "sc_rna.meta.rds", sep="")
write_rds(meta2, file=opfn)


###
###
#################################
### visulization data 
################################


outdir <- "./2_Integrate_output/scHOLD_RNA_all.outs/"


### input

infn <- paste(outdir, "sc_rna.meta.rds", sep="")
x <- read_rds(infn)



### 1
res_ls <- c(0.1, 0.15, 0.2, 0.3, 0.4)
fig_ls <- lapply(res_ls, function(ii){
   ###
   ### 
   ii2 <- paste("RNA_snn_res.", ii, sep="")
   nnSel <- c(ii2, "umap_1", "umap_2") 
   plotDF <- x%>%dplyr::select(all_of(nnSel))
   names(plotDF)[1] <- "Cluster" 

   cat(ii2, "\n")
    
   plotDF <- plotDF%>%mutate(cl_val = as.numeric(as.character(Cluster)),
       Cluster2 = forcats::fct_reorder(Cluster, cl_val))
   ###
   ### 
   p0 <- ggplot(plotDF, aes(x=umap_1, y=umap_2, colour=Cluster2))+
   rasterise(geom_point(size=0.1), dpi=300)+
   ## facet_wrap(~factor(EXP_sort), ncol=8)+
   guides(col=guide_legend(title=paste("#", ii, sep=""),
          override.aes=list(size=1.2), ncol=1))+
   theme_bw()+
   theme(###strip.text=element_text(size=9),
         legend.title=element_text(size=9),
         ###legend.position.inside=c(0.8, 0.3),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.4,"cm"),
         legend.text=element_text(size=9),
         axis.title=element_text(size=9),
         axis.text=element_text(size=9))
p0

})

### h
## nres <- length(res_ls)
## pcomb <- plot_grid(plotlist=fig_ls, ncol=7) 
## ggsave(figfn, pcomb, width=1800, height=280, units="px", dpi=100)

###
### v
figfn <- paste(outdir, "Figure1.2v_rna_umap.cl.png", sep="")
nres <- length(res_ls)
pcomb <- plot_grid(plotlist=fig_ls, nrow=nres, byrow=F) 
ggsave(figfn, pcomb, width=320, height=1200, units="px", dpi=100)

### pdf 
figfn2 <- paste(outdir, "Figure1.2v_rna_umap.cl.pdf", sep="")
nres <- length(res_ls)
pcomb <- plot_grid(plotlist=fig_ls, nrow=nres, byrow=F) 
ggsave(figfn, pcomb, width=3.2, height=12)




######################
### 0.1 resolution
#######################

outdir <- "./2_Integrate_output/scHOLD_RNA_all.outs/"

col2 <- hue_pal()(10)
col_rna <- c(col2[c(1:5, 9)],"grey")

### input

infn <- paste(outdir, "sc_rna.meta.rds", sep="")
x <- read_rds(infn)

ii <- 0.1
ii2 <- paste("RNA_snn_res.", ii, sep="")
nnSel <- c(ii2, "umap_1", "umap_2") 
plotDF <- x%>%dplyr::select(all_of(nnSel))
names(plotDF)[1] <- "Cluster" 


plotDF <- plotDF%>%mutate(cl_val = as.numeric(Cluster),
       Cluster2 = fct_reorder(Cluster, cl_val))
###
### 
p0 <- ggplot(plotDF, aes(x=umap_1, y=umap_2, colour=Cluster2))+
   rasterise(geom_point(size=0.1), dpi=300)+
   ## facet_wrap(~factor(EXP_sort), ncol=8)+
   scale_color_manual(values = col_rna,
          guide=guide_legend(title=paste("#", ii, sep=""),
          override.aes=list(size=1.2), ncol=1))+
   theme_bw()+
   theme(###strip.text=element_text(size=9),
         legend.title=element_text(size=8),
         ###legend.position.inside=c(0.8, 0.3),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.4,"cm"),
         legend.text=element_text(size=9),
         axis.title=element_text(size=9),
         axis.text=element_text(size=9))

###
figfn <- paste(outdir, "Figure1_res0.1.umap2.png", sep="")
ggsave(figfn, p0, width=500, height=420, units="px", dpi=120)

## pdf
figfn <- paste(outdir, "Figure1_res0.1.umap2.pdf", sep="")
ggsave(figfn, p0, width=5, height=4.20)





###
### END



## ###############################
## ### annotate cell type 
## ###############################


## fn <- paste(outdir, "scHOLD_rna_LogNorm.seurat.rds", sep="")
## query <- read_rds(fn)

## ## query <- SCTransform(query)


## ###
## ### reference data
## fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/pbmc_multimodal.h5seurat"
## ref <- LoadH5Seurat(fn)


## meta <- ref@meta.data
## count <- ref[["SCT"]]$counts

## ref2 <- CreateSeuratObject(counts=count, project="ref-pbmc", meta.data=meta) 
## ref2 <- ref2%>% 
##     NormalizeData()%>%
##     ScaleData()%>%
##     FindVariableFeatures()

## ref2 <- ref2%>%RunPCA(npcs=100)             
             

## ###
## ### integrate RNA with reference 

## anchors <- FindTransferAnchors(reference=ref2, query=query,
##      reference.assay="RNA", reference.reduction="pca", dims=1:50, query.assay ="RNA")

## opfn <- paste(outdir, "2.0_anchors.rds", sep="")
## write_rds(anchors, file=opfn)


## ###
## ### Using predicted.celltype.l1
## pred <- TransferData(anchorset=anchors, refdata=ref2$celltype.l1, dims=1:50)
## opfn <- paste(outdir, "2.1_pred1.rds", sep="")
## write_rds(pred, file=opfn)

## ###
## ### cell-type.l2
## pred2 <- TransferData(anchorset=anchors, refdata=ref2$celltype.l2, dims=1:50)
## opfn <- paste(outdir, "2.1_pred2.rds", sep="")
## write_rds(pred2, file=opfn)


## ### Add celltype in meta.data and update 
## meta <- query@meta.data

## identical(rownames(meta), rownames(pred))
## identical(rownames(meta), rownames(pred2))

## query$predicted.id.l1 <- pred$predicted.id
## query$predicted.id.l2 <- pred2$predicted.id

## meta <- query@meta.data
## table(meta$customclassif, meta$predicted.id.l1)
 
## opfn <- paste(outdir, "scHOLD_rna_LogNorm.seurat.rds", sep="")
## write_rds(query, file=opfn)


## pred <- pred%>%mutate(NEW_BARCODE=rownames(pred))
## opfn <- paste(outdir, "2.1_pred.rds", sep="")
## write_rds(pred, file=opfn)




#######################################
### add umap 
#######################################


## ### input and output file 

## infn_czi2 <- paste(outdir, "scHOLD_rna_LogNorm.seurat.rds", sep="")
## infn_cindy <- "/rs/rs_grp_schold/CZI/RNA/analysis/5b_IdenCelltype_fastdemux/ALL.seuratObj-.harmony-sctype-2024-05-15.rds"

## opfn_czi2 <- infn_czi2


## sc_cindy <- read_rds(infn_cindy)
## ## sc_cindy <- RenameCells(sc_cindy, new.names=gsub("-1", "", Cells(sc_cindy)))

## sc <- read_rds(infn_czi2)
## ## cellSel <- colnames(sc)
 
## ## sc2 <- subset(sc_cindy, cells=cellSel)

## identical(Cells(sc_cindy), Cells(sc))

## sc[["umap"]] <- sc_cindy[["umap"]]
## sc[["harmony"]] <- sc_cindy[["harmony"]]

## write_rds(sc, file=opfn_czi2)



 






##########################################
### visulization scRNA data 
#########################################


## ### input data
## fn  <- "/rs/rs_grp_schold/CZI/RNA/analysis/5b_IdenCelltype_fastdemux/ALL.seuratObj-.harmony-sctype-2024-05-15.rds"
## sc <- read_rds(fn)



## ###
## ### add new column into meta data
## meta <- sc@meta.data
## x <- as.data.table(table(meta$customclassif))
## names(x) <- c("MCls", "ncell")
## x <- x%>%arrange(dplyr::desc(ncell))%>%mutate(MCl_val=(1:nrow(x))-1)

## MCls_new <- c("CD4+", "CD8+", "NK", "Macro", "B", "Mono-CD16+", "DC", "Mono-CD14+")
## names(MCls_new) <- x$MCls
## MCls_val <- x$MCl_val
## names(MCls_val) <- x$MCls

## ###
## meta$MCl2 <- MCls_new[meta$customclassif]
## meta$MCl_value <- MCls_val[meta$customclassif]
## meta <- meta%>%mutate(MCl2_sort=fct_reorder(MCl2, MCl_value))
## sc <- AddMetaData(sc, metadata=meta)


## ###
## ### umap 1
## p1 <- DimPlot(sc, reduction="umap", group.by="seurat_clusters", repel=F,
##    pt.size=0.3, label.size=3, label=T, raster=F)+
##    guides(color=guide_legend(override.aes=list(size=2), ncol=1))+ 
##    ## ggtitle("Seurat clusters (res=0.3)")+ 
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.2,"cm"),
##          ## legend.position="none",
##          plot.title=element_blank(),
##          ## plot.title=element_text(hjust=0.5, size=9),
##          axis.text=element_text(size=9),
##          axis.title=element_text(size=9))


## ###
## ### umap-2
## p2 <- DimPlot(sc, reduction="umap", group.by="MCl2_sort", repel=T,
##    pt.size=0.3, label.size=3, label=T, raster=F)+
##    guides(color=guide_legend(override.aes=list(size=2), ncol=1))+ 
##    ## ggtitle("Seurat clusters (res=0.3)")+ 
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.2,"cm"),
##          ## legend.position="none",
##          plot.title=element_blank(),
##          ## plot.title=element_text(hjust=0.5, size=9),
##          axis.text=element_text(size=9),
##          axis.title=element_text(size=9))

## ###
## pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, rel_widths=c(1, 1.25))

## ### save figures
## figfn <- paste(outdir, "Figure0_scRNA_umap.png", sep="")
## ggsave(figfn, pcomb, width=950, height=400, units="px", dpi=120)


## ### summ statistics

## tmp <- as.data.frame(table(sc$seurat_clusters))
## names(tmp) <- c("seurat_clusters", "ncell")
## opfn <- paste(outdir, "Table1_cl.tsv", sep="")
## write_tsv(tmp, file=opfn)  
## ##
 
## tmp <- as.data.frame(table(sc$MCl2))
## names(tmp) <- c("MCls", "ncell")
## tmp <- tmp%>%arrange(desc(ncell))
## opfn2 <- paste(outdir, "Table1.2_MCls.tsv", sep="")
## write_tsv(tmp, file=opfn2)  
 




## ##################################
## ### summary annottaion results ###
## ##################################

## fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
## atac <- read_rds(fn)

## p0 <- DimPlot(atac, reduction="umap.atac", label=T, raster=F)+
##    theme_bw()## +
##    ## ## ## guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
##    ## theme(legend.title=element_blank(),
##    ##       legend.key.size=grid::unit(0.8,"lines"))
## figfn <- "./4.2_Integrate.outs/Figure1.0_atac.cluster.png"
## png(figfn, width=500, height=500, res=120)
## print(p0)
## dev.off()

## ###
## p1 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.celltype.l1",
##    label=T,  raster=F, repel=T)+
##    theme_bw()+
##    theme(plot.title=element_blank(),
##          legend.title=element_blank(),
##          ## legend.text=element_text(size=8),
##          legend.key.size=grid::unit(1,"lines"))
## figfn <- "./4.2_Integrate.outs/Figure1.1_atac.pred1.png"
## png(figfn, width=500, height=450, res=120)
## print(p1)
## dev.off()

## ###
## p2 <- DimPlot(atac, reduction="umap.atac", group.by="predicted.celltype.l2",
##    label=T, label.size=2.5, raster=F, repel=T)+
##    theme_bw()+
##    theme(plot.title=element_blank(),
##          legend.title=element_blank(),
##          legend.text=element_text(size=8),
##          legend.key.size=grid::unit(0.8,"lines"))

## figfn <- "./4.2_Integrate.outs/Figure1.2_atac.pred2.png"
## png(figfn, width=650, height=400, res=120)
## print(p2)
## dev.off()


## ###
## ### heatmap
## fn <- "./4.2_Integrate.outs/3_scATAC.annot.rds"
## atac <- read_rds(fn)
## meta <- atac@meta.data


## df1 <- meta%>%
##    group_by(predicted.celltype.l1, seurat_clusters)%>%
##    summarize(Freq=n(),.groups="drop")%>%
##    group_by(seurat_clusters)%>%
##    mutate(Perc=Freq/sum(Freq))

## p1 <- ggplot(df1, aes(x=seurat_clusters, y=predicted.celltype.l1, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Fraction of cells",
##       low="#ffffc8", high="#7d0025", na.value=NA)+
##    xlab("Cluster")+ylab("predicted.celltype.l1")+
##    theme_bw()+theme(axis.text.x=element_text(hjust=0.5, size=10))

## ###
## figfn <- "./4.2_Integrate.outs/Figure2.1_atac.heatmap.png"
## png(figfn, width=800, height=600, res=120)
## print(p1)
## dev.off()


## ###
## df2 <- meta%>%
##    group_by(predicted.celltype.l2, seurat_clusters)%>%
##    summarize(Freq=n(),.groups="drop")%>%
##    group_by(seurat_clusters)%>%
##    mutate(Perc=Freq/sum(Freq))

## p2 <- ggplot(df2, aes(x=seurat_clusters, y=predicted.celltype.l2, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Fraction of cells",
##       low="#ffffc8", high="#7d0025", na.value=NA)+
##    xlab("Cluster")+ylab("predicted.celltype.l2")+
##    theme_bw()+theme(axis.text.x=element_text(hjust=0.5))

## ###
## figfn <- "./4.2_Integrate.outs/Figure2.2_atac.heatmap.png"
## png(figfn, width=600, height=600, res=120)
## print(p2)
## dev.off()


####
####
## ref <- LoadH5Seurat("../pbmc_multimodal.h5seurat")

## atac <- read_rds("./4.2_Integrate.outs/1.1_scATAC.cicero.rds")
## DefaultAssay(atac) <- "ACTIVITY"





 
