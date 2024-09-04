library(Seurat)
library(future)
library(tidyverse)
library(harmony)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","ALL","fastdemux",0.1) #for testing
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","ALL","demux",0.1)

base <- args[1]
project <- args[2]
method <- args[3]

outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_hg38/",base)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

resset <- as.numeric(args[4])

#future::plan(strategy = 'multicore', workers = 10) #had an issue: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers
#options(future.globals.maxSize = 30 * 1024 ^ 3)

########################

opfn_i <- file.info(dir(outFolder, full.names=T, pattern=paste0(project,".seuratObj-post-umap.")))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

#for(resset in c(0.1,0.15,0.2,0.3,0.4)){
sc <- sc %>% FindClusters(resolution = resset) %>% 
    identity()

opfn <- paste0(outFolder,project,".seuratObj-post-clustering-res",resset,".",Sys.Date(),".rds")
write_rds(sc, opfn)

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern=paste0(project,".seuratObj-post-clustering-res",resset)))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)

# make initial umap group by cluster
fname=paste0(figuredir,project,".Figure5.4_UMAP_Harmony-res",resset,"_group_seurat_cluster_with_names",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=T, group.by="seurat_clusters", label.size=15,pt.size=0.5)+ #, cols=col0)+
  theme(legend.position = "none",legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()

# make initial umap group by treatment
fname=paste0(figuredir,project,".Figure5.1_UMAP_Harmony-res",resset,"_group_treatment",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", group.by = "treats", pt.size = .5)+
        ggtitle("")+
  theme(legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    legend.title=element_blank(),axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()

# make initial umap group by batch
fname=paste0(figuredir,project,".Figure5.2_UMAP_Harmony-res",resset,"_group_batch",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="BATCH",pt.size=0.5)+ #, cols=col0)+
  theme(legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()

# make initial umap group by library
fname=paste0(figuredir,project,".Figure5.3_UMAP_Harmony-res",resset,"_group_Library",Sys.Date(),".png");
png(fname,width=5000,height=4500, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="Library",pt.size=0.5)+ #, cols=col0)+
  theme(legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()

#rownames(sc[["umap"]]@cell.embeddings) <- Cells(sc)
aa <- FetchData(sc,c("umap_1","umap_2","BATCH","EXP","treats","Sample_ID", "seurat_clusters"))

fname=paste0(figuredir,project,".Figure6.1_UMAP_Harmony-res",resset,".grid-batch",Sys.Date(),".png");
png(fname,width=5000,height=4000, res=240)
    p2 <- ggplot(aa,aes(umap_1,umap_2,color=seurat_clusters)) +
                geom_point(size=0.1) +
                facet_wrap(~BATCH) +
  theme(legend.position = "none",legend.key.size = unit(50,"point"),legend.title=element_text(size = rel(1.8)),
  	panel.background = element_rect(fill="white",colour = "black"), strip.text.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
    p2
dev.off()

              
fname=paste0(figuredir,project,".Figure6.2_UMAP_Harmony-res",resset,".grid-treats",Sys.Date(),".png");
png(fname,width=5000,height=3000, res=240)
    p2 <- ggplot(aa,aes(umap_1,umap_2,color=seurat_clusters)) +
                geom_point(size=0.1) +
                facet_wrap(~treats) +
  theme(legend.position = "none",legend.key.size = unit(50,"point"),legend.title=element_text(size = rel(1.8)),
  	panel.background = element_rect(fill="white",colour = "black"), strip.text.x = element_text(size = rel(1.8)),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
    p2
    ##    theme_black()
dev.off()
