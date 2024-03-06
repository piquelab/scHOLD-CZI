#preharmony umap plotting. Done in case harmony removes effects of interest

library(Seurat)
library(future)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/") #for testing
base <- args[1]
outFolder=paste0(base,"2.1_mergeCellRangerAndDemuxlet_renamed/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_hg38/",base)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

resset <- as.numeric(args[2])

#################

opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-preharmony-post-umap."))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

sc <- sc %>% FindClusters(resolution = resset) %>% 
    identity()

opfn <- paste0(outFolder,"seuratObj-preharmony-post-clustering-res",resset,".",Sys.Date(),".rds")
write_rds(sc, opfn)

#preharmony umap plotting. Done in case harmony removes effects of interest
# make initial umap group by cluster
fname=paste0(figuredir,"Figure5.4_UMAP_preharmony-res",resset,"_group_seurat_cluster_with_names",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=T, group.by="seurat_clusters", label.size=15,pt.size=0.5)+ #, cols=col0)+
  theme(legend.position = "none",legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()

# make initial umap group by treatment
fname=paste0(figuredir,"Figure5.1_UMAP_preharmony-res",resset,"_group_treatment",Sys.Date(),".png");
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
fname=paste0(figuredir,"Figure5.2_UMAP_preharmony-res",resset,"_group_batch",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="BATCH",pt.size=0.5)+ #, cols=col0)+
  theme(legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()

# make initial umap group by library
fname=paste0(figuredir,"Figure5.3_UMAP_preharmony-res",resset,"_group_Library",Sys.Date(),".png");
png(fname,width=5000,height=4500, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="Library",pt.size=0.5)+ #, cols=col0)+
  theme(legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()