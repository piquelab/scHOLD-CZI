library(Seurat)
library(future)
library(tidyverse)
library(harmony)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",10,0.5) #for testing
base <- args[1]
outFolder=paste0(base,"2.1_mergeCellRangerAndDemuxlet_renamed/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_hg38/",base)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

########################

opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-afterPCA."))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

dimset <- as.numeric(args[2])
resset <- as.numeric(args[3])

sc <- RunHarmony(sc,c("Library"),reduction="pca")
sc <- sc %>% RunUMAP(reduction = "harmony", dims = 1:dimset) 
sc <- sc %>% FindNeighbors(reduction = "harmony", dims = 1:dimset) %>% 
    FindClusters(resolution = resset) %>% 
    identity()

opfn <- paste0(outFolder,"seuratObj-post-clustering-res",resset,".",Sys.Date(),".rds")
write_rds(sc, opfn)

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-post-clustering-res"))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)

fname=paste0(figuredir,"UMAP_Harmony-res",resset,".",Sys.Date(),".png")
png(fname,width=5000,height=5000, res=240)
ggdp = DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) #+ NoLegend()
print(ggdp)
dev.off()


# make initial umap group by treatment
fname=paste0(figuredir,"Figure5.1_UMAP_Harmony-res",resset,"_group_treatment",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="treats")+ #, cols=col0)+
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
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()


# make initial umap group by batch
fname=paste0(figuredir,"Figure5.2_UMAP_Harmony-res",resset,"_group_batch",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="BATCH")+ #, cols=col0)+
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
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()

# make initial umap group by library
fname=paste0(figuredir,"Figure5.3_UMAP_Harmony-res",resset,"_group_Library",Sys.Date(),".png");
png(fname,width=5000,height=6000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="Library")+ #, cols=col0)+
        theme_bw()+
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
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()

# make initial umap group by library
fname=paste0(figuredir,"Figure5.4_UMAP_Harmony-res",resset,"_group_seurat_cluster_with_names",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=T, group.by="seurat_clusters", label.size=15)+ #, cols=col0)+
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

aa <- FetchData(sc,c("UMAP_1","UMAP_2","BATCH","EXP","treats","Sample_ID", "seurat_clusters"))

#fname=paste0(outFolder,"UMAP_Harmony-res0.5.grid",".pdf");
#pdf(fname,width=14,height=5)
fname=paste0(figuredir,"Figure6.1_UMAP_Harmony-res",resset,".grid-batch",Sys.Date(),".png");
png(fname,width=5000,height=4000, res=240)
    #fname=paste0(outFolder,"UMAP_LocationHarmony.Origin.png");
    #png(fname,width=1600,height=1200)
    p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
                geom_point(size=0.1) +
                ##    scale_color_manual(values=group.colors) +
##                guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
##                scale_color_manual("Origin",values=c("M"="#D1D1D1","F"="#A61BB5"))+
                facet_wrap(~BATCH) +
                theme_bw() +
                theme(legend.position="right",
                axis.title=element_text(size=25),
                axis.text=element_text(size=25),
                legend.text=element_text(size=25), 
                strip.text.x = element_text(size = 30),
                #legend.position=c(0.1, 0.85),
                legend.title=element_text(size=25),
                #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
                legend.key=element_blank(), #legend.key=element_rect(fill=NA),
                #legend.box.background=element_blank(), 
                plot.title=element_text(size=35),           
                panel.border=element_rect(colour="black", fill=NA))
                #guides(shape = guide_legend(override.aes = list(size = 35)))
    p2
    ##    theme_black()
dev.off()

#fname=paste0(outFolder,"UMAP_Harmony-res0.5.grid",".pdf");
#pdf(fname,width=14,height=5)
fname=paste0(figuredir,"Figure6.2_UMAP_Harmony-res",resset,".grid-treats",Sys.Date(),".png");
png(fname,width=5000,height=4000, res=240)
    #fname=paste0(outFolder,"UMAP_LocationHarmony.Origin.png");
    #png(fname,width=1600,height=1200)
    p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
                geom_point(size=0.1) +
                ##    scale_color_manual(values=group.colors) +
##                guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
##                scale_color_manual("Origin",values=c("M"="#D1D1D1","F"="#A61BB5"))+
                facet_wrap(~treats) +
                theme_bw() +
                theme(legend.position="right",
                axis.title=element_text(size=25),
                axis.text=element_text(size=25),
                legend.text=element_text(size=25), 
                strip.text.x = element_text(size = 30),
                #legend.position=c(0.1, 0.85),
                legend.title=element_text(size=25),
                #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
                legend.key=element_blank(), #legend.key=element_rect(fill=NA),
                #legend.box.background=element_blank(), 
                plot.title=element_text(size=35),           
                panel.border=element_rect(colour="black", fill=NA))
                #guides(shape = guide_legend(override.aes = list(size = 35)))
    p2
    ##    theme_black()
dev.off()
