library(Seurat)
library(future)
library(tidyverse)
library(harmony)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","ALL","fastdemux",0.2,11) #for testing
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","ALL","demux",0.2,50)

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
dimset <- as.numeric(args[5])

#future::plan(strategy = 'multicore', workers = 10) #had an issue: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers
#options(future.globals.maxSize = 30 * 1024 ^ 3)

########################
opfn <- paste0(outFolder,project,".seuratObj-post-umap.",dimset,".rds")
sc <- read_rds(opfn)

#CZI load in is this:
opfn_i <- file.info(dir(paste0(base,"2.1_mergeCellRangerAnd",method,"/"), full.names=T, pattern=paste0(project,".seuratObj-post-umap.2024-05-08.rds")))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

cat("runing", dimset, resset, "\n")

sc <- subset(sc,subset=treats != "RNA-LPS-DEX") #2 singletons identified. 7 final clusters.
scres <- sc %>% FindClusters(resolution = resset) %>% 
    identity()

#opfn <- paste0(outFolder,project,".seuratObj-post-clustering-res",resset,".",dimset,".rds")
opfn <- paste0(outFolder,project,".seuratObj-post-clustering-nodex-res",resset,".",dimset,".rds")
write_rds(scres, opfn)

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern=paste0(project,".seuratObj-post-clustering-res",resset)))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)

# make initial umap group by cluster
fname=paste0(figuredir,project,".Figure5.4_UMAP_Harmony-nodex-res",resset,".",dimset,"_group_seurat_cluster_with_names",".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(scres, reduction = "umap", label=T, group.by="seurat_clusters", label.size=15,pt.size=0.5)+ #, cols=col0)+
  theme(legend.position = "none",legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()

####

library(tidyverse)
library(data.table)
library(future)
library(HGNChelper)
library(Seurat)
library(plyr);library(dplyr)

####### cell typing 
args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",0.2,"ALL","fastdemux",11) #for testing
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/",0.2,"ALL","demux",13) 

base <- args[1]
resset <- as.numeric(args[2])
project <- args[3]
method <- args[4]
dimset=args[5]
cat("resolution=",resset,"\nproject=",project,"\n","method=",method,"\n")

outdir=paste0(base,"5b_IdenCelltype_",method,"/")
#outdir=paste0(outdir,"indwavecellcountfilt/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outdir,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

harmonyfolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
#harmonyfolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/indwavecellcountfilt/")

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

opfn <- paste0(harmonyfolder,project,".seuratObj-post-clustering-nodex-res",resset,".",dimset,".rds")
sc <- read_rds(opfn)

cat("starting sctype \n")
  es.max = sctype_score(scRNAseqData = sc[["RNA"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)  #sc[["RNA"]]@scale.data for seurat V<5
  cL_resutls = do.call("rbind", lapply(unique(sc@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sc@meta.data[sc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% dplyr::group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  sc@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sc@meta.data$customclassif[sc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  png(width = 9, height = 8, file=paste0(figuredir,project,".",resset,".",dimset,".harmony_nodex_umap_QC_wSCtypecelltype.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
  p <- DimPlot(sc, reduction = "umap", label = TRUE, repel=TRUE, group.by = 'customclassif', pt.size = .1)
  print(p)
  dev.off()
  cat("saving \n")
opfn <- paste0(outdir,project,".seuratObj-.harmony_nodex-sctype-",resset,".",dimset,".rds")
write_rds(sc, opfn)

#plot which barcodes were assigned to the same cluster as previosuly
sc_wdex <- read_rds(paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds"))
sc_wdex@meta.data$letter_clusters <- paste0("C",sc_wdex@meta.data$seurat_clusters)
#subset for 20cell count filt
counts <- plyr::count(sc_wdex@meta.data,c("Library","letter_clusters","Sample_ID"))
names(counts)[4] <- "btic_cellcounts"
sc_wdex@meta.data <- left_join(sc_wdex@meta.data,counts,by=c("Library","letter_clusters","Sample_ID"))
rownames(sc_wdex@meta.data) <- sc_wdex@meta.data$NEW_BARCODE
sc_wdex <- subset(sc_wdex, subset=btic_cellcounts>=20)
wdex <- data.frame(barcode=rownames(sc_wdex@meta.data), seurat_clusters=sc_wdex@meta.data$seurat_clusters)
rm(counts)
gc(reset=T)

sc <- read_rds(paste0(nodexoutdir,project,".seuratObj-.harmony_nodex-sctype-",resset,".",dimset,".rds"))
sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
#subset for 20cell count filt
counts <- plyr::count(sc@meta.data,c("Library","letter_clusters","Sample_ID"))
names(counts)[4] <- "btic_cellcounts"
sc@meta.data <- left_join(sc@meta.data,counts,by=c("Library","letter_clusters","Sample_ID"))
rm(counts)
rownames(sc@meta.data) <- sc@meta.data$NEW_BARCODE
sc <- subset(sc, subset=btic_cellcounts>=20)
nodex <- data.frame(barcode=rownames(sc@meta.data), seurat_clusters=sc@meta.data$seurat_clusters)

scheck <- merge(wdex,nodex,by="barcode")
scheck <- subset(scheck, seurat_clusters.x %in% c("0","1","2","3","4","5") & seurat_clusters.y %in% c("0","1","2","3","4","5"))
identical(scheck$seurat_clusters.x,scheck$seurat_clusters.y) #false as expected
scheck <- transform(scheck, matching_type=ifelse(as.character(seurat_clusters.x)==as.character(seurat_clusters.y),"matching","not_matching"))
table(scheck$matching_type)
    matching not_matching
      174424       381775
fwrite(scheck,file=paste0(outdir,"scheck.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#plot bar, match cluster vs not match cluster
celltype <- data.frame(cluster=c("0","1","2","3","4","5"),celltype=c("CD4_T","CD4_T","CD8_T","NK","Monocyte","B"))
#heatmap : rows=barcode, col=old and new cluster
scheck <- fread(paste0(outdir,"scheck.txt"))
scheck <- transform(scheck, seurat_clusters.y=ifelse(seurat_clusters.y==4,as.integer(3),ifelse(seurat_clusters.y==3,as.integer(4),seurat_clusters.y)))
scheckm <- melt(scheck)
scheckm <- transform(scheckm, cluster_assigned=paste0(variable,"_",value))
scheckm$cluster_assigned <- reorder(scheckm$cluster_assigned, scheckm$value)
scheckm$barcode <- reorder(scheckm$barcode, scheckm$value)

    p <- ggplot(scheckm, aes(x = barcode, y = cluster_assigned)) +
      geom_tile() +
      #scale_fill_gradient(low = "white", high = "steelblue") + # Adjust colors as needed
      theme_minimal()

  png(width = 9, height = 8, file=paste0(figuredir,"dex_vs_nodex_clusters.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
print(p)
dev.off()

scheckl <- ldply(lapply(split(scheck,scheck$seurat_clusters.x),function(c){
  xclus=unique(c$seurat_clusters.x)
  cdf <- subset(celltype, cluster==xclus)
  yc <- plyr::count(c,c("seurat_clusters.y"))
  yc <- merge(yc, celltype, by.x="seurat_clusters.y",by.y="cluster")
  yc <- transform(yc, matchingcluster=ifelse(seurat_clusters.y==cdf$cluster,"match","not_match"),
    matchingcelltype=ifelse(celltype==cdf$celltype,"match","not_match"))
  mcluster_ratio=sum(yc[yc$matchingcluster=="match",]$freq,na.rm=T)/sum(yc$freq,na.rm=T)
  mcell_ratio=sum(yc[yc$matchingcelltype=="match",]$freq,na.rm=T)/sum(yc$freq,na.rm=T)
  return(data.frame(cdf,mcluster_ratio,mcell_ratio)) #cluster and celltype from old clustering (ie.including dex)
  }),data.frame)[,-1]

schecklm <- melt(scheckl,id.vars=c("cluster","celltype"))
p <- ggplot(schecklm,aes(x=paste0(cluster,"_",celltype),y=value,fill=factor(variable)))+
  geom_bar(stat="identity",position="dodge")+
      theme_minimal()
  png(width = 9, height = 8, file=paste0(figuredir,"dex_vs_nodex_ratios.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
print(p)
dev.off()

scheck_count <- data.frame(dexcount=plyr::count(scheck,"seurat_clusters.x"),nodexcount=plyr::count(scheck,"seurat_clusters.y"))
scheck_count <- transform(scheck_count, perc=nodexcount.freq/(nodexcount.freq+dexcount.freq))
    p <- ggplot(scheck_count, aes(x = dexcount.seurat_clusters.x, y = nodexcount.seurat_clusters.y,fill=perc)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high =  "steelblue", limits = c(0, 1))+
      labs(fill="percent without dex cells",x="with dex clusters",y="without dex clusters")+
      scale_x_continuous(breaks = scheck_count$dexcount.seurat_clusters.x)+
      scale_y_continuous(breaks = scheck_count$nodexcount.seurat_clusters.y)+
      theme_minimal()
  png(width = 9, height = 8, file=paste0(figuredir,"dex_vs_nodex_clusters_heat.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
print(p)
dev.off()

scheck_c <- ldply(lapply(split(scheck,scheck$seurat_clusters.x),function(c){
  xclus=unique(c$seurat_clusters.x)
  cdf <- subset(celltype, cluster==xclus)
  yc <- plyr::count(c,c("seurat_clusters.y"))
  xc <- plyr::count(c,c("seurat_clusters.x"))
  perc=yc$freq/(yc$freq+xc$freq)
  df <- data.frame(cdf,xc,yc,perc)
  return(df)
  }), data.frame)
    p <- ggplot(scheck_c, aes(x = paste0(seurat_clusters.x,"_",celltype), y = seurat_clusters.y,fill=perc)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high =  "steelblue",limits=c(0,1))+
      labs(x="with dex clusters",y="without dex clusters")+
      #scale_x_continuous(breaks = scheck_c$.id)+
      scale_y_continuous(breaks = scheck_c$seurat_clusters.y)+
      theme_minimal()
  png(width = 9, height = 8, file=paste0(figuredir,"dex_vs_nodex_clusters_heat_sumy.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
print(p)
dev.off()


scheck_countm <- melt(scheck_count,id.vars=c("dexcount.seurat_clusters.x","nodexcount.seurat_clusters.y"))
top10m <- reshape2::dcast(seurat_clusters.x  ~ seurat_clusters.y, data=scheck,fun.aggregate=sum,na.rm=T)
#rownames(top10m) <- top10m$cryptic_exon  
myMat <- top10m[,-1]

  png(width = 9, height = 8, file=paste0(figuredir,"dex_vs_nodex_clusters_heat.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
p <- ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw()) +
scale_fill_gradient2(low = "white", high =  "darkblue", mid = "blue", midpoint = 0.5)+
labs(fill = "sum BC")
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
    print(p)
dev.off()

######################################################
### visulization-option 2 for cell type annotation from Julong
#######################################################

library(Matrix)
library(data.table)
library(plyr);library(dplyr)
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(scales)
library(openxlsx)
library(tidyverse)
library(Seurat)

#args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",0.2,"ALL","fastdemux",11) #for testing
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/",0.2,"ALL","demux",13) 

base <- args[1]
resset <- as.numeric(args[2])
project <- args[3]
method <- args[4]
dimset=args[5]
cat("resolution=",resset,"\nproject=",project,"\n","method=",method,"\n")
outdir=paste0(base,"5b_IdenCelltype_",method,"/")
figuredir=paste0(outdir,"figures/")
harmonyfolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")

###
scheck <- fread(file=paste0(outdir,"scheck.txt"))
celltype <- data.frame(cluster.withdex=c("0","1","2","3","4","5","6"),celltype.withdex=c("CD4_T","CD4_T","CD8_T","NK","Monocyte","B","DC"),
  cluster.withoutdex=c("0","1","2","3","4","5","6"),celltype.withoutdex=c("CD4_T","CD8_T","NK","Monocyte1","B","DC","Monocyte2"))

###
### plot data From Julong's script: https://github.com/piquelab/sc-atac-cziHOLD/blob/main/analyses_2024-09-13/1.2_ArchR_process/2_summary_annot.R lines 295-411

plotDF  <- scheck%>%
   group_by(seurat_clusters.x,seurat_clusters.y)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters.y)%>%
   mutate(nt=sum(Freq), Perc=Freq/sum(Freq))%>%
   ungroup()

mat <- plotDF%>%
    pivot_wider(id_cols=seurat_clusters.y, names_from=seurat_clusters.x, values_from=Perc, values_fill=NA)%>%
    column_to_rownames(var="seurat_clusters.y")
mat <- as.matrix(mat)
mat2 <- t(mat)

## sort
rnSel <- data.frame(rn=rownames(mat2))%>%mutate(rn_val=as.numeric(rn))%>%arrange(rn_val)%>%pull(rn)
colSel <- data.frame(cn=colnames(mat2))%>%mutate(cn_val=as.numeric(cn))%>%arrange(cn_val)%>%pull(cn)

mat2 <- mat2[rnSel,colSel]    
colnames(mat2) <- paste0(colnames(mat2),"_",celltype$celltype.withoutdex)
rownames(mat2) <- paste0(rownames(mat2),"_",celltype$celltype.withdex)

### set color
mycol <- colorRamp2(seq(0, 1, length.out=20), colorRampPalette(brewer.pal(n=7, name="YlGnBu"))(20))

### figures
p <- Heatmap(mat2, col=mycol, cluster_rows=F, cluster_columns=F, na_col="white", 
  rect_gp=gpar(col="grey", lwd=0.7),          
  show_row_names=T, row_names_gp=gpar(fontsize=10), 
  show_column_names=T, column_names_gp=gpar(fontsize=9),
  #column_split=split,
  column_title="without dex clusters", column_title_gp=gpar(fontsize=10),
  row_title="with dex clusters", row_title_gp=gpar(fontsize=10), 
  #top_annotation=ha, 
  heatmap_legend_param=list(title="Percent", title_gp=gpar(fontsize=9),
      at=seq(0, 1, by=0.25), labels_gp=gpar(fontsize=8),
      grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(5, "cm")),
  use_raster=T, raster_device="png")

### save figures
  png(width = 9, height = 8, file=paste0(figuredir,"dex_vs_nodex_clusters_heat_norm.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
print(p)
dev.off()

#########################################################################################
#version with no filtering
outdir=paste0(base,"5b_IdenCelltype_",method,"/")
resset=0.2
dimset=11
sc_wdex <- read_rds(paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds"))
sc_wdex@meta.data$letter_clusters <- paste0("C",sc_wdex@meta.data$seurat_clusters)
wdex <- data.frame(barcode=sc_wdex@meta.data$NEW_BARCODE, seurat_clusters=sc_wdex@meta.data$seurat_clusters)
rm(sc_wdex)
gc(reset=T)

nodexoutdir=paste0(base,"5b_IdenCelltype_",method,"/nodex/")
resset=0.15
dimset=13
sc <- read_rds(paste0(nodexoutdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds"))
sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
nodex <- data.frame(barcode=sc@meta.data$NEW_BARCODE, seurat_clusters=sc@meta.data$seurat_clusters)
rm(sc)
gc(reset=T)

scheck <- merge(wdex,nodex,by="barcode")
fwrite(scheck,file=paste0(nodexoutdir,"scheck.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

scheck <- fread(file=paste0(nodexoutdir,"scheck.txt"))
figuredir=paste0(nodexoutdir,"figures/")
#run previous section
  png(width = 9, height = 8, file=paste0(figuredir,"nofilt_dex_vs_nodex_clusters_heat_norm.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
print(p)
dev.off()

###################################################
################################################
#ALOFT ###########################################
library(Matrix)
library(data.table)
library(plyr);library(dplyr)
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(scales)
library(openxlsx)
library(tidyverse)
library(Seurat)

##version with no filtering
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/",0.3,"ALL","demux",50) 
base <- args[1]
project <- args[3]
method <- args[4]

outdir=paste0(base,"5b_IdenCelltype_",method,"/")
resset=0.2
dimset=50
sc_wdex <- read_rds(paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds"))
sc_wdex@meta.data$letter_clusters <- paste0("C",sc_wdex@meta.data$seurat_clusters)
wdex <- data.frame(barcode=sc_wdex@meta.data$NEW_BARCODE, seurat_clusters=sc_wdex@meta.data$seurat_clusters)
rm(sc_wdex)
gc(reset=T)

nodexoutdir=paste0(base,"5b_IdenCelltype_",method,"/nodex/")
resset=0.3
dimset=50
sc <- read_rds(paste0(nodexoutdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds"))
sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
nodex <- data.frame(barcode=sc@meta.data$NEW_BARCODE, seurat_clusters=sc@meta.data$seurat_clusters)
rm(sc)
gc(reset=T)

scheck <- merge(wdex,nodex,by="barcode")
fwrite(scheck,file=paste0(nodexoutdir,"scheck.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

scheck <- fread(file=paste0(nodexoutdir,"scheck.txt"))
figuredir=paste0(nodexoutdir,"figures/")

plotDF  <- scheck%>%
   group_by(seurat_clusters.x,seurat_clusters.y)%>%
   summarize(Freq=n(),.groups="drop")%>%
   group_by(seurat_clusters.y)%>%
   mutate(nt=sum(Freq), Perc=Freq/sum(Freq))%>%
   ungroup()

mat <- plotDF%>%
    pivot_wider(id_cols=seurat_clusters.y, names_from=seurat_clusters.x, values_from=Perc, values_fill=NA)%>%
    column_to_rownames(var="seurat_clusters.y")
mat <- as.matrix(mat)
mat2 <- t(mat)

## sort
rnSel <- data.frame(rn=rownames(mat2))%>%mutate(rn_val=as.numeric(rn))%>%arrange(rn_val)%>%pull(rn)
colSel <- data.frame(cn=colnames(mat2))%>%mutate(cn_val=as.numeric(cn))%>%arrange(cn_val)%>%pull(cn)

mat2 <- mat2[rnSel,colSel]    

celltype <- data.frame(cluster.withdex=c("0","1","2","3","4","5","6","7","8","9","10","11"),celltype.withdex=c("CD4_T","CD4_T","CD8_T","NK","Monocyte","B","DC"),
  cluster.withoutdex=c("0","1","2","3","4","5","6","7","8","9","10","11"),celltype.withoutdex=c("CD4_T","CD8_T","NK","Monocyte1","B","DC","Monocyte2"))

colnames(mat2) <- paste0(colnames(mat2),"_",celltype$celltype.withoutdex)
rownames(mat2) <- paste0(rownames(mat2),"_",celltype$celltype.withdex)

### set color
mycol <- colorRamp2(seq(0, 1, length.out=20), colorRampPalette(brewer.pal(n=7, name="YlGnBu"))(20))

### figures
p <- Heatmap(mat2, col=mycol, cluster_rows=F, cluster_columns=F, na_col="white", 
  rect_gp=gpar(col="grey", lwd=0.7),          
  show_row_names=T, row_names_gp=gpar(fontsize=10), 
  show_column_names=T, column_names_gp=gpar(fontsize=9),
  #column_split=split,
  column_title="without dex clusters", column_title_gp=gpar(fontsize=10),
  row_title="with dex clusters", row_title_gp=gpar(fontsize=10), 
  #top_annotation=ha, 
  heatmap_legend_param=list(title="Percent", title_gp=gpar(fontsize=9),
      at=seq(0, 1, by=0.25), labels_gp=gpar(fontsize=8),
      grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(5, "cm")),
  use_raster=T, raster_device="png")

### save figures
  png(width = 9, height = 8, file=paste0(figuredir,"nofilt_dex_vs_nodex_clusters_heat_norm.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
print(p)
dev.off()
