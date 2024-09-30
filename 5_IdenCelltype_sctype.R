library(tidyverse)
library(data.table)
library(future)
library(HGNChelper)
library(Seurat)
library(plyr)

#######alt cell typing 
args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",0.3,"ALL","fastdemux") #for testing
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/",0.2,"ALL","demux") 

base <- args[1]
resset <- as.numeric(args[2])
project <- args[3]
method <- args[4]
#dimset=50
cat("resolution=",resset,"\nproject=",project,"\n","method=",method,"\n")

outdir=paste0(base,"5b_IdenCelltype_",method,"/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outdir,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

#preharmony
opfn_i <- file.info(dir(paste0(base,"2.1_mergeCellRangerAnd",method,"/"), full.names=T, pattern=paste0(project,".seuratObj-preharmony-post-clustering-res",resset)))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)


  es.max = sctype_score(scRNAseqData = sc[["RNA"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)  #sc[["RNA"]]@scale.data for seurat V<5
  cL_resutls = do.call("rbind", lapply(unique(sc@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sc@meta.data[sc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  sc@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sc@meta.data$customclassif[sc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  png(width = 8, height = 8, file=paste0(figuredir,project,".preharmony_umap_QC_wSCtypecelltype.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
  p <- DimPlot(sc, reduction = "umap", label = TRUE, repel=TRUE, group.by = 'customclassif', pt.size = .1)
  print(p)
  dev.off()
opfn <- paste0(outdir,project,".seuratObj-.preharmony-sctype-",Sys.Date(),".rds")
write_rds(sc, opfn)

#harmony

#opfn_i <- file.info(dir(paste0(base,"2.1_mergeCellRangerAnd",method,"/"), full.names=T, pattern=paste0(project,".seuratObj-post-clustering-res",resset)))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
for (resset in c(0.1, 0.15, 0.2, 0.3, 0.4)){
  cat("running ",resset)
opfn <- paste0(base,"2.1_mergeCellRangerAnd",method,"/",project,".seuratObj-post-clustering-res",resset,".",dimset,".rds")
sc <- read_rds(opfn)

  es.max = sctype_score(scRNAseqData = sc[["RNA"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)  #sc[["RNA"]]@scale.data for seurat V<5
  cL_resutls = do.call("rbind", lapply(unique(sc@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sc@meta.data[sc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sc@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  sc@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sc@meta.data$customclassif[sc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  png(width = 9, height = 8, file=paste0(figuredir,project,".",resset,".",dimset,".harmony_umap_QC_wSCtypecelltype.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
  p <- DimPlot(sc, reduction = "umap", label = TRUE, repel=TRUE, group.by = 'customclassif', pt.size = .1)
  print(p)
  dev.off()
opfn <- paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
write_rds(sc, opfn)

cL_resutls_perc <- ldply(lapply(split(cL_resutls,cL_resutls$cluster), function(i){
  #i <- split(cL_resutls,cL_resutls$cluster)[[1]]
  sumscore=sum(i$scores,na.rm=T)
  i <- transform(i, perc_score=scores/sumscore)
  i <- transform(i, perc_score=ifelse(perc_score>=0,perc_score,0))
  ii <- merge(i,unique(sc@meta.data[,c("seurat_clusters","customclassif")]), by.x="cluster",by.y="seurat_clusters")
  ii <- transform(ii, matching_type=ifelse(type==customclassif,TRUE,FALSE))
  return(ii)
  }),data.frame)[,-1]
opfn <- paste0(outdir,project,".perc_scores.harmony-sctype-",resset,".",dimset,".rds")
write_rds(cL_resutls_perc, opfn)

rm(sc,cL_resutls,es.max)
gc(reset=T)
}

Bcell <- unique(cL_resutls_perc[grep("B cell",cL_resutls_perc$type),])$type
Tcell <- unique(cL_resutls_perc[grep("T cell",cL_resutls_perc$type),])$type
NKcell <- unique(cL_resutls_perc[grepl("Natural killer|NKT",cL_resutls_perc$type),])$type
Monocyte <- unique(cL_resutls_perc[grepl("Monocyte|monocyte|Macrophage",cL_resutls_perc$type),])$type
Dendritic <- unique(cL_resutls_perc[grep("Dendritic",cL_resutls_perc$type),])$type
other <- unique(cL_resutls_perc$type)[!unique(cL_resutls_perc$type) %in% c(Bcell,Tcell,NKcell,Monocyte,Dendritic)]

majorcelltype <- data.frame(majorcelltype=c(rep("Bcell",length(Bcell)),rep("Tcell",length(Tcell)),rep("NKcell",length(NKcell)),rep("Monocyte",length(Monocyte)),rep("Dendritic",length(Dendritic)),rep("other",length(other))),celltype=c(Bcell,Tcell,NKcell,Monocyte,Dendritic,other))

l_perc <- ldply(lapply(c(0.1, 0.15, 0.2, 0.3, 0.4),function(resset){
  #resset=0.2
  opfn <- paste0(outdir,project,".perc_scores.harmony-sctype-",resset,".",dimset,".rds")
  cL_resutls_perc <- read_rds(opfn)
  cL_resutls_perc <- transform(cL_resutls_perc, res=as.factor(resset),perc_score=ifelse(perc_score>=1,1,perc_score))
  cL_resutls_perc_cl <- ldply(lapply(split(cL_resutls_perc,cL_resutls_perc$cluster),function(c){
    #c <- split(cL_resutls_perc,cL_resutls_perc$cluster)[[8]]
    majorcelltype_c <- subset(majorcelltype, celltype %in% c$type)
    summajor <- ldply(lapply(split(majorcelltype_c,majorcelltype_c$majorcelltype),function(i){
    #i <- split(majorcelltype_c,majorcelltype_c$majorcelltype)[[1]]
    df <- subset(c, type %in% i$celltype)
    score=data.frame(majorcelltype=unique(i$majorcelltype), res=as.factor(resset), cluster=unique(c$cluster),sumscore=sum(df$perc_score,na.rm=T))
    return(score)
    }),data.frame)[,-1]
    #topcell=summajor[which.max(summajor$sumscore),]$majorcelltype
    totalsumscore=sum(summajor$sumscore,na.rm=T)
    summajor=transform(summajor,correctness=sumscore/totalsumscore)
    summajor <- transform(summajor, assign=ifelse(correctness>0.5, "TRUE", "FALSE"))
    fullassign=ifelse(length(which(summajor$assign==TRUE))==1,"single", ifelse(length(which(summajor$assign==TRUE))>1,"multiple","mixed"))
    summajor <- transform(summajor, assignedcelltype=unique(c$customclassif),fullassign=fullassign)
    summajor <- transform(summajor, assignedcelltype=ifelse(assignedcelltype=="Macrophages","Macrophage (Monocyte)",assignedcelltype))
    return(summajor)
    }), data.frame)[,-1]
  return(cL_resutls_perc_cl)
}),data.frame)
l_perc <- transform(l_perc, letter_clusters=paste0("C",cluster))

test <- unique(l_perc[,c("res","cluster","fullassign")])
table(test$res,test$fullassign)

top10m <- reshape2::dcast(res ~ letter_clusters, data=l_perc, value.var="perc_score",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$res
myMat <- top10m[,-1]
opfn <- paste0(outdir,project,".perc_scores.harmony-sctype-heatmap.png")
png(width = 15000, height = 15000, file=opfn, pointsize=12, 
      bg = "transparent", res = 800)
par(mar=c(5,5,4,2)+0.5) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2,digits = 1) +
scale_fill_gradient2(low = "white", high =  "darkblue", mid = "lightblue", midpoint = 0.5)
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
dev.off()

for (resset in c(0.1, 0.15, 0.2, 0.3, 0.4)){
  #resset=0.2
p <- ggplot(l_perc[l_perc$res==resset,], aes(x=paste0(cluster,"\n",customclassif), y=type, fill=perc_score))+
  geom_tile()+
  xlab("cluster and assigned cell type")+
  ylab("cell type")+
  ggtitle(paste0("res=", resset))+
  scale_fill_gradient("perc_score", low="#ffffc8", high="#7d0025", na.value=NA)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
opfn <- paste0(outdir,project,".",resset,".",dimset,".perc_scores.harmony-sctype-heatmap.png")
png(width = 11, height = 7, file=opfn, pointsize=12, 
      bg = "transparent", units = "in", res = 800)
print(p)
dev.off()
}


unique(sc@meta.data[,c("seurat_clusters","customclassif")])

