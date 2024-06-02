library(tidyverse)
library(data.table)
library(future)
library(HGNChelper)
library(Seurat)

#######alt cell typing 
args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",0.3,"ALL","fastdemux") #for testing
base <- args[1]
resset <- as.numeric(args[2])
project <- args[3]
method <- args[4]

cat("resolution=",resset,"\nproject=",project,"\n","method=",method,"\n")

outdir=paste0(base,"5b_IdenCelltype_",method,"/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outdir,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

opfn_i <- file.info(dir(paste0(base,"2.1_mergeCellRangerAnd",method,"/"), full.names=T, pattern=paste0(project,".seuratObj-preharmony-post-clustering-res",resset)))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

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

opfn_i <- file.info(dir(paste0(base,"2.1_mergeCellRangerAnd",method,"/"), full.names=T, pattern=paste0(project,".seuratObj-post-clustering-res",resset)))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

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
  png(width = 9, height = 8, file=paste0(figuredir,project,".harmony_umap_QC_wSCtypecelltype.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
  p <- DimPlot(sc, reduction = "umap", label = TRUE, repel=TRUE, group.by = 'customclassif', pt.size = .1)
  print(p)
  dev.off()
opfn <- paste0(outdir,project,".seuratObj-.harmony-sctype-",Sys.Date(),".rds")
write_rds(sc, opfn)
