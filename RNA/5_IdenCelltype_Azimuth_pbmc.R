##############################################
### cell type classification using Azimuth ###
##############################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(pheatmap)
library(ggrepel)

library(Azimuth)
library(sctransform)

#################

args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",0.1,"ALL","fastdemux",13)

future::plan(strategy = 'multicore', workers = 7)
options(future.globals.maxSize = 100 * 1024 ^ 3)

base <- args[1]
resset <- as.numeric(args[2])
project <- args[3]
method <- args[4]
dimset=args[5]
cat("resolution=",resset,"\nproject=",project,"\n","method=",method,"\n")
filter <- "noDEX"
##filter <- "CTRLonly" #ALOFT used

indir=paste0(base,"5b_IdenCelltype_",method,"/",filter,"/")


outdir=paste0(base,"5b_IdenCelltype_",method,"Azimuth_PBMC/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)
outdir=paste0(base,"5b_IdenCelltype_",method,"Azimuth_PBMC/",filter,"/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outdir,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#anno <- read_rds("/wsu/home/groups/prbgenomics/preeclampsia/preeclampsia_sc/3_runHarmony/anno-after-mt-filtering.2023-05-07.rds")

#cell.type.anno <- read_tsv("cell.type.anno.res0p5_2023-05-07.tsv")
#cell.type.anno$seurat_clusters <- factor(cell.type.anno$seurat_clusters)

opfn <- paste0(indir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
sc2 <- readRDS(opfn)

# Set default to RNA and merge all layers into one (fixes the warning)
DefaultAssay(sc2) <- "RNA"
sc2 <- JoinLayers(sc2)  # Merges data.* layers into single 'data'

#sc2 <- SCTransform(sc2, verbose = TRUE)

## Run Azimuth may need the following line to run. 
#httr::set_config(httr::config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

sc2 <- RunAzimuth(sc2, reference = "pbmcref")

md <- sc2@meta.data

head(md)

tt <-  table(md$predicted.celltype.l1,md$seurat_clusters)

tt


write.csv(md,paste0(outdir,"sc2_metadata_Azimuth_PBMC.csv"))

write.csv(tt,paste0(outdir,"Azimuth_PBMC_celltype_counts_per_clusters.csv"))

saveRDS(sc2,paste0(outdir,"Azimuth_PBMC.Seurat.rds"))


png(width = 8, height = 8, file=paste0(figuredir,"heatmap_celltype_Azimuth_PBMC.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
p <- pheatmap(tt,cluster_rows=TRUE,cluster_cols=FALSE,scale = "column")
print(p)
dev.off()

##########################################################################
