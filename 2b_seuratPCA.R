library(Seurat)
library(future)
library(tidyverse)
library(harmony)

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

##############################

### load the seurat object
opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-postmerge-after-mt-filtering."))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

#sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5) # could try mt<10 to increase data

sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)
sc <- ScaleData(sc, features = rownames(sc)) #can regress out highly variable input with flag: vars.to.regress = "percent.mt"
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
#sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

png(width=1000, height=1000, res=120, file=paste0(figuredir,"pca1-15_heatmap_QC.png"), bg = "transparent")
p <- DimHeatmap(sc, dims = 1:15, cells = 500, balanced = TRUE) # for multiple PCs
print(p)
dev.off()
png(width = 1000, height = 1000, file=paste0(figuredir,"elbowplot_QC.png"), bg = "transparent", res = 120)
p <- ElbowPlot(sc)
print(p)
dev.off()

opfn <- paste0(outFolder,"seuratObj-afterPCA.",Sys.Date(),".rds") 
write_rds(sc, opfn)
