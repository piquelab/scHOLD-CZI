library(Seurat)
library(future)
library(tidyverse)
library(harmony)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","ALL","fastdemux","CZI") #for testing
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","ALL","demux","ALOFT")
base <- args[1]
project <- args[2]
method <- args[3]
job <- args[4]
outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 7)
options(future.globals.maxSize = 100 * 1024 ^ 3)

##############################

### load the seurat object
opfn_i <- file.info(dir(outFolder, full.names=T, pattern=paste0(project,".seuratObj-postmerge-after-mt-filtering.")))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- readRDS(opfn)

#sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5) # could try mt<10 to increase data
if(job=="CZI"){
sc <- subset(sc,subset=treats != "RNA-LPS-DEX") 
filter="noDEX"
} else if(job=="ALOFT"){
	#sc <- subset(sc,subset=treats != "LPS-DEX") 
	sc <- subset(sc,subset=treats == "CTRL") 
	filter="CTRLonly"
}

sc <- NormalizeData(sc)
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)
sc <- ScaleData(sc, features = rownames(sc)) #can regress out highly variable input with flag: vars.to.regress = "percent.mt"
sc <- RunPCA(sc, features = VariableFeatures(object = sc))
#sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)
opfn <- paste0(outFolder,project,".seuratObj-afterPCA.",filter,".rds") 
write_rds(sc, opfn)

png(width=1000, height=1000, res=120, file=paste0(figuredir,project,".",filter,".pca1-15_heatmap_QC.png"), bg = "transparent")
p <- DimHeatmap(sc, dims = 1:15, cells = 500, balanced = TRUE) # for multiple PCs
print(p)
dev.off()
png(width = 1000, height = 1000, file=paste0(figuredir,project,".",filter,".elbowplot_QC.png"), bg = "transparent", res = 120)
p <- ElbowPlot(sc)
print(p)
dev.off()

gc()

sc <- RunHarmony(sc,"Library")#had to remove reduction="pca" with SeuratV5 #kmeans_init_nstart=100, kmeans_init_iter_max=5000 did not remove quick transfe warning

opfn <- paste0(outFolder,project,".seuratObj-afterharmony.",filter,".rds") 
write_rds(sc, opfn)
