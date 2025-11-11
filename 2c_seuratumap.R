library(Seurat)
library(future)
library(tidyverse)
library(harmony)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",13,"ALL","fastdemux","CZI") #for testing
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/",12,"ALL","demux","ALOFT")

base <- args[1]
project <- args[3]
method <- args[4]
job <- args[5]
outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_2024-04-19/",base)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

dim <- as.numeric(args[2])
filter <- "noDEX"
#filter <- "CTRLonly" #ALOFT used

future::plan(strategy = 'multicore', workers = 7)
options(future.globals.maxSize = 100 * 1024 ^ 3)

########################

#test both dim picked from elbowplot and default
for (dimset in c(dim,50)){ 
	cat("running",dimset,"\n")
#opfn_i <- file.info(dir(outFolder, full.names=T, pattern=paste0(project,".seuratObj-afterharmony.")))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
opfn <- paste0(outFolder,project,".seuratObj-afterharmony.",filter,".rds")
sc <- read_rds(opfn)

sc <- sc %>% RunUMAP(reduction = "harmony", dims = 1:dimset) 
sc <- sc %>% FindNeighbors(reduction = "harmony", dims = 1:dimset)

opfn <- paste0(outFolder,project,".seuratObj-post-umap.",dimset,".",filter,".rds")
write_rds(sc, opfn)
}

