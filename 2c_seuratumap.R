library(Seurat)
library(future)
library(tidyverse)
library(harmony)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/",11,"ALL","fastdemux") #for testing
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/",13,"ALL","demux")

base <- args[1]
project <- args[3]
method <- args[4]
outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_hg38/",base)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

dimset <- as.numeric(args[2])

future::plan(strategy = 'multicore', workers = 4)
options(future.globals.maxSize = 30 * 1024 ^ 3)

########################

opfn_i <- file.info(dir(outFolder, full.names=T, pattern=paste0(project,".seuratObj-afterharmony.")))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

sc <- sc %>% RunUMAP(reduction = "harmony", dims = 1:dimset) 
sc <- sc %>% FindNeighbors(reduction = "harmony", dims = 1:dimset)

opfn <- paste0(outFolder,project,".seuratObj-post-umap.",dimset,".rds")
write_rds(sc, opfn)
