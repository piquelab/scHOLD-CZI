#preharmony umap plotting. Done in case harmony removes effects of interest

library(Seurat)
library(future)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","ALL","fastdemux",11) #for testing
base <- args[1]
project <- args[2]
method <- args[3]

outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")

if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_hg38/",base)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

dimset <- as.numeric(args[4])

#################

opfn_i <- file.info(dir(outFolder, full.names=T, pattern=paste0(project,".seuratObj-afterPCA.")))
opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
sc <- read_rds(opfn)

sc <- sc %>% RunUMAP(dims = 1:dimset) 
sc <- sc %>% FindNeighbors(dims = 1:dimset)

opfn <- paste0(outFolder,project,".seuratObj-preharmony-post-umap.",Sys.Date(),".rds")
write_rds(sc, opfn)
