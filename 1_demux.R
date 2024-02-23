######################################
### demux2 for downstream analysis ###
######################################
library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/CZI/RNA/counts_cellranger_hg38/demuxlet/demuxlet/","CZ1_group.txt") #for testing
analysis=args[1]
basefolder=args[2]
if(!is.na(args[3])){
samples=fread(args[3],header=F)$V1
}
outdir <- paste0(analysis,"1_demux_output/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

#################################
### 1. folders of counts data ###
#################################
demuxfn <- list.files(basefolder,pattern="*.out.best")

expNames_in <- gsub(".out.best", "", demuxfn) 

if(!is.na(args[3])){
expNames <- expNames_in[expNames_in %in% samples] #run subset of all output in the folder
} else {
expNames <- expNames_in
}

############################
### 2, read demuxlet data ###
############################

# new column "NEW_BARCODE" was added to the data frame associated to each experiment (experiment name+barcode) "-1" was removed from end of barcode
demux <- mclapply(expNames,function(ii){
  cat("#Loading ", ii, "\n")
  fn <- paste0(basefolder, ii,".out.best")
  dd <- data.frame(fread(fn,header=T))
  dd <- dd%>%mutate(NEW_BARCODE=paste0(ii,"_", gsub("-1","",BARCODE)),EXP=ii)
},mc.cores=10)

# merging all the experiments 
demux <- do.call(rbind,demux)

### output
opfn <- paste0(outdir,"1_demux_New.ALL.rds");
write_rds(demux, opfn)

# DROPLET.TYPE: "SNG" "AMB"
opfn <- paste0(outdir,"1_demux_New.SNG.rds");
demux %>% filter(DROPLET.TYPE=="SNG") %>%
    write_rds(opfn)


