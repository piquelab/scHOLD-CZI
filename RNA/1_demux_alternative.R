######################################
### demux2 starting with alternative program for downstream analysis ###
######################################
library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/CZI/RNA/counts_cellranger_2024-04-19/fastdemux/fdout/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_HOLD01-HOLD14_dbgap.ID_cziexp_fixed_05_03_2024.txt","CZI2_group.txt") #for testing
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/CZI/RNA/counts_cellranger_2024-04-19/fastdemux/fdout/","/rs/rs_grp_schold/covariates/dbgap/HOLD_library_metadata_prelim_n165_dbgapIDs_batch_09_17_2024.txt") ## RPR run 2025-11-10
analysis=args[1]
basefolder=args[2]
if(!is.na(args[4])){
samples=read.table(args[4],header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
project=sapply(strsplit(args[4],"_"),function(y)y[1])
} else {
   project="ALL"
}
cat("project=",project,"\n")

cov_file=args[3]
exp <- read.table(cov_file, row.names=NULL,header=T)
exp$Batch <- gsub("HOLD0","HOLD",exp$Batch)
exp <- transform(exp,ID_Batch= paste0(dbgap.ID, "_", Batch))
outdir <- paste0(analysis,"1_demux_alt_output/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

#################################
### 1. folders of counts data ###
#################################
demuxfn <- list.files(basefolder,pattern="*.info.2nd.txt.gz")
expNames_in <- gsub(".fdout.raw.info.2nd.txt.gz", "", demuxfn) 

if(!is.na(args[4])){
expNames <- expNames_in[expNames_in %in% samples$V1] #run subset of all output in the folder
exp <- exp %>% dplyr::filter(Batch %in% samples$Batch)
} else {
expNames <- expNames_in
}

############################
### 2, read demuxlet data ###
############################

# new column "NEW_BARCODE" was added to the data frame associated to each experiment (experiment name+barcode) "-1" was removed from end of barcode
demux <- mclapply(expNames,function(ii){
  #ii <- expNames[1]
  cat("#Loading ", ii, "\n")
  fn <- paste0(basefolder, ii,".fdout.raw.info.2nd.txt.gz")
  dd <- data.frame(fread(fn,header=T))
  dd <- dd%>%mutate(NEW_BARCODE=paste0(ii,"_", gsub("-1","",BARCODE)),EXP=ii,Sample_ID=BestSample,NUM.READS=Numi,NUM.SNPS=Nsnp)
  missing <- setdiff(exp$dbgap.ID, dd$Sample_ID)
  IDmiss <- exp %>% dplyr::filter(dbgap.ID %in% missing) %>% mutate(Sample_ID=dbgap.ID) %>% select(dbgap.ID, Batch, Sample_ID) 
  write.csv(IDmiss, paste0(outdir,"missing_samples_in_demultiplexing_results.csv"), row.names=F)
  dd <- transform(dd, BATCH=gsub("-.*", "", EXP), treats=gsub(".*[0-9].{,2}-","",EXP))
  dd <- dd[dd$Sample_ID %in% exp$dbgap.ID,] #run subset of all output in the folder
  dd <- dd %>% select(NEW_BARCODE,NUM.READS,NUM.SNPS,DropType,bcnum,EXP,BATCH,treats,Sample_ID)
  return(dd)
},mc.cores=10)

# merging all the experiments 
demux <- do.call(rbind,demux)

### output
opfn <- paste0(outdir,project,".1_demux_alt_New.ALL.rds")
write_rds(demux, opfn)

# DROPLET.TYPE: "SNG" "AMB"
opfn <- paste0(outdir,project,".1_demux_alt_New.SNG.rds")
demux <- demux %>% dplyr::filter(DropType==1,NUM.READS>100,NUM.SNPS>100)
write_rds(demux,opfn)
