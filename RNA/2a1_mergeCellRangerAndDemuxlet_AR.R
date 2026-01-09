library(Seurat)
library(Matrix)
library(future)
library(readr)
library(tidyverse)
library(harmony)


#####################################################################
### 07/25/2023, Ali R                                          ###### 
###  SCAIP7-18 CellRanger alignment and demultiplexing         ######
###  comparison of CellRanger and Kallisto                     ######
###  include post renaming EtOH -> CTRL barcod reruns          ######
###  modified from pre-eclampsia paper:                        ######
###                    2_mergeCellRangerAndDemuxlet.R          ######
#####################################################################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","fastdemux") #for testing "CZI2_group.txt"
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/", "demux")
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","fastdemux")

base <- args[1]
method <- args[2]
#read in samples file (just list of samples to run, each sample on newline)
if(!is.na(args[3])){
samples=read.table(paste0(base,args[3]),header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
project=sapply(strsplit(args[3],"_"),function(y)y[1])
cat("samplefile= ",args[3])
}else{
   project="ALL"
}
if(method=="demux"){
   outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
   demux_in <- paste0(base,"1_demux_output/")
   kallisto_in <- paste0("2b_mergeKallistoAnd",method,"/")
   opfn <- paste0(demux_in,project,".1_demux_filt.SNG.rds")
   demux <- read_rds(opfn)
} else {
   outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
   demux_in <- paste0(base,"1_demux_alt_output/")
   opfn <- paste0(demux_in,project,".1_demux_alt_filt.SNG.rds")
   demux <- read_rds(opfn)   
   kallisto_in <- paste0("2b_mergeKallistoAnd",method,"/")
}

if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_2024-04-19/",base)  ##

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

libList <- scan(paste0(basefolder,"libList.txt"),what=character(0)) #was libList.txt, not sure why that only has 4 samples
if(!is.na(args[3])){
  libList <- libList[libList %in% samples$V1]
}


##Removig DEX here:
libList <- libList[grep("DEX",libList,invert=TRUE)]

cat("creating seurat object with ",libList)

sc_list<-sapply(libList, function(x){
##  x<-libList[1]
  cat("#Processing: ",x,"\n")
  gp.data<- Read10X(data.dir = paste0(basefolder,x,"/filtered_feature_bc_matrix"))
  #################################################################################
  # creating seurat object
  #################################################################################
  sc <- CreateSeuratObject(counts = gp.data, project = paste0("cellranger-CZI.",project),min.cells = 3, min.features=200)
  sc@meta.data$Library<-rep(x,nrow(sc@meta.data))
  sc
} )

## find matching barcodes demuxlet and the sc object. 

opfn <- paste0(outFolder,project,".seuratObj-merge.",Sys.Date(),".rds") 
write_rds(sc_list, opfn)

# read it again if needed
#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-merge."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc_list <- read_rds(opfn)

sc <- merge(sc_list[[1]],sc_list[-1],add.cell.ids = libList, project=paste0("cellranger-CZI.",project))

rm(sc_list)


sc[["RNA"]] <- JoinLayers(sc[["RNA"]])

opfn <- paste0(outFolder,project,".seuratObj-all-ulnist-prior-to-demux.",Sys.Date(),".rds") 
write_rds(sc, opfn)

