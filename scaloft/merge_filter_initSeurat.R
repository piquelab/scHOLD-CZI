#merge and filter ALOFT seurat objects
library(Seurat)
library(future)
library(tidyverse)
library(plyr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt") #for testing "CZI2_group.txt"
base <- args[1]
method <- args[2]
sample_batch <- args[3]

#read in samples file (just list of samples to run, each sample on newline)
project="ALL"

outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_hg38/",base)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

scALOFT1="/rs/rs_grp_scaloft/scALOFT_2024/2.1_mergeCellRangerAndDemuxlet/scALOFT1/seuratObj-after-mt-filtering.2024-04-10.rds"
scALOFT2="/rs/rs_grp_scaloft/scALOFT_2024/2.1_mergeCellRangerAndDemuxlet/scALOFT2/seuratObj-after-mt-filtering.2024-04-10.rds"
sc_1 <- read_rds(scALOFT1)
sc_2 <- read_rds(scALOFT2)
sc_list <- list(scALOFT1=sc_1,scALOFT2=sc_2)

sc <- merge(sc_list[[1]],sc_list[-1])
#sc[["RNA"]] <- JoinLayers(sc[["RNA"]]) #dont seem to need to do this here?

opfn <- paste0(outFolder,project,".seuratObj-sc1and2_prefilter.rds") 
write_rds(sc, opfn)

# filter adapted from 1b_demux_anal3-AR.R
#this 10 filter was used for dumuxlet input. unecessary for alternative input as already filtered for >100
aa <- sc@meta.data %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID) 

cell.counts <- sc@meta.data %>% group_by(EXP,Sample_ID, BATCH, treats) %>% summarize(n=n()) 
cell.counts <- cell.counts%>%group_by(Sample_ID)%>%
   mutate(Perc=n/sum(n), comb=paste(EXP, Sample_ID, sep="_"))%>%as.data.frame()

cell.counts.filt <- cell.counts %>% filter(n>100) 
n100toremove <- cell.counts %>% filter(!n>100) 

cat("sum cell counts filtered n>100= ",sum(cell.counts.filt$n))
cat("unique cell counts filtered n>100 samples= ",length(unique(cell.counts.filt$Sample_ID)))

cell.counts.filt <- cell.counts.filt %>% mutate(ID_BATCH=paste(Sample_ID,BATCH, sep="_"))

write_tsv(cell.counts.filt,paste0(outFolder,project,".cell.count.filter.100min.v1.tsv"))

IDmiss <- n100toremove %>% select(comb,n) 
write.csv(IDmiss, paste0(outFolder,project,".missing_samples_in_demultiplexing_results_duetoFilter.csv"), row.names=F)

##################################################################
############ experimental data sheet to find unmatched ###########
##################################################################

cat("using experimental data sheet to find unmatched")

# load experimental data sheet
#covariate file numbers are in format 01,02,etc which does not match file naming labels of 1,2,etc
exp <- read.table(sample_batch, row.names=NULL,header=T)
exp <- exp %>% mutate(Sample_ID=dbgap.ID)

# merge experimental data and filtered data
#expsub <- exp %>% mutate(Sample_ID=dbgap.ID) %>% select(dbgap.ID, Batch, Sample_ID) #removed Participant.ID as not a variable
d <- subset(plyr::count(exp,"Sample_ID"),freq>1)
dp <- subset(plyr::count(exp,"Sample_ID"),freq<2)


		merge <- merge(cell.counts.filt, d, by.x = c("Sample_ID","BATCH"),by.y=c("Sample_ID","Batch2"))
		A_sum=sum(merge[grep("A$",merge$BATCH),"n"],na.rm=T)
		B_sum=sum(merge[grep("B$",merge$BATCH),"n"],na.rm=T)
		if(A_sum>B_sum){d_sum <- (unique(merge[grep("A$",merge$BATCH),]))} else if(B_sum>A_sum){d_sum <- (unique(merge[grep("B$",merge$BATCH),]))} else if(isTRUE(unique(merge$Batch)=="SCAIP5")){d_sum <- (unique(merge[grep("ChemV3",merge$chemistry),]))} 



merge <- rbind(dl, transform(merge(cell.counts.filt[,c("Sample_ID","EXP","BATCH","treats","n","Perc","comb","ID_BATCH")], subset(exp[,c(1:3,5)], Sample_ID %in% dp$Sample_ID), by = "Sample_ID",all.x=T),ID_Batch=paste0(dbgap.ID, "_", Batch2)))

#SCAIP5V3 - SCAIP18 is Chemistry V3
V3chem <- c("SCAIP5V3","SCAIP6","SCAIP18","SCAIP13","SCAIP9","SCAIP8","SCAIP15","SCAIP10","SCAIP11","SCAIP12","SCAIP6A","SCAIP6B","SCAIP14","SCAIP17","SCAIP16")
#merge <- subset(merge,BATCH %in% V3chem) #first run removed all V2, now want to only remove SCAIP5V2
merge <- subset(merge,!BATCH=="SCAIP5V2") #first run removed all V2, now want to only remove SCAIP5V2

###
#sample batch file with subset samples:
mergeexp <- rbind(dl[,c("Sample_ID","Batch")], subset(exp[,c("Sample_ID","Batch")], Sample_ID %in% dp$Sample_ID))
mergeexp <- subset(mergeexp,Batch %in% V3chem)
fwrite(mergeexp,sep='\t',quote=F,row.names=F,file=paste0(base,"filtered.",gsub(".*/","",sample_batch)))
#first run removed all V2, now want to only remove SCAIP5V2
mergeexp <- rbind(dl[,c("Sample_ID","Batch")], subset(exp[,c("Sample_ID","Batch")], Sample_ID %in% dp$Sample_ID))
mergeexp <- subset(mergeexp,!Batch=="SCAIP5V2")
fwrite(mergeexp,sep='\t',quote=F,row.names=F,file=paste0(base,"filtered.incV2.",gsub(".*/","",sample_batch)))

###


#there were a bunch of summary items here, kept only QC
cat("checking NA ids and mismatches of ids")
sum(is.na(merge$Sample_ID))
sum(is.na(merge$dbgap.ID))
sum(merge$Sample_ID != merge$dbgap.ID,na.rm=T)
sum(merge$Batch != merge$BATCH,na.rm=T)

# samples where IDs dont match
unmatch_ID <- merge %>%
  dplyr::filter(Sample_ID != dbgap.ID)
print(unmatch_ID)

# samples where ID_batch dont m

# samples where demultiplexing BATCH does not match to the experimental batch for a sampleIDs
#unmatch_batch <- merge %>% dplyr::filter(BATCH != Batch2) # -- at least in current data seems to only be instances of NA in exp batch
unmatch_comb <- subset(merge, ID_BATCH != ID_Batch)
unmatch_comb <- subset(unmatch_comb, !ID_Batch=="NA_NA")
print(unmatch_comb)

write.csv(unmatch_comb, paste0(outFolder,project,".umatched_batches_samples_IDs_compared_to_exp_file.csv"), row.names=F)

### remove the unmatched from the >100 filtered data. 
dim(cell.counts.filt)
unmatch_comb$comb

# how many cells are unmatched
cell.unm <- cell.counts.filt %>% dplyr::filter(comb %in% unmatch_comb$comb)
dim(cell.unm)
sum(cell.unm$n)

# how many cells remain
rem.unm <- cell.counts.filt %>% dplyr::filter(!comb %in% unmatch_comb$comb)
dim(rem.unm)
sum(rem.unm$n)

write_tsv(rem.unm,paste0(outFolder,project,".cell.count.filter.100min.v2_unmatched_cells_removed.tsv"))

rem.unm <- fread(paste0(outFolder,project,".cell.count.filter.100min.v2_unmatched_cells_removed.tsv"))

opfn <- paste0(outFolder,project,".seuratObj-sc1and2_prefilter.rds") 
sc <- read_rds(opfn)
sc@meta.data$comb <- paste0(sc@meta.data$EXP,"_",sc@meta.data$Sample_ID)
sc <- subset(sc, subset=comb %in% rem.unm$comb)

opfn <- paste0(outFolder,project,".seuratObj-postmerge-after-mt-filtering.",Sys.Date(),".rds") 
write_rds(sc, opfn)

#I CANT BELIEVE GITHUB DIDNT SAVE MY NEW SCRIPT
############ OLD METHOD #########################
################################################
dl <- ldply(lapply(unique(d$Sample_ID),function(i){
		df <- subset(exp, Sample_ID==i)
		merge <- merge(cell.counts.filt, df, by.x = c("Sample_ID","BATCH"),by.y=c("Sample_ID","Batch2"))
		A_sum=sum(merge[grep("A$",merge$BATCH),"n"],na.rm=T)
		B_sum=sum(merge[grep("B$",merge$BATCH),"n"],na.rm=T)
		if(A_sum>B_sum){d_sum <- (unique(merge[grep("A$",merge$BATCH),]))} else if(B_sum>A_sum){d_sum <- (unique(merge[grep("B$",merge$BATCH),]))} else if(isTRUE(unique(merge$Batch)=="SCAIP5")){d_sum <- (unique(merge[grep("ChemV3",merge$chemistry),]))} 
		return(d_sum)
		}), data.frame)
dl <- transform(dl, Batch2=BATCH, BATCH=gsub("[AB]$","",dl$BATCH),ID_BATCH=gsub("[AB]$","",dl$ID_BATCH),ID_Batch=paste0(dbgap.ID, "_", Batch))
dl <- transform(dl, BATCH=gsub("V3","",dl$BATCH),ID_BATCH=gsub("V3","",dl$ID_BATCH))
dl <- dl[,c("Sample_ID","EXP","BATCH","treats","n","Perc","comb","ID_BATCH","Batch","dbgap.ID","Batch2","ID_Batch")]


##############################
#table of number of cells per library
cell.counts.filt <- fread(paste0(outFolder,project,".cell.count.filter.100min.v1.tsv"))
n100filtered <- fread(paste0(outFolder,project,".missing_samples_in_demultiplexing_results_duetoFilter.csv"))
n100filtered <- transform(n100filtered, EXP=sapply(strsplit(comb,"_"),function(y) y[1]))
allcellcounts <- rbind(n100filtered,cell.counts.filt,fill=T)
cellperlib <- ddply(rem.unm, c("EXP"), plyr::summarize,
	ncell=sum(n,na.rm=T))
fwrite(cellperlib,sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,project,".cell.count.nofilt_cellsperlibrary.txt"))


#final after filtering
rem.unm <- fread(paste0(outFolder,project,".cell.count.filter.100min.v2_unmatched_cells_removed.tsv"))

cellperlib <- ddply(rem.unm, c("EXP"), plyr::summarize,
	ncell=sum(n,na.rm=T))

fwrite(cellperlib,sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,project,".cell.count.filter.100min.v2_unmatched_cells_removed_cellsperlibrary.txt"))

############
#original cell counts from cellranger summary
summarycellranger <- fread("/rs/rs_grp_scaloft/scALOFT_2024/counts_cellranger_hg38/summary.tsv")
summarycellranger_cellperlibrary <- unique(summarycellranger[,c(1:2)])
fwrite(summarycellranger_cellperlibrary,sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,project,"summarycellranger_cellsperlibrary.txt"))
