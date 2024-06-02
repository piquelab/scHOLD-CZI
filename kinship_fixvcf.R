library(tidyverse)
library(parallel)
library(viridis)
library(pheatmap)

### 
args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_HOLD01-HOLD14_dbgap.ID_cziexp_fixed_05_03_2024.txt","alternative","CZI2_group.txt") #for testing
cov_file=args[2]

#read in samples file (just list of samples to run, each sample on newline)
if(!is.na(args[4])){
samples=read.table(args[4],header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
cat("samplefile=",args[4],"\n")
project=sapply(strsplit(args[4],"_"),function(y)y[1])
} else{
   project="ALL"
}

cat("project=",project,"\n")

if(args[3]=="demux"){
   outFolder=paste0(args[1],"1_demux_output/")
   opfn <- paste0(outFolder,project,".1_demux_New.SNG.rds")
} else {
   outFolder=paste0(args[1],"1_demux_alt_output/")
   opfn <- paste0(outFolder,project,".1_demux_alt_New.SNG.rds")
}
demux <- read_rds(opfn)
demux_c <- ddply(demux, c("Sample_ID"),plyr::summarize,
   bcnum_t=sum(bcnum,na.rm=T))

kin <- fread(paste0(args[1],"kinship.txt"))
colnames(kin) <- c("ID1","ID2","NSNP","HETHET","IBS0","KINSHIP")

kin1 <- merge(kin,demux_c, by.x="ID1",by.y="Sample_ID")
colnames(kin1) <- c("ID1","ID2","NSNP","HETHET","IBS0","KINSHIP","bcnum_ID1")
kin2 <- merge(kin,demux_c, by.x="ID2",by.y="Sample_ID")
colnames(kin2) <- c("ID1","ID2","NSNP","HETHET","IBS0","KINSHIP","bcnum_ID2")
kin <- cbind(kin1,kin2[,"bcnum_ID2"])
kin <- transform(kin, indpick=ifelse(bcnum_ID1>bcnum_ID2, ID1, ID2))

toremove <- data.frame(ind=c("HO-047","HO-110","HO-116","HO-117","HO-118"))
fwrite(toremove,quote=F,row.names=F, col.names=F,file=paste0(args[1],"kinshipindtoremove.txt"))

#this is run on instance without removing gnu
module load bcftools/1.9
bcftools view -Oz -S ^/rs/rs_grp_schold/CZI/RNA/analysis/kinshipindtoremove.txt /rs/rs_grp_schold/reimpute/ref.ac1.vcf.gz > /rs/rs_grp_schold/CZI/RNA/analysis/ref.ac1.removekin.vcf.gz