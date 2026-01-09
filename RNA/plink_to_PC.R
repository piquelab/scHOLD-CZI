library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_HOLD01-HOLD14_dbgap.ID_cziexp_fixed_05_03_2024.txt","CZ1_group.txt") #for testing
base <- args[1]
cov_file=args[2]

#genofolder=paste0(base,"genotypePC/")
genofolder=paste0(base,"genotypePCnokin/")

# set new output dir for filtered out unmatched figures
figuredir=paste0(genofolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#read in samples file (just list of samples to run, each sample on newline)
if(!is.na(args[3])){
samples=read.table(args[3],header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
project=sapply(strsplit(args[3],"_"),function(y)y[1])
cat("samplefile= ",args[3])
}else{
   project="ALL"
}

exp <- read.table(cov_file, row.names=NULL,header=T)
exp <- transform(exp, Sample_ID=dbgap.ID)
exp$Batch <- gsub("HOLD0","HOLD",exp$Batch)
if(!is.na(args[3])){
  exp <- exp %>% dplyr::filter(Batch %in% samples$Batch)
}

eigen <- fread(paste0(genofolder,"ref.eigenval"))$V1
eigen/sum(eigen)

eigenvec <- fread(paste0(genofolder,"ref.eigenvec"), header=T)
names(eigenvec)[2] <- "Sample_ID"
eigenvec <- eigenvec %>%
  dplyr::select(PC1, PC2, PC3, Sample_ID)

eigenvec2 <- eigenvec %>% inner_join(exp[,-1], by="Sample_ID") #removing column 1 aka dbgap id dup of sample id col

fwrite(eigenvec2, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(genofolder,project,".eigenvec_pc.txt"))

p1 <- ggplot(eigenvec2, aes(x=PC1, y=PC2, color=factor(Batch)))+
   geom_point(shape=0.8, size=1)+
   theme_bw()+
   theme(legend.title=element_blank())

figfn <- paste0(figuredir,project,".1.1_pca.png")
png(figfn, width=520, height=420, res=120)
print(p1)
dev.off()

p1 <- ggplot(eigenvec2, aes(x=PC2, y=PC3, color=factor(Batch)))+
   geom_point(shape=0.8, size=1)+
   theme_bw()+
   theme(legend.title=element_blank())

figfn <- paste0(figuredir,project,".2.3_pca.png")
png(figfn, width=520, height=420, res=120)
print(p1)
dev.off()

p1 <- ggplot(eigenvec2, aes(x=PC1, y=PC3, color=factor(Batch)))+
   geom_point(shape=0.8, size=1)+
   theme_bw()+
   theme(legend.title=element_blank())

figfn <- paste0(figuredir,project,".1.3_pca.png")
png(figfn, width=520, height=420, res=120)
print(p1)
dev.off()
