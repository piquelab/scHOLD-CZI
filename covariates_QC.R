#a variation of the following script was used to convert spss file to R file and convert participant id to dbgap id
#/rs/rs_grp_scaloft/scALOFT_2024/covariates/make_cov/make_covariate_03-30-24.R
#Please note that subject 6055 should be dropped from the analyses (final N = 210) - was not in current list of dbgap id possibilities and so is not in the file. may need to be removed in the future if that changes

library(data.table)
library(ggplot2)
library(plyr)
#library(grid)
#library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/","/rs/rs_grp_schold/covariates/","HOLD-CZI_covariates_converted_dbgapID_key_n183_AR.txt","CZI_expectedrange_042024.txt") #for testing
base <- args[1]
analysis=paste0(args[1],"RNA/analysis/")
covfolder=args[2]
cov_file=paste0(covfolder,args[3])
cov_expected_file=paste0(covfolder,args[4])

figuredir=paste0(base,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

cov <- fread(cov_file)

#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
eigenvec2 <- fread(file=paste0(analysis,"genotypePC/eigenvec_pc.txt")) #will use col PC1

cov_sub <- cov[cov$dbgap.ID %in% eigenvec2$Sample_ID,]

col_names <- colnames(cov)
col_names <- col_names[-1]

for (i in col_names){
	cat("plotting ",i)
	png(paste0(figuredir,"hist_",i,".png"), width=500, height=500, res=120)
	plot <- hist(na.omit(cov[[i]]),main="",xlab=i,breaks=15)
    dev.off()
}

for (i in col_names){
	cat("plotting ",i)
	png(paste0(figuredir,"CZI1_hist_",i,".png"), width=500, height=500, res=120)
	plot <- hist(na.omit(cov_sub[[i]]),main="",xlab=i,breaks=15)
    dev.off()
}

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

CZI1_table <- fread(cov_expected_file)
min <- list()
max <- list()
for (i in col_names){
	min[[i]] <- round(min(cov_sub[[i]],na.rm=T),3)
	max[[i]] <- round(max(cov_sub[[i]],na.rm=T),3)
}
df <- cbind(ldply(min, data.frame),ldply(max, data.frame)[,2])
colnames(df) <- c("variable","data_min","data_max")

CZI1_table_c <- merge(CZI1_table,df,by="variable")

fwrite(CZI1_table_c, sep='\t', quote=F, row.names=F, col.names=T, paste0(base,"CZI1_range.txt"),append=T)

#################################3
####################################
#plot_list <- list()

#for (i in col_names){
#	plot <- histylim(na.omit(cov_sub[[i]]),main="",xlab=i,breaks=15)
#    plot_list[[i]] <- plot
#}

#png(paste0(figuredir,"CZI1_hist_allcols.png"), width=1000, height=1000, res=120)
#plot_grob <- arrangeGrob(grobs=plot_list)
#grid.arrange(plot_grob)
#dev.off()