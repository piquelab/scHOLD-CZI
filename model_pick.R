#selecting variables for CZI model, originally picked based on correlation and then logic. but failed to take into account number of individuals missing data.
#now subsetting for number missing samples, then correlation, then logic

library(sva)
require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(tidyverse)
library(data.table)
library(plyr)
library(cowplot)
library(parallel)
library(ggpubr)
future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)
args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_converted_dbgapID_key_n183_AR_18newvar_updated_08-05.txt","ALL","fastdemux") #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip")
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]

opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".DESeq_countlists.RData")
load(opfn)

cluster="C0"
shortvars <- c("FCDEM_12", "smoke", "CVDRISK", "PSQI_total", "SES", "cytocomp", "pr_comp", "StressSev", "SNI_NoP", 
	"DED_all_mean", "Chol", "BPd_avg")
cluster_metadata_sce <- metadata_ls[[cluster]]
cluster_metadata <- data.frame(cluster_metadata_sce)
toremove <- data.frame(ind=c("HO-047","HO-110","HO-116","HO-117","HO-118"))
cluster_metadata <- subset(cluster_metadata, treats=="RNA-CTRL" & !Sample_ID %in% toremove$ind)

##c <- apply(cluster_metadata[,c("treats","sex_alph","age","FCDEM_12", "smoke", "CVDRISK", "PSQI_total", "SES", "cytocomp", "pr_comp", "StressSev")], MARGIN = 1, function(x) sum(is.na(x)))
#c <- apply(cluster_metadata[, c("treats","sex_alph","age",shortvars)], MARGIN = 2, function(x) sum(is.na(x)))
#c2 <- apply(cluster_metadata[, c("treats","sex_alph","age",shortvars)], MARGIN = 2, function(x) sum(!is.na(x)))
#cdf <- data.frame(variable=rownames(as.data.frame(c)),number_samples=as.data.frame(c2)$c,number_missing_samples=as.data.frame(c)$c)

c <- apply(cluster_metadata[,c("treats","sex_alph","age",psychvarstorun)], MARGIN = 2, function(x) sum(is.na(x)))
c2 <- apply(cluster_metadata[,c("treats","sex_alph","age",psychvarstorun)], MARGIN = 2, function(x) sum(!is.na(x)))
cdf <- data.frame(variable=rownames(as.data.frame(c)),number_samples=as.data.frame(c2)$c,number_missing_samples=as.data.frame(c)$c)
cdf_s <- cdf[cdf$number_missing_samples<5,]
head(cdf_s[order(cdf_s$number_missing_samples),])
cdf_nomiss <- cdf[cdf$number_missing_samples<1,]
cdf_nomiss <- subset(cdf_nomiss, !variable %in% c("DSES_07","sex_alph","treats","age"))
fwrite(as.data.frame(unique(cdf_nomiss$variable)), file=paste0(base,"COMBAT_modelvars_nomiss.txt"), sep='\t', quote=F, row.names=F, col.names=F)


load("/rs/rs_grp_schold/covariates/other_covariates/correlation/variable_score_correlations_n51.RData")
m <- reshape2::melt(cor_matrix)
m <- subset(m, !Var1==Var2)
m_h <- subset(m, !Var1 %in% c(unique(cdf_s$variable),notrun_var) & !Var2 %in% c(unique(cdf_s$variable),notrun_var))
high_cor <- subset(m_h, value>0.5)
corremove <- c("LogCRP","EDS_mean","Chol_LDL","Chol_Ratio","StressSev")
varremove <- c("Logifny","Logil10","Logil12","Logil13","Logil1b","Logil2","Logil4","Logil6","Logil8","Logtnfa","DSES_07","DSES_09","NAIscr2019","BPs_avg")
keep <- c("HS_CRP","HVS_mean","Chol","Trig","StressCount")
vars <- psychvarstorun[psychvarstorun %in% cdf_s$variable & !psychvarstorun %in% varremove]
fwrite(as.data.frame(vars), file=paste0(base,"COMBAT_modelvars.txt"), sep='\t', quote=F, row.names=F, col.names=F)

#check numbers
cluster_metadata_var1 <- cluster_metadata[,c("Sample_ID","PC1","PC2","sex_alph","age","treats",shortvars)]
cluster_metadata_var1 <- cluster_metadata_var1[complete.cases(cluster_metadata_var1), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
cluster_metadata_var <- cluster_metadata[,c("Sample_ID","PC1","PC2","sex_alph","age","treats",vars)]
cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs

dim(cluster_metadata[,c("Sample_ID","PC1","PC2","sex_alph","age","treats")]) #161
dim(cluster_metadata_var) #149
dim(cluster_metadata_var1) #117