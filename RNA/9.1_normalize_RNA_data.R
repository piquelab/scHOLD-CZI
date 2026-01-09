####
####
library(tidyverse)
library(Matrix)
## library(DESeq2)
library(limma)
library(edgeR)
## library(biobroom)
library(data.table)

##
library(cowplot)
library(RColorBrewer)
library(scales)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(ggrastr)

###

rm(list=ls())

args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.1) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
outFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")

## change cluster 0-5 and run for each 
cluster <- "C5"


###
### passing arguments 
#args <- commandArgs(trailingOnly=T)
#if ( length(args) > 0){
#   ##
#   cluster <- as.character(args[1])
#}else{
#   cluster <- "C0"
#}    


outdir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"
#outdir <- "./0_normalized_data/"
#if ( !file.exists(outdir)) dir.create(outdir, showWarnings = F, recursive = T)

opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt_proteincoding.RData")
load(opfn)

count <- counts(counts_ls[[cluster]])
meta <- metadata_ls[[cluster]]
sampleIDsel <- unique(meta$Sample_ID)

# make sure its raw counts
#mean(as.matrix(counts) == 0)
#any(counts != floor(counts)) 

###
### gene counts data 
## rna_path <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/"
## rna_path <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/"
## fn_ls <- list.files(rna_path, "ALL.0.15.13.ComBat_seq")
## fn_ls <- fn_ls[grepl("_adjusted.RData", fn_ls)]

## ## df_file <- data.frame(file_name=fn_ls, cluster=gsub(".*ComBat_seq\\.|\\.SES.*", "", fn_ls))

## ## write.table(df_file, file="rna_cluster.infor.txt", quote=F, sep="\t", row.names=FALSE)
 
## df_file <- read.table("rna_cluster.infor.txt", header=T, sep="\t")


###
### covariates
#cv <- read.table("0.1_covariates_raw.txt", header=T, sep="\t")
cv <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/old/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0.1_covariates_raw.txt", header=T, sep="\t")
sampleIDsel <- cv$sampleID

#cv <- as.data.frame(cov_file)
#cv$sampleID <- cv$dbgap.ID
#cv$Sample_ID <- cv$dbgap.ID

# read in genotype PCs: 
#eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
#eigenvec2 <- merge(eigenvec2_o[,-c("sex","sex_alph","age")],meta,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
#meta_pc <- merge(cv, eigenvec2_o, by="Sample_ID")


#treat0 <- "CTRL" 
#treat0 <- "LPS" 

cl_rna <- cluster

###
### load counts data

#rna_path <- "noDex_raw_counts/"
#fn0 <- paste(rna_path, cl_rna, "_raw_counts.rds", sep="")
#count <- read_rds(fn0)



###
### column information 
col_RNA <- str_split(colnames(count), "_", simplify=T)
col_RNA <- data.frame(comb = colnames(count), cluster = col_RNA[, 1], #Batch = col_RNA[, 2],
      sampleID = col_RNA[,3], treat = gsub("RNA-", "", col_RNA[,4]))

#col_RNA <- col_RNA%>%filter(treat==treat0, sampleID%in%sampleIDsel)%>%inner_join(cv, by="sampleID")
col_RNA <- col_RNA%>%filter(sampleID%in%sampleIDsel)%>%inner_join(cv, by="sampleID")
summ_batch <- as.data.frame(table(col_RNA$Batch))
batch_sel <- summ_batch%>%filter(Freq>1)%>%pull(Var1)

col_RNA <- col_RNA%>%filter(Batch%in%batch_sel)



###
### filter genes and column 

count2 <- count[, col_RNA$comb]
nsample <- ncol(count2)

dge <- DGEList(count2)
cpm <- cpm(dge)
keep.exprs <- rowSums(cpm>0.1) >= 0.2*nsample
dge <- dge[keep.exprs, ]
    
dge <- calcNormFactors(dge, method = "TMM")
### cpm <- cpm(dge, log=T)
### For plotting and descriptive purposes we recommend cpm(x, log=TRUE) from Gordon Smyth


###
### normalize procedure and using two pcs. Previously I used the three PCs
### We didn't remove confounder effects

design <- model.matrix(~col_RNA$sex_alph + col_RNA$Batch + col_RNA$age + col_RNA$PC1 + col_RNA$PC2)
v <- voom(dge, design, plot=F)
geneExpr <- v$E


###
### save 
#opfn <- paste("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/norm_counts/", cluster, "_", treat0, "_rna.normal.rds", sep="")
opfn <- paste("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/norm_counts/alltreats/", cluster, "_rna.normal.rds", sep="")
write_rds(geneExpr, file=opfn)

# make sure its raw counts
#mean(as.matrix(geneExpr) == 0)
#any(geneExpr != floor(geneExpr)) 
#any(counts != floor(counts)) 

###
### END
