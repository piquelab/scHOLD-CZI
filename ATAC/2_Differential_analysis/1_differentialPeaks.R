### 
library(tidyverse)
library(Matrix)
library(biobroom)
library(DESeq2)
library(limma)
library(edgeR)
library(zingeR)
library(sva)
library(data.table)
library(BiocParallel)
library(argparse)
library(SummarizedExperiment)
rm(list=ls())



## time0 <- Sys.time()

###
### passing arguments 
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
                        dest="verbose", help="Print little output")

parser$add_argument("-cl", "--cluster", type="character", default="C0",
                        help="Cell type or cluster",
                        metavar="character")

parser$add_argument("-adj", "--adj_batch", type="character", default="no_adjust",
                        help="Adjust batch effects for counts data",
                        metavar="character")

parser$add_argument("-opt", "--option", type="character", default="option_nFeature15K_cluster_res0.12",
                        help="Different approach to define cluster",
                        metavar="character")
parser$add_argument("-th", "--filter_th", type="character", default="20",
                    help="filtering conditions,#cells/combination",
                    metavar="number")
parser$add_argument("-cond", "--condition", type="character", default="CTRL",
                    help="which group used for analysis",
                    metavar="character")
parser$add_argument("-psycho", "--psycho_file", type="character", default="psychosocial_top10.txt",
                    help="psychosocial variable",
                    metavar="character")

parser$add_argument("-core", "--paralle_core", type="character", default="1",
                    help="cores required in parallel setting",
                    metavar="number") 


##
args <- parser$parse_args()
cluster0 <- args$cluster
combat_adjust <- args$adj_batch
option <- args$option
th <- as.numeric(args$filter_th)
treat0 <- args$condition
psycho_file <- args$psycho_file

nc <- as.numeric(args$paralle_core)
register(MulticoreParam(nc))

## psycho_file <- "psychosocial_top10.txt"

 
## args <- commandArgs(trailingOnly=T)
## if ( length(args)>0){
##    ###
##    option <- args[1]
##    th <- args[2]
##    cluster0 <- args[3]
##    psycho_file <- args[4] 
##    treat0 <- args[5]
##    combat_adjust <- args[6] 
## }else{
##   option <- "Cluster_res0.07"
##   th <- 20
##   cluster0 <- "C0"
##   psycho_file <- "psychosocial_variables"
##   treat0 <- "CTRL"
##   combat_adjust <- "combat_seq"
## }    

## option <- "Cluster_res0.07"
node_grid <- "/rs/rs_grp_scatac/" 
path_wk <- paste(node_grid, "schold/ATAC/sc-atac-cziHOLD/", sep="")




###
###
outdir <- paste("./1_DiffPeak.outs/", option, "_th", th, "_psycho/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
### count file
infn_bulk <- paste("./0_pseudobulk.outs/", option, "/2_YtX_comb.th", th, ".clean.rds", sep="")
### infn_combat <- paste("./0_pseudobulk.outs/", option, "_combat_th", th, "_allBatch/", cluster0, "_combat_seq.adjusted.rds", sep="")
infn_ncell <- paste("./0_pseudobulk.outs/", option, "/1_dd.ncell.rds", sep="")


###
### convariate files 
infn_cv <- "./psycho_var_dir/0.1_covariates_varSel.txt"
infn_var <- paste("./psycho_var_dir/", psycho_file, sep="") 


###
### counts data

if ( combat_adjust=="combat_seq"){
###    
if ( file.exists(infn_combat)) {
   counts <- read_rds(infn_combat)
   cat(cluster0, "Combat-seq adjust counts", "\n")
}else{
   cat("Combat-seq file doesn't exist and please use combat-seq to correct counts data", "\n")
   stop("combat-seq File not exist")
}
}
###
###
    
if ( combat_adjust=="no_adjust"){
###    
if ( file.exists(infn_bulk)) {
   counts <- read_rds(infn_bulk)
   cat(cluster0, "use the raw counts data", "\n")
}else{
   cat("the raw count file doesn't exist and please create it first", "\n")
   stop("raw File not exist")
}
}



  
cvt <- as.data.frame(str_split(colnames(counts), "_", simplify=T))
names(cvt) <- c("Cluster", "treat", "sampleID")
cvt$bti <- colnames(counts)


cv <- read.table(infn_cv, header = T, sep = "\t")
psycho_vars <- read.table(infn_var)$V1


## batch_df <- read.table(infn_batch, header=T)%>%dplyr::select(sampleID=dbgap.ID, Batch)
## cvt <- cvt%>%left_join(batch_df, by="sampleID")

## ###
## ### genotype pcs
## pcs <- read.table(infn_pc, header=T)
## pcs <- pcs%>%dplyr::select(sampleID=Sample_ID, PC1, PC2, PC3)
 

## ###
## ### covariates
## cv <- read.table(infn_cv, header=T)

## ### psychosocial variables
## psycho_vars <- read.table(infn_var)$V1

## cv <- cv%>%dplyr::select(sampleID=dbgap.ID, age, sex_alph, all_of(psycho_vars))

## ## mm <- sapply(1:ncol(cv), function(i) mode(cv[,i])) 

## ###
## ### combine covariates
## cvt_comb <- cvt%>%inner_join(pcs, by="sampleID")%>%inner_join(cv, by="sampleID")

cvt_comb <- cvt%>%inner_join(cv, by = "sampleID")

dd_ncell <- read_rds(infn_ncell)
bti2 <- dd_ncell%>%dplyr::filter(ncell>th)%>%pull(bti)



### get cluster
## opfn <- paste("./Cluster/", option, ".cl.txt", sep="")
## write.table(sort(unique(cvt$Cluster)), file=opfn, quote=F, row.names=F, col.names=F)
 



#######################
### run DESseq2
######################

###
### filtering conditions
cvt_comb <- cvt_comb%>%dplyr::filter(Cluster==cluster0, bti%in%bti2)

max_grid <- .Machine$integer.max
### DESeq model~PC1+PC2+sex_alph+age+var
 
## res_all <- NULL
summ_all <- NULL
for (vari in psycho_vars){
   ### 

   ###
   ### apply filter  
   cvt0 <- cvt_comb%>%drop_na(all_of(c(c("PC1", "PC2", "sex_alph", "age"), vari)))
   ###
   ## Batch_filter=TRUE ### default TRUE 
   ## if ( Batch_filter ){
   ##    ### 
   ##    summ_batch <- cvt0%>%group_by(Batch)%>%summarize(n0=length(unique(treat)), .groups="drop")%>%ungroup()
   ##    batch2 <- summ_batch%>%dplyr::filter(n0>1)%>%pull(Batch)%>%unique()
   ## }else{
   ##    summ_batch <- cvt0%>%dplyr::filter(treat==treat0)%>%group_by(Batch)%>%summarize(n0=n(), .groups="drop")%>%ungroup()
   ##    batch2 <- summ_batch%>%dplyr::filter(n0>1)%>%pull(Batch)%>%unique() 
   ## }
    
   ###
   ### Filter Batch  
   summ_batch <- cvt0%>%dplyr::filter(treat==treat0)%>%group_by(Batch)%>%summarize(n0=n(), .groups="drop")%>%ungroup()
   batch2 <- summ_batch%>%dplyr::filter(n0>1)%>%pull(Batch)%>%unique()
   cvt0 <- cvt0%>%dplyr::filter(treat==treat0, Batch%in%batch2) 
   btiSel <- cvt0$bti
   
   count0 <- counts[, btiSel]

   ### 
   ### we will update this filtering threshold in later analysis 
   ### filter cpm final use 
   ## th_cpm <- 0.5 
   ## dge <- DGEList(count0)
   ## cpm <- cpm(dge) 
   ## rnz <- rowMeans(cpm>th_cpm)
   ## count0 <- count0[rnz>0.5, ] 

    
   ## filter peak>0  
   ## rnz <- rowSums(count0)
   ## count0 <- count0[rnz>0, ]
    
   ### filter peaks expressed 25% individuals, the same to RNA   
   rnz <- rowMeans(count0>0) 
   count0 <- count0[rnz>0.25, ]
    

   rmax <- rowSums(count0>max_grid)
   count0 <- count0[rmax==0,]

   cat("drop features because of large reads", sum(rmax>0), "\n")
   cat("variable", vari, "Features", nrow(count0), "Individuals", ncol(count0), "\n") 
   
   ## design <- as.formula(paste0("~sex_alph + PC1 + PC2 + age + ", vari))
   design <- as.formula(paste0("~sex_alph + Batch + PC1 + PC2 + age + ", vari))

   ###
   ### running DESeq
   time0 <- Sys.time()
   ##
   dds <- DESeqDataSetFromMatrix(count0, colData=cvt0, design=design)       
   dd_fn <- paste(outdir, cluster0, "_", treat0, "_", vari, "_DESeq.dds.rds", sep="")

      
   ##
   if ( file.exists(dd_fn) & file.size(dd_fn)>0){
       ##
       cat("Already have DESeq object", "\n")
       dds <- read_rds(dd_fn)
   }    

    
   ### DESeq2 
   if ( ! ( file.exists(dd_fn) & (file.size(dd_fn)>0) )){
      ###
      cat("Run DESeq", "\n")
      dds <- try(DESeq(dds, parallel=T, BPPARAM=MulticoreParam(nc)), silent=T)
      write_rds(dds, file=dd_fn)
   }

       
    
   if ( class(dds)=="try-error" ) next
    

###
   ### results
   res <- try(results(dds), silent=T)
   if ( class(res)=="try-error") next
    
   res <- tidy(res)
   res$psycho_variable <- vari
   res$Cluster <- cluster0
   res$treat <- treat0
   res$nind <- length(unique(cvt0$sampleID))
    
   opfn <- paste(outdir, cluster0, "_", treat0, "_", vari, "_DESeq.results.txt.gz", sep="")
   fwrite(res, file=opfn, sep="\t", row.names=F, col.names=T, quote=F, na=NA) 
 
   ### summary
   summ <- data.frame(Cluster=cluster0, variable=vari, treat=treat0,
       nind=nrow(cvt0),  ngene_test=nrow(res), ngene_ind=sum(!is.na(res$p.adjusted)),
       nsig_fdr0.1=sum(res$p.adjusted<0.1, na.rm=T))
 

   time1 <- Sys.time()
   elapsed <- difftime(time1, time0, units="mins")
    
   cat("DESeq", vari, ncol(count0), nrow(count0), elapsed, "\n")
   summ_all <- rbind(summ_all, summ)
}

###
### save object
opfn2 <- paste(outdir, cluster0, "_", treat0, "_", gsub("\\.txt", "", psycho_file), "_DESeq.summ.tsv", sep="")
write_tsv(summ_all, opfn2)

###
### End DESeq


##############################
### visualization results
##############################





                              

                         
