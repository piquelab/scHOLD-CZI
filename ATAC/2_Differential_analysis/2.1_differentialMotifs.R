###
library(tidyverse)
library(Matrix)
## library(biobroom)
## library(DESeq2)
library(sva)
library(data.table)
library(BiocParallel)
library(argparse)

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

## parser$add_argument("-adj", "--adj_batch", type="character", default="combat_seq",
##                         help="Adjust batch effects for counts data",
##                         metavar="character")

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
## combat_adjust <- args$adj_batch
option <- args$option
th <- as.numeric(args$filter_th)
treat0 <- args$condition
psycho_file <- args$psycho_file

nc <- as.numeric(args$paralle_core)
register(MulticoreParam(nc))



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


node_grid <- "/rs/rs_grp_scatac/" 
path_wk <- paste(node_grid, "schold/ATAC/sc-atac-cziHOLD/", sep="")


## path_wk <- "/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/"
### working path
### `/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis`



###
###
outdir <- paste("./2_motif.outs/", option, "_th", th, "_psycho/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
### count file
infn_bulk <- paste("./2_motif.outs/", option, "/1.3_YtX_ave.th", th, ".clean.rds", sep="")
## infn_combat <- paste("./0_pseudobulk.outs/", option, "_combat/", cluster0, "_combat_seq.adjusted.rds", sep="")
infn_ncell <- paste("./2_motif.outs/", option, "/0_ncell.rds", sep="")


###
### convariate files 
## infn_cv <- paste(path_wk, "HOLD-CZI_covariates_Ali_updated.txt", sep="")
## infn_batch <- paste(path_wk, "HOLD-CZI_covariates_Ali.dbgap.txt", sep="")
## infn_pc <- paste(path_wk, "HOLD-CZI_geno_pc_nokin.txt", sep="")
infn_cv <- "./psycho_var_dir/0.1_covariates_varSel.txt" 
infn_var <- paste("./psycho_var_dir/", psycho_file, sep="") 


###
### counts data

counts <- read_rds(infn_bulk)


###
### convariate data

cvt <- as.data.frame(str_split(colnames(counts), "_",simplify=T))
names(cvt) <- c("Cluster", "treat", "sampleID")
cvt$bti <- colnames(counts)

cv <- read.table(infn_cv, header=T, sep="\t")
psycho_vars <- read.table(infn_var)$V1


cvt_comb <- cvt%>%inner_join(cv, by="sampleID")

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


## mm <- sapply(1:ncol(cv), function(i) mode(cv[,i])) 

###
### combine covariates
## cvt_comb <- cvt%>%inner_join(pcs, by="sampleID")%>%inner_join(cv, by="sampleID")


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
cvt_comb <- cvt_comb%>%dplyr::filter(Cluster==cluster0, bti%in%bti2) ###, Batch!="HOLD07")

 
res_all <- NULL
summ_all <- NULL
for (vari in psycho_vars){
   ### 

   ###
   ### apply filter  
   cvt0 <- cvt_comb%>%drop_na(all_of(c(c("PC1", "PC2", "sex_alph", "age"), vari)))
   ###
   ## Batch_filter=FALSE ## default TRUE 
   ## if ( Batch_filter ){
   ##    ### 
   ##    summ_batch <- cvt0%>%group_by(Batch)%>%summarize(n0=length(unique(treat)), .groups="drop")%>%ungroup()
   ##    batch2 <- summ_batch%>%dplyr::filter(n0>1)%>%pull(Batch)%>%unique()
   ## }else{
   ##    ###
   ##    summ_batch <- cvt0%>%filter(treat==treat0)%>%group_by(Batch)%>%summarize(n0=n(), .groups="drop")%>%ungroup()
   ##    batch2 <- summ_batch%>%dplyr::filter(n0>1)%>%pull(Batch)%>%unique() 
   ## }
   summ_batch <- cvt0%>%filter(treat==treat0)%>%group_by(Batch)%>%summarize(n0=n(), .groups="drop")%>%ungroup()
   batch2 <- summ_batch%>%dplyr::filter(n0>1)%>%pull(Batch)%>%unique() 
   ##
   ## 
   cvt0 <- cvt0%>%dplyr::filter(treat==treat0, Batch%in%batch2)  
   btiSel <- cvt0$bti
   
   count0 <- counts[, btiSel]
   motifs <- rownames(count0)

   
   ## design <- as.formula(paste0("~sex_alph + PC1 + PC2 + age + ", vari))
 
   ###
   ### running DESeq
   time0 <- Sys.time()
    
   cvti <- cvt0%>%dplyr::select(all_of(c("PC1", "PC2", "sex_alph", "Batch", "age")))
   cvti$z <- cvt0%>%dplyr::pull(vari)        
   res <- NULL
   for ( k in 1:nrow(count0)){
       ###
       y <- as.numeric(count0[k,])   
       cvti$y <- as.numeric(count0[k,])

       lm0 <- try(lm(y ~ 0 + factor(sex_alph) + factor(Batch) + PC1 + PC2 + age + z, data=cvti), silent=T)

       if ( class(lm0) == "try-error") next
       
       ###
       bhat <- coef(lm0)["z"]
       sdhat <- sqrt(diag(vcov(lm0))["z"])
       zval <- bhat/sdhat
       p <- 2*pnorm(-abs(zval))


       ##
       coef_df <- summary(lm0)$coefficients["z",]
       ## tval <- coef_df[3] 
       ## pval_t <-  coef_df[4]

       dd0 <- data.frame(Cluster = cluster0,
                         gene = motifs[k], estimate = bhat, stderror = sdhat, statistic = zval, p.value = p,
                         pval_t = coef_df[4],               
                         psycho_variable = vari,
                         treat = treat0,
                         nind = nrow(cvti))
                         
       res <- rbind(res, dd0)
    }
   
   res <- res%>%mutate(p.adjusted=p.adjust(p.value, method="BH"),
                       padj_t = p.adjust(pval_t, method="BH"))
   res_all <- rbind(res_all, res) 
   
 
   ### summary
   summ <- data.frame(Cluster=cluster0, variable=vari, treat=treat0,
       nind=nrow(cvt0),  ngene_test=nrow(res),
       nsig_fdr0.1=sum(res$p.adjusted<0.1, na.rm=T),
       nsig_fdr0.1_t=sum(res$padj_t<0.1, na.rm=T))
    
   summ_all <- rbind(summ_all, summ) 

   time1 <- Sys.time()
   elapsed <- difftime(time1, time0, units="mins")
    
   cat("lm motif", vari, ncol(count0), nrow(count0), elapsed, "\n")

}
 
###
### save object
opfn2 <- paste(outdir, cluster0, "_", treat0, "_", gsub("\\.txt", "", psycho_file), "_motifs.summ.tsv", sep="")
write_tsv(summ_all, opfn2)


###
###
opfn <- paste(outdir, cluster0, "_", treat0, "_", gsub("\\.txt", "", psycho_file), "_comb_motifs.results.txt.gz", sep="")
fwrite(res_all, file=opfn, sep="\t", row.names=F, col.names=T, quote=F, na=NA) 


###
### End DESeq


##############################
### visualization results
##############################





                              

                         
