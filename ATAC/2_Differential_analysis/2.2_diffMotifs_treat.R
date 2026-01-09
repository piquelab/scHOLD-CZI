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
## parser$add_argument("-con", "--condition", type="character", default="CTRL",
##                     help="which group used for analysis",
##                     metavar="character")
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

## treat0 <- args$condition
psycho_file <- args$psycho_file

nc <- as.numeric(args$paralle_core)
register(MulticoreParam(nc))



## option <- "Cluster_res0.07"
node_grid <- "/rs/rs_grp_scatac/" 
path_wk <- paste(node_grid, "schold/ATAC/sc-atac-cziHOLD/", sep="")




###
###
outdir <- paste("./2_motif.outs/", option, "_th", th, "_treat/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)

###
### count file. we use th=0 and we will filter by th=20 in the later procedures  
infn_bulk <- paste("./2_motif.outs/", option, "/1.3_YtX_ave.th20.clean.rds", sep="")
## infn_combat <- paste("./0_pseudobulk.outs/", option, "_combat/", cluster0, "_combat_seq.adjusted.rds", sep="")
infn_ncell <- paste("./2_motif.outs/", option, "/0_ncell.rds", sep="")


###
### convariate files 
## infn_cv <- paste(path_wk, "HOLD-CZI_covariates_Ali_updated.txt", sep="")
## infn_batch <- paste(path_wk, "HOLD-CZI_covariates_Ali.dbgap.txt", sep="")
## infn_pc <- paste(path_wk, "HOLD-CZI_geno_pc_nokin.txt", sep="")
## infn_var <- paste("./psycho_var_dir/", psycho_file, sep="")
infn_cv <- "./psycho_var_dir/0.1_covariates_varSel.txt"


###
### counts data
counts <- read_rds(infn_bulk)


###
### convariate data

cvt <- as.data.frame(str_split(colnames(counts), "_",simplify=T))
names(cvt) <- c("Cluster", "treat", "sampleID")
cvt$bti <- colnames(counts)

cv <- read.table(infn_cv, header=T)
## psycho_vars <- read.table(infn_var)$V1

cvt_comb <- cvt%>%inner_join(cv, by="sampleID")


dd_ncell <- read_rds(infn_ncell)
bti2 <- dd_ncell%>%dplyr::filter(ncell>th)%>%pull(bti)
 

 



#######################
### run DESseq2
######################

###
### filtering conditions

cvt_comb <- cvt_comb%>%dplyr::filter(Cluster==cluster0, bti%in%bti2)



res_all <- NULL
summ_all <- NULL
treat0 <- c("CTRL", "LPS")



###
### apply filter  
cvt0 <- cvt_comb%>%drop_na(all_of(c("PC1", "PC2", "sex_alph", "age")))
###
summ_batch <- cvt0%>%group_by(Batch)%>%summarize(n0=length(unique(treat)), .groups="drop")%>%ungroup()
batch2 <- summ_batch%>%dplyr::filter(n0>1)%>%pull(Batch)%>%unique()
cvt0 <- cvt0%>%dplyr::filter(treat%in%treat0, Batch%in%batch2)
    
btiSel <- cvt0$bti
   
count0 <- counts[,btiSel]
motifs <- rownames(count0)


## design <- as.formula(paste0("~sex_alph + PC1 + PC2 + age + ", vari))
## identical(colnames(count0), cvt0$bti)
    
###
### running DESeq
time0 <- Sys.time()
    
cvti <- cvt0%>%dplyr::select(all_of(c("sampleID", "PC1", "PC2", "sex_alph", "Batch", "age", "treat")))

res <- NULL
nmotif_total <- nrow(count0) 
for ( k in 1:nmotif_total){
    ###
    y <- as.numeric(count0[k,])   
    cvti$y <- as.numeric(count0[k,])
       
    ## lm_full <- try(lm(y ~ 0 + factor(sex_alph) + factor(Batch) + PC1 + PC2 + age + factor(treat), data=cvti),
    ##                silent=F)
    lm_full <- try(lm(y ~ 0 + factor(sampleID) + factor(treat), data=cvti), silent=F)
           
    if ( class(lm_full) == "try-error") next 

    coef_df <- summary(lm_full)$coefficients
       
       ### slight difference between anova and loglik, detail explanations 
       ### https://stats.stackexchange.com/questions/155474/why-does-lrtest-not-match-anovatest-lrt
    irow <- nrow(coef_df)
    dd0 <- data.frame(Cluster = cluster0,
                         gene = motifs[k], 
                         beta = coef_df[irow, 1],
                         std = coef_df[irow, 2],
                         tval = coef_df[irow, 3],
                         pval = coef_df[irow, 4],
                         nind = length((unique(cvti$sampleID))),
                         treat = "LPS")
 
    res <- rbind(res, dd0)
    cat(cluster0, motifs[k], "\n")
    ###
    ###
}

###
###  padj 
res <- res%>%mutate(padj = p.adjust(pval, method="BH"))

    
time1 <- Sys.time()
elapsed <- difftime(time1, time0, units="mins")
    
cat("lm motif", ncol(count0), nrow(count0), elapsed, "mins\n")


###
### output 

opfn <- paste(outdir, cluster0, "_motifs.results.txt.gz", sep="")
fwrite(res, file=opfn, sep="\t", row.names=F, col.names=T, quote=F, na=NA) 


###
### End DESeq






                              

                         
