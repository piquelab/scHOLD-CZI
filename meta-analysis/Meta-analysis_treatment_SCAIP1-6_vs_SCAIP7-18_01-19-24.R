###
###
library(Matrix)
library(MASS)
library(scales)
library(tidyverse)
library(parallel)
library(data.table)
library(Rcpp)
library(reshape)
library(qqman)
library(qvalue)
##
library(Seurat)
library(DESeq2)
library(annotables)
library(biobroom)
###
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(gtable)
library(ggsignif)
library(pheatmap)
library(corrplot)
library(gtable)
library(RColorBrewer)
library(viridis)
library(ggrastr)
library(dplyr)
library(tidyr)
library(annotables)


##

rm(list=ls())

setwd("/rs/rs_grp_scaloft/8_meta_analysis/06-SCAIP16_SCAIP718_treatment/")

platePrefix <- "Meta_treatment_SCAIP16_SCAIP718_01-19-24"


#####################################################
####   meta analysis                              ###
###    SCAIP1-6 treatment results vs SCAIP7-18    ###
####    script by JW mdified by AR                ###
#####################################################

#############################
### self-defined function ###
#############################

###my function (1)
myMeta <- function(dx){
   sub0 <-  !(is.na(dx[,"log2FoldChange"])|is.na(dx[,"lfcSE"])) 
   if (any(sub0)){
            baseMean <- dx[sub0,"baseMean"]
            b <- dx[sub0,"log2FoldChange"]
            std <- dx[sub0,"lfcSE"]
      ###
            vb_i <- 1/std^2
            w <- vb_i/sum(vb_i)
            baseMeanHat <- sum(w*baseMean)
            bhat <- sum(w*b)
            sdhat <- sqrt(1/sum(vb_i))
            z <- bhat/sdhat
            p <- 2*pnorm(abs(z),lower.tail=F)
            tmp <- c(baseMeanHat, bhat, sdhat, p)
   }else{
            tmp <- c(NA, NA, NA, NA)
   }
   tibble(baseMeanHat=tmp[1], beta=tmp[2],stderr=tmp[3],pval=tmp[4])
}


### my function (2)
myqval <- function(pval){
   qval <- pval
   ii0 <- !is.na(pval)
   qval[ii0] <- qvalue(pval[ii0],pi0=1)$qvalues
   qval
}

#### assign data name, directories, and prefixes
data1 <- "SCAIP1-6-trt"
data2 <- "SCAIP7-18-trt"

dir1 <- ""
dir2 <- "/rs/rs_grp_scaloft/DGE_analysis/16_Treatment_SCAIP7-18/"

prefix1 <- ""
prefix2 <- "SCAIP7-18_treatment_01_17_24_"


#### colnames for the function: baseMean    log2FoldChange  lfcSE   stat    pvalue  padj

# SCAIP1-6 treatment results: 

### load SCAIP1-6 meta.data results 
# /nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
oldres <- read_rds(fn)

names(oldres)[names(oldres) == 'baseMeanHat'] <- 'baseMean'
names(oldres)[names(oldres) == 'beta'] <- 'log2FoldChange'
names(oldres)[names(oldres) == 'stderr'] <- 'lfcSE'
names(oldres)[names(oldres) == 'pval'] <- 'pvalue'
names(oldres)[names(oldres) == 'qval'] <- 'padj'

oldres$ensgene <- oldres$gene
oldres <- oldres %>% filter(contrast != "LPS-DEX")
oldres <- oldres %>% dplyr::select(MCls, contrast, rn, baseMean, log2FoldChange, lfcSE, pvalue, padj, ensgene)


# load SCAIP7-18 results: 
newres <- read.csv("/rs/rs_grp_scaloft/DGE_analysis/16_Treatment_SCAIP7-18/DESeq_results_all_treatment_SCAIP7-18_treatment_01_17_24.csv", row.names=NULL)

#names(newres)[names(newres) == 'baseMeanHat'] <- 'baseMean'
names(newres)[names(newres) == 'estimate'] <- 'log2FoldChange'
names(newres)[names(newres) == 'stderror'] <- 'lfcSE'
names(newres)[names(newres) == 'p.value'] <- 'pvalue'
names(newres)[names(newres) == 'p.adjusted'] <- 'padj'


# get ensemble geneIDs
grch38_unq <- grch38%>%
              distinct(ensgene,.keep_all=T)%>%
              dplyr::select(ensgene, symbol)

# get ensemble geneIDs if data1 or data2 is SCAIP7-18 (where its by gene symb)
            newres$symbol <- newres$gene
            newres <- newres %>% left_join(grch38_unq, by="symbol")
            #newres <- newres %>% distinct(symbol, .keep_all=TRUE)
            #newres <- newres[-7]
            newres$gene <- newres$ensgene
            newres$rn <- paste0(newres$MCls, "_", newres$contrast, "_", newres$gene)

newres <- newres %>% dplyr::select(MCls, contrast, rn, baseMean, log2FoldChange, lfcSE, pvalue, padj, ensgene)

colnames(oldres)
colnames(newres)
setdiff(colnames(oldres), colnames(newres))

############## assing treatment cell type, wave (if necessary before the loop), and filter covariate if necessary ####################
cell_types <- c("NKcell","Tcell","Bcell","Monocyte")

contrasts <- c("LPS", "PHA", "PHA-DEX")

#writing R output to this textfile
sink(paste0(platePrefix,"_output_1.txt"))

timestamp()




################################################# write the loop for meta-analysis 


for(i in 1:length(cell_types)){
    
    cell <- cell_types[i]

    res1Cell <- oldres %>% filter(MCls == cell)
    res2Cell <- newres %>% filter(MCls == cell)

            
    cat(" \n")
    cat("######################################################################################################\n")
    cat("################# Processing Cell Type: ", cell, "#######################################","\n")
    cat("################# Meta analysis for: ", data1, " and ", data2, "#################","\n")
    cat("################# plate prefeix: ", platePrefix, "##########","\n")

    for(i in 1:length(contrasts)){
            contr <- contrasts[i]
            #contr <- variable_names[i]

            ## assign variables, load data, and load experiment information
            topDirectory <- 'out_data_'
            outDir <- paste(topDirectory, platePrefix,"/",contr, sep='')
            system(paste("mkdir -p",outDir))
            ##
            plotsDir <- paste(outDir, '/plots', sep='')
            system(paste("mkdir -p", plotsDir))
            ##
            statsDir <- paste(outDir, '/stats', sep='')
            system(paste("mkdir -p", statsDir))
            ##
            dataDir <- paste(outDir, '/data_objects', sep='')
            system(paste("mkdir -p", dataDir))

            
            cat("########################################################\n")
            cat("##Processing: ", contr, "for", cell, "\n")

            # load DESeq results genes for data 1
            #resall1 <- read.table(paste0(dir1, contr, '/stats/', prefix1, '_', contr, '_', cell, '_DESeq_results.txt', sep=''), header=T)
            #resall2 <- read.table(paste0(dir2, contr, '/stats/', prefix2, '_', contr, '_', cell, '_DESeq_results.txt', sep=''), header=T)

            res1contr <- res1Cell %>% filter(contrast == contr)
            res2contr <- res2Cell %>% filter(contrast == contr)

            cat("## number of genes in ", data1, " for ", cell, " for cotrast ", contr, ": ",  nrow(res1contr),"\n")
            cat("## number of genes in ", data2, " for ", cell, " for cotrast ", contr, ": ",  nrow(res2contr),"\n")

            cat("\n")
            cat("# number of unique genes in ", data1, " for ", cell, " for cotrast ", contr, ": ",  length(unique(res1contr$ensgene)),"\n")
            cat("# number of unique genes in ", data2, " for ", cell, " for cotrast ", contr, ": ",  length(unique(res2contr$ensgene)),"\n")

            
            # concatenate the DESEq results
            rescontr <- rbind(res1contr, res2contr)

            # meta-analysis
            rescontrA <- rescontr%>%group_by(ensgene)%>%
                    nest()%>%
                    mutate(outlist=mclapply(data, myMeta, mc.cores=1))%>%  #mutate(outlist=future_map(data,myMeta))%>%
                    dplyr::select(-data)%>%unnest(c(outlist))%>%as.data.frame()           
                       
            # add qvalue
            rescontrB <- rescontrA%>%
                    mutate(qval=myqval(pval))

            rescontrB_sig <- rescontrB %>% filter(qval < 0.1)

            cat("Meta_analysis number of genes across both data in ", cell, " ", contr, ": ", nrow(rescontrB), "\n") 
            #cat("Meta_analysis DEGs: ", nrow(resB_NK %>% filter(qval < 0.1)), "\n") 
            cat("Meta_analysis DEGs: ", nrow(rescontrB_sig), "\n") 

            fname=paste(statsDir, '/', contr, "_", cell, "_", platePrefix, "_meta_analysis_results", ".txt", sep="")
            write.table(rescontrB, file=fname, quote=F, sep="\t")

            fname=paste(statsDir, '/', contr, "_", cell, "_", platePrefix, "_meta_analysis_significant_results", ".txt", sep="")
            write.table(rescontrB_sig, file=fname, quote=F, sep="\t")
    }
}






############################# write a loop to save the results in a table

sink()
sink()
sink()

#setwd("/rs/rs_grp_scaloft/DGE_analysis/11_DESeq_loop_SCAIP7-18_W-A/")

############## assing treatment cell type, wave (if necessary before the loop), and filter covariate if necessary ####################
cell_types <- c("NKcell","Tcell","Bcell","Monocyte")
contrasts <- c("LPS", "PHA", "PHA-DEX")


#platePrefix <- "SCAIP7-18_W-A_DESeq_new_cov_01_03_24"

sink(paste0(platePrefix,"_stats_smmary_out.txt"))

timestamp()

tbl <- data.frame(matrix(ncol=0, nrow=3))

for(i in 1:length(cell_types)){
    
    cell <- cell_types[i]
            
    cat(" \n")
    cat("######################################################################################################\n")
    cat("############################### Processing Cell Type: ", cell, "#######################################","\n")

   
    cat("Now selecting count data for cell type: ", cell, "!\n")

   
    table <- data.frame(matrix(ncol=3, nrow=0))
    x <- c("contrast", paste(cell, "_ngenes", sep=''), paste(cell, "_DEGs", sep=''))
    colnames(table) <- x


    # DESeq loop for each variable
        for(i in 1:length(contrasts)){
            contr <- contrasts[i]
            
            cat("########################################################\n")
            cat("##Processing: ", contr, "\n")

            topDirectory <- 'out_data_'
            outDir <- paste(topDirectory, platePrefix,"/",contr, sep='')
            ##
            statsDir <- paste(outDir, '/stats', sep='')
            ##
            dataDir <- paste(outDir, '/data_objects', sep='')

            cat("Now testing ", contr, " for cell type: ", cell, "!\n")
            cat(" \n")


          #load result
          fname=paste(statsDir, '/', contr, "_", cell, "_", platePrefix, "_meta_analysis_results", ".txt", sep="")
          resVar <- read.table(fname, header=T)

          resVar_sig <- resVar %>% filter(qval <0.1,abs(beta)>0.5) %>%drop_na(beta)         
          resVar_sig <- resVar_sig[order(resVar_sig$qval),]

          #cat("Number of individuals in variable", contr , " in cell type ", cell," :", length(unique(ddsVar@colData$dbgap.ID)), "\n")
          cat("# of genes tested for - ", contr, " in cell type ", cell," :", nrow(resVar)," \n")
          cat(" \n")

          cat("# of DEGs for - ", contr,  " in celltype ", cell," (FDR < 0.1) : ", nrow(resVar_sig),"\n")

            tableVar <- data.frame(matrix(ncol=3, nrow=1))
            x <- c("contrast", paste(cell, "_ngenes", sep=''), paste(cell, "_DEGs", sep=''))
            colnames(tableVar) <- x
            tableVar[1] <- paste(contr)
            tableVar[2] <- paste(nrow(resVar))
            tableVar[3] <- paste(nrow(resVar_sig))


           table <- rbind(table, tableVar)
        }

        # write table for the cell type
        tbl <- cbind(tbl, table)

    }         



    fname=paste("stats_table_all_treatments_per_cell_type_genes_degs", "_", platePrefix,".txt", sep="")
    write.table(tbl, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

   