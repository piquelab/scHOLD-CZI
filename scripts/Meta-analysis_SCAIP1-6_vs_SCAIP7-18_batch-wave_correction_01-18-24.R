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

setwd("/rs/rs_grp_scaloft/8_meta_analysis/03-SCAIP16_SCAIP718_bw_corrected/")

platePrefix <- "Meta_SCAIP16_SCAIP718_bw_correct_01-18-24"


#####################################################
####   meta analysis                              ###
###    SCAIP1-6 results vs SCAIP7-18              ###
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
data1 <- "SCAIP1-6-bw-correct"
data2 <- "SCAIP7-18-bw-correct"

dir1 <- "/rs/rs_grp_scaloft/DGE_analysis/10_DESeq_loop_SCAIP1-6/cov_fixed_age_genPC/out_data_SCAIP1-6_DESeq_new_cov_01_10_24/"
dir2 <- "/rs/rs_grp_scaloft/DGE_analysis/15_DESeq_loop_SCAIP7-18/03_batch-wave_corrected/out_data_SCAIP7-18_batch-wave_corrected_01_17_24/"

prefix1 <- "SCAIP1-6_DESeq_new_cov_01_10_24"
prefix2 <- "SCAIP7-18_batch-wave_corrected_01_17_24"

############## assing treatment cell type, wave (if necessary before the loop), and filter covariate if necessary ####################

# assign variabeles: needs modifications
variables <- c("pedu", "pincme", "psesl", "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasgm", "csinhl", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq" 
                )
 
variable_names <- c("Parental Education", "Parental Income","Subjective SES", "Neighborhood Stress", "Self-disclosure", 
                "Perceived responsiveness", "Parental Warmth", "YR Parent Child Conflict",  "Risky Family", "Youth Depression", 
                "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Stimulated IL-5 cortisol trt", "Stimulated IL-13 cortisol trt", "Stimulated IFNy cortisol trt",    
                "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", 
                "DD Asthma symptoms", "Nightly Asthma", "Nightly inhaler use", "Asthma Frequency", "Asthma Severity", 
                "Nightly Awakenings",  "Sleep Duration", "Sleep quality"
                )


#variables <- c("IFNG_co", "pnsi")
#variable_names <- c("Stimulated IFNy", "Neighborhood Stress")

#treatment <- "CTRL"
cell_types <- c("NKcell","Tcell","Bcell","Monocyte")


# get ensemble geneIDs
grch38_unq <- grch38%>%
              distinct(ensgene,.keep_all=T)%>%
              dplyr::select(ensgene, symbol)

getwd()
platePrefix


#writing R output to this textfile
sink(paste0(platePrefix,"_output_1.txt"))

timestamp()




################################################# write the loop for meta-analysis 


for(i in 1:length(cell_types)){
    
    cell <- cell_types[i]
            
    cat(" \n")
    cat("######################################################################################################\n")
    cat("################# Processing Cell Type: ", cell, "#######################################","\n")
    cat("################# Meta analysis for: ", data1, " and ", data2, "#################","\n")
    cat("################# plate prefeix: ", platePrefix, "##########","\n")

    for(i in 1:length(variables)){
            myvar <- variables[i]
            var_name <- variable_names[i]

            ## assign variables, load data, and load experiment information
            topDirectory <- 'out_data_'
            outDir <- paste(topDirectory, platePrefix,"/",myvar, sep='')
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
            cat("##Processing: ", var_name, "for", cell, "\n")

            # load DESeq results genes for data 1
            resall1 <- read.table(paste0(dir1, myvar, '/stats/', prefix1, '_', myvar, '_', cell, '_DESeq_results.txt', sep=''), header=T)
            resall2 <- read.table(paste0(dir2, myvar, '/stats/', prefix2, '_', myvar, '_', cell, '_DESeq_results.txt', sep=''), header=T)

            cat("## number of genes in ", data1, " for ", cell, " for variable ", var_name, ": ",  nrow(resall1),"\n")
            cat("## number of genes in ", data2, " for ", cell, " for variable ", var_name, ": ",  nrow(resall2),"\n")

            # remove the decimal from ensgene on data1 or data2
            resall1$ensgene <- substr(rownames(resall1), 0, 15)

            # get ensemble geneIDs if data1 or data2 is SCAIP7-18 (where its by gene symb)
            resall2$symbol <- rownames(resall2)
            resall2 <- resall2 %>% left_join(grch38_unq, by="symbol")
            resall2 <- resall2 %>% distinct(symbol, .keep_all=TRUE)
            resall2 <- resall2[-7]
            cat("## post conversion, number of genes in ", data2, " for ", cell, " for variable ", var_name, ": ",  nrow(resall2),"\n")

            # concatenate the DESEq results
            resall <- rbind(resall1, resall2)

            # meta-analysis
            resallA <- resall%>%group_by(ensgene)%>%
                    nest()%>%
                    mutate(outlist=mclapply(data, myMeta, mc.cores=1))%>%  #mutate(outlist=future_map(data,myMeta))%>%
                    dplyr::select(-data)%>%unnest(c(outlist))%>%as.data.frame()           
                       
            # add qvalue
            resallB <- resallA%>%
                    mutate(qval=myqval(pval))

            resallB_sig <- resallB %>% filter(qval < 0.1)

            cat("Meta_analysis number of genes across both data in ", cell, " ", var_name, ": ", nrow(resallB), "\n") 
            #cat("Meta_analysis DEGs: ", nrow(resB_NK %>% filter(qval < 0.1)), "\n") 
            cat("Meta_analysis DEGs: ", nrow(resallB_sig), "\n") 

            fname=paste(statsDir, '/', myvar, "_", cell, "_", platePrefix, "_meta_analysis_results", ".txt", sep="")
            write.table(resallB, file=fname, quote=F, sep="\t")

            fname=paste(statsDir, '/', myvar, "_", cell, "_", platePrefix, "_meta_analysis_significant_results", ".txt", sep="")
            write.table(resallB_sig, file=fname, quote=F, sep="\t")
    }
}






############################# write a loop to save the results in a table

sink()
sink()
sink()

#setwd("/rs/rs_grp_scaloft/DGE_analysis/11_DESeq_loop_SCAIP7-18_W-A/")

# assign variabeles: needs modifications
variables <- c("pedu", "pincme", "psesl", "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasgm", "csinhl", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq" 
                )
 
variable_names <- c("Parental Education", "Parental Income","Subjective SES", "Neighborhood Stress", "Self-disclosure", 
                "Perceived responsiveness", "Parental Warmth", "YR Parent Child Conflict",  "Risky Family", "Youth Depression", 
                "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Stimulated IL-5 cortisol trt", "Stimulated IL-13 cortisol trt", "Stimulated IFNy cortisol trt",    
                "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", 
                "DD Asthma symptoms", "Nightly Asthma", "Nightly inhaler use", "Asthma Frequency", "Asthma Severity", 
                "Nightly Awakenings",  "Sleep Duration", "Sleep quality"
                )


#variables <- c("IFNG_co", "pnsi")
#variable_names <- c("Stimulated IFNy", "Neighborhood Stress")

############## assing treatment cell type, wave (if necessary before the loop), and filter covariate if necessary ####################
treatment <- "CTRL"
cell_types <- c("NKcell","Tcell","Bcell","Monocyte")
#Wave <-"a"
#covariates <- covariates %>% filter(!LibBatch %in% c('SCAIP2', 'SCAIP3'))
#length(unique(covariates$sample.ID))

#platePrefix <- "SCAIP7-18_W-A_DESeq_new_cov_01_03_24"

sink(paste0(platePrefix,"_stats_smmary_out.txt"))

timestamp()

tbl <- data.frame(matrix(ncol=0, nrow=2))

for(i in 1:length(cell_types)){
    
    cell <- cell_types[i]
            
    cat(" \n")
    cat("######################################################################################################\n")
    cat("############################### Processing Cell Type: ", cell, "#######################################","\n")

   
    cat("Now selecting count data for cell type: ", cell, "!\n")

   
    table <- data.frame(matrix(ncol=4, nrow=0))
    x <- c("symb", "variable", paste(cell, "_ngenes", sep=''), paste(cell, "_DEGs", sep=''))
    colnames(table) <- x


    # DESeq loop for each variable
        for(i in 1:length(variables)){
            myvar <- variables[i]
            var_name <- variable_names[i]
            
            cat("########################################################\n")
            cat("##Processing: ", var_name, "\n")

            topDirectory <- 'out_data_'
            outDir <- paste(topDirectory, platePrefix,"/",myvar, sep='')
            ##
            statsDir <- paste(outDir, '/stats', sep='')
            ##
            dataDir <- paste(outDir, '/data_objects', sep='')

            cat("Now testing ", var_name, " for cell type: ", cell, "!\n")
            cat(" \n")


          #load result
          fname=paste(statsDir, '/', myvar, "_", cell, "_", platePrefix, "_meta_analysis_results", ".txt", sep="")
          resVar <- read.table(fname, header=T)

          resVar_sig <- resVar %>% filter(qval <0.1)
          resVar_sig <- resVar_sig[order(resVar_sig$qval),]

          #cat("Number of individuals in variable", var_name , " in cell type ", cell," :", length(unique(ddsVar@colData$dbgap.ID)), "\n")
          cat("# of genes tested for - ", myvar, " in cell type ", cell," :", nrow(resVar)," \n")
          cat(" \n")

          cat("# of DEGs for - ", myvar,  " in celltype ", cell," (FDR < 0.1) : ", nrow(resVar_sig),"\n")

            tableVar <- data.frame(matrix(ncol=4, nrow=1))
            x <- c("symb", "variable", paste(cell, "_ngenes", sep=''), paste(cell, "_DEGs", sep='')) 
            colnames(tableVar) <- x
            tableVar[1] <- paste(myvar)
            tableVar[2] <- paste(var_name)

            #tableVar[3] <- paste(length(unique(ddsVar@colData$dbgap.ID)))
            tableVar[3] <- paste(nrow(resVar))
            tableVar[4] <- paste(nrow(resVar_sig))


           table <- rbind(table, tableVar)
        }

        # write table for the cell type
        tbl <- cbind(tbl, table)

    }         



    fname=paste("stats_table_all_variables_per_cell_type_genes_degs", "_", platePrefix,".txt", sep="")
    write.table(tbl, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

