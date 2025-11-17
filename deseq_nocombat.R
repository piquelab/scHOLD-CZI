require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
library(annotables)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(plyr); library(dplyr)
library(parallel)
library("AnnotationHub")
library(ggseurat)
library(cowplot)
library(sva)
library(ggpubr)
future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 15 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.15) #for testing
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.1) #for testing

#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",11,0.2) #for testing

base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
outFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
#outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)
leadvar <- fread("/rs/rs_grp_schold/covariates/other_covariates/HOLD LEAD 5.27.25.csv")
cov_pluslead <- merge(cov_file,leadvar,by="pID",all=T)
cov_pluslead <- transform(cov_pluslead, Lead=ifelse(Lead==-99,NA,Lead))
#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o[,-c("sex","sex_alph","age")],cov_pluslead,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip","SNI_NoP","age","sex","sex_alph","isel","pID") #SNI_NoP is the only variable that should be excluded based on the observed issues in score distributions that don’t make sense (negative values and extreme outliers)
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
old_cytokines <- psychvarstorun[c(8,14:24)] #old cytokines
psychvarstorun <- psychvarstorun[!psychvarstorun %in% old_cytokines]
firstrunvars=c("SES","pr_comp","ISEL_Mean","PSS_all_mean","BPd_avg","Chol_HDL","Chol_LDL","nii_mean","SNI_NumPeople_r","BPs_avg","DED_all_mean","LogCRP","HVS_mean","LivingAlone")
reordered_psychvarstorun <- psychvarstorun[c(which(psychvarstorun %in% firstrunvars),which(!psychvarstorun %in% firstrunvars))]
secondrunvars <- c(reordered_psychvarstorun[c(1:20)],"age","Lead")
zcytokines <- c("cytocomp","z_ifny_0_log_w",   
"z_il10_0_log_w","z_il12p70_0_log_w",
"z_il13_0_log_w",   
"z_il1b_0_log_w","z_il2_0_log_w",    
"z_il4_0_log_w",    
"z_il6_0_log_w","z_il8_0_log_w",    
"z_tnfa_0_log_w")
treatments=c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")
treatmentsfirst=c("RNA-CTRL","RNA-LPS")
allPFAS <- psychvarstorun[c(34:54)]
combPFAS <- allPFAS[seq(1,21,3)]
PFASvars <- c(allPFAS,"sumPFAS","log_sumPFAS")
allvars <- c(reordered_psychvarstorun[c(1:83)],"age","Lead","sumPFAS","log_sumPFAS")
varinterest <- c("ISEL_Mean","PSS_all_mean")#,"cytocomp")

# Variables to focus on
    variables <- c("DSES_07","DSES_09","FCDEM_12","chronic_sum","smoke","CVDRISK",
                "HS_CRP","il6","PSQI_total","HVS_mean","EDS_mean","SES","LogCRP","Logifny","Logil10","Logil12",
                "Logil13","Logil1b","Logil2","Logil4","Logil6","Logil8","Logtnfa","cytocomp","NAIscr2019","NAIscr2017",
                "LQ2019","StressSev","StressCount","nii_mean","NII_fam","SLS","LivingAlone","SNI_HCG","SNI_NoP","PSS_all_mean",
                "DED_all_mean","Chol","Trig","Chol_HDL","Chol_Ratio","Chol_LDL","HbA1C","BPs_avg","BPd_avg","ISEL_Mean", "pr_comp",
                "SNI_NumPeople_r","LogDED","LogTrig","chronic_sum_categ",
                "age","factor_HS_CRP"
                )
variable_names <- c("Participant's highest level of education","Pre-tax household income","Are you currently taking any prescription medications",
                    "Total number of chronic conditions","Smoking status","Cardiovascular disease risk",
                    "High Sensitivity CRP (mg/L) raw values","IL-6 raw values","Total sleep score (PSQI)","Vigilance Scale - Mean","Everyday Discrimination Scale - Mean",
                    "SES composite [mean(ZDSES_07, ZDSES_09)]","LogCRP","Logifny","Logil10","Logil12",
                    "Logil13","Logil1b","Logil2","Logil4","Logil6",
                    "Logil8","Logtnfa","Cytokines composite [mean(ZLogtnfa, ZLogi6, ZLogifny)]","Neighborhood Adversity Score","Neighborhood Adversity Score",
                    "Neighborhood Segregation Score","Core: Total Severity of Stressors","Core: Total Count of Stressors","Negative interactions mean","Negative interactions with family","Avg on loneliness scale",
                    "Dichotomous variable measuring living alone (1 = living alone)","High contact groups: number of categories of social groups p interacted with at least once every two weeks",
                    "This is the mean score representing the number of people with whom the respondent has regular contact","Mean of daily pss measures","Daily experience of discrimination cumulative average",
                    "Cholesterol in mg/dL","Triglycerides in mg/dL","HDL Cholesterol in mg/dL","CHOL/HDL ratio","LDL Cholesterol in mg/dL",
                    "Glycated HGB Affinity HPLC HbA1C","Average systolic blood pressure (wonky BP3_s measure corrected)","Average diastolic blood pressure","social support","psychological resources composite",
                    "number of pople in contact with","Log of Daily Discrimination", "LogTrig", "Chronic conditions transformed",
                    "age","High Sensitivity CRP (mg/L) factorized 1 to 5"
                    )
variables_df <- data.frame(variable=variables, description=variable_names)
additionalvars <- fread("/rs/rs_grp_schold/covariates/other_covariates/chronic_chond_lables_attributes_varnames.txt")
variables_df <- merge(variables_df,additionalvars[!additionalvars$Variable %in% c("HS_CRP"),-1],by.x=c("variable","description"),by.y=c("Variable","Column_Label"),all=T)
allvarsdf <- data.frame(variable=c(allvars,PFASvars,"Lead"))
variables_df <- merge(variables_df,allvarsdf,by="variable",all=T)

cat_cov <- fread(paste0(base,"categories_cov.txt"))
cat_cov_m <- merge(variables_df,cat_cov[,-2],by.y="Column_Name",by.x="variable",all.y=T)

## RPR AR:
## this input was changed to the input with the suffix: .DESeq_countlists.bticfilt_proteincoding.RData
## .DESeq_countlists.bticfilt.RData was generated by filtering genes that was used in Combat, which we are not currently using
#opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")

opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt_proteincoding.RData")
load(opfn)

baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
run="noCombat_DESeq"
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

runPFAS=FALSE
zingeRrun=FALSE
iselcov=FALSE
if(iselcov){
    run=paste0(run,"_iselcov")
}
SEScov=FALSE
if(SEScov){
    run=paste0(run,"_SEScov")
}
WHRcov=TRUE
if(WHRcov){
    run=paste0(run,"_WHRcov")
}
PSScov=FALSE
if(PSScov){
    run=paste0(run,"_PSScov")
}
if(zingeRrun){
    run=paste0(run,"_zingeR")
}

lapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    #cluster_metadataL <- left_join(cluster_metadata,unique(eigenvec2[,c("Sample_ID","Lead","WHR")]),by="Sample_ID")
    #rownames(cluster_metadataL) <- cluster_metadataL$rowid
    #cluster_metadata <- cluster_metadataL

    if(runPFAS){
        cluster_metadata_pfas <- cluster_metadata %>% 
        rowwise() %>% 
        mutate(
        sumPFAS = sum(!!!syms(combPFAS),na.rm=T)
        )
    cluster_metadata <- as.data.frame(cluster_metadata_pfas)
    cluster_metadata <- transform(cluster_metadata, log_sumPFAS=log2(sumPFAS+1))
    rownames(cluster_metadata) <- cluster_metadata$rowid
    #adjusted_countsP <- adjusted_counts[rownames(adjusted_counts) %in% genes$hgnc_symbol, ]
    #keep <- rowSums(adjusted_countsP,na.rm=T) >= (ncol(adjusted_countsP) / 2) #filter to keep genes expressed in 50% of samples
    #adjusted_counts <- unlist(adjusted_countsP[keep, ])
    }

    for(i in treatmentsfirst){
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

## AR: script only ran for zcytokines. changed this to the top three variables instead:  varinterest
        #lapply(zcytokines,function(var){
        lapply(varinterest,function(var){

            #c("ISEL_Mean","PSS_all_mean","cytocomp")
            #var="pr_comp"
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
                cat("running deseq ",cluster, var,i," \n")
        #if(var=="Lead"){
        #    adjusted_counts <- adjusted_counts[rownames(adjusted_counts) %in% genes$hgnc_symbol, ]
        #}

        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + age + ",var)
        if(var=="factor_HS_CRP"){
        cluster_metadata_t <- subset(cluster_metadata_t, HS_CRP<10) #advised to remove as likely an infection
        cluster_metadata_t <- transform(cluster_metadata_t, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
        cluster_metadata_t <- within(cluster_metadata_t, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
        }
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age",var)]
        if(SEScov){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + age + SES + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age","SES",var)]
        }
        if(iselcov){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + age + ISEL_Mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age","ISEL_Mean",var)]
        }
        if(WHRcov){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + age + WHR + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age","WHR",var)]
        }
        if(PSScov){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + age + PSS_all_mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age","PSS_all_mean",var)]
        }
        if(var=="CVDRISK"){
        design <-  paste0("~ BATCH + PC1 + PC2 + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2",var)]
        if(SEScov){
        design <-  paste0("~ BATCH + PC1 + PC2 + SES + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","SES",var)]
        }
        if(WHRcov){
        design <-  paste0("~ BATCH + PC1 + PC2 + WHR + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","WHR",var)]
        }
        if(iselcov){
        design <-  paste0("~ BATCH + PC1 + PC2 + ISEL_Mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","ISEL_Mean",var)]
        }
        } #end of cvdrisk
        if(var=="age"){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph",var)]
        if(SEScov){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + SES + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","SES",var)]
        }
        if(WHRcov){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + WHR + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","WHR",var)]
        }
        if(iselcov){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + ISEL_Mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","ISEL_Mean",var)]
        }
        if(PSScov){
        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + PSS_all_mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","PSS_all_mean",var)]
        }
        } # end of age
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- try(DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        )
        if(inherits(dds, "try-error")){
          print("An error occurred")
        } else {
        if(zingeRrun){
                    designb = model.matrix( as.formula(design) , cluster_metadata_var)
        weights <- zeroWeightsLS(counts=cluster_counts_t, design=designb, maxit=200, normalization="DESeq2_poscounts", colData=cluster_metadata_var, designFormula=as.formula(design))
        assays(dds)[["weights"]]=weights
        dds = DESeq2::estimateSizeFactors(dds, type="poscounts")
        dds = estimateDispersions(dds)
        dds = nbinomWaldTest(dds, betaPrior=FALSE, useT=TRUE, df=rowSums(weights)-2)
        } else {
        dds <- DESeq(dds,parallel=TRUE)
        }

        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
        save(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'baseMean', res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier','baseMean','padj','pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$tested_genes <- paste(nrow(sub.table))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
        rm(table,sigDEGs,sub.table,res,dds,cluster_metadata_var)
        gc()
            }
            }
        })
    }
})

## AR: script only ran for zcytokines. changed this to the top three variables instead: varinterest
#for (var in c(zcytokines)){
 for (var in c(varinterest)){
    #for (var in c("PSS_all_mean","ISEL_Mean","cytocomp")){
    for (i in c(treatmentsfirst)){
    if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    }
    }
}

#no independent filtering
for (i in c(treatmentsfirst)){
for(run in c("noCombat_DESeq","noCombat_DESeq_PSScov","noCombat_DESeq_iselcov")){
    for (var in c("ISEL_Mean","PSS_all_mean")){
    lapply(names(counts_ls),function(cluster) {
  cat("running",i, run, var, cluster,"\n")
    opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
    load(opfn)
    res <- results(dds)
    deseqres <- data.frame(res@rownames, res$'baseMean', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
    names(deseqres) <- c('identifier', 'baseMean', 'pvalue', 'logFC','SE')
    deseqres <- transform(deseqres, padj=p.adjust(pvalue,method="BH"),var=var,cluster=cluster,treats=i)
    fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".noindfilt.txt"))
    })
    }
}
}

#################################################### 
#LRT test 
#Perform likelihood ratio tests by incorporating both psychological stress and social support variables simultaneously into the statistical models. 
#This approach allows us to rigorously assess the independent effects of each psychosocial factor on gene expression patterns, while controlling for potential confounding influences.
#also: model with both, and we contrast to the model without any of the two

lapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    for(i in treatmentsfirst){
            #i <- treatmentsfirst[1]
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))
        
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age","PSS_all_mean","ISEL_Mean")]
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + age + PSS_all_mean + ISEL_Mean")
        reduced_design <- paste0("~ BATCH + PC1 + PC2 + sex_alph + age")
        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))

        dds <- DESeq(dds, test = "LRT", reduced = as.formula(reduced_design))

        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output_LRT-stress_socialsupport-",i,"-",cluster,".",run,".RData")
        save(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'baseMean', res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier','baseMean','padj','pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var="PSS_isel_LRT"
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_LRT-stress_socialsupport-",i,".",run,".txt"))
        table <- data.frame(variable= "PSS_isel_LRT", cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$tested_genes <- paste(nrow(sub.table))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types_LRT-stress_socialsupport-",i,".",run,".txt"))
        rm(table,sigDEGs,sub.table,res,dds,cluster_metadata_var)
        gc()
    }
})

for (var in c("LRT-stress_socialsupport")){
    #for (var in c("PSS_all_mean","ISEL_Mean","cytocomp")){
    for (i in c(treatmentsfirst)){
    if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types_",var,"-",i,".",run,".txt")) > 0)){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types_",var,"-",i,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    }
    }
}