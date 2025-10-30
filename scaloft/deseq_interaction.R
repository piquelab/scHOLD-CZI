library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(plyr)
library(parallel)
library(ggseurat)
library(sva)
library(data.table)
library(flextable)
library(ftExtra)
library(rlist)
library(officer)
library(zingeR)

#future::plan(strategy = 'multicore', workers = 10)
#options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
#job="ALOFT"
job="CZI"

if(job=="ALOFT"){
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_uniq-n265_psesl-a2_fixed_12-19-2024.txt","ALL","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt",0.2) 
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
sample_batch <- args[5]
resset <- args[6]
dimset=50
outFolder=paste0(base,method,"_pseudobulk_ctrl/")
#for resolution 0.2 and including V2 chem
outFolder=paste0(outFolder,"lessfilt/")
figuredir=paste0(outFolder,"figures/")
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
exp <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt")
eigenvec2 <- merge(exp,cov_file,by="dbgap.ID",all.x=T)
eigenvec2 <- eigenvec2 %>% dplyr::filter(SCAIP7_18 == 1 | SCAIP1_6 == 1)
eigenvec2$Sample_ID <- eigenvec2$dbgap.ID
eigenvec2 <- transform(eigenvec2, SCAIP1_6=ifelse(is.na(SCAIP1_6),0,SCAIP1_6)) #weird case where some are NA ... Ali not sure why
psychvarstorun <- colnames(eigenvec2[,24:59])
treatments=c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")
treatmentsfirst=c("CTRL","PHA")
contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA-DEX")) #no lps vs lps-dex as too few ind
opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)
combatrun="income_PCs_sex_age_and_treats_adjusted"

puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
allvars <- psychvarstorun[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]
firstrunvars <- allvars[c(1:10)]
secondrunvars <- c(allvars[c(1:20)],"cage1")
} else if(job=="CZI"){
proteincoding=FALSE
withCOMBAT=FALSE
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.15) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
contrastdf <- data.frame(control=c("RNA-CTRL"),treatment=c("RNA-LPS"))
#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
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
if(proteincoding){
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt_proteincoding/")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt_proteincoding.RData")
#load(opfn)
} else{
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
load(opfn)
}
#counts_ls$C6 <- NULL
#metadata_ls$C6 <- NULL
combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
}
contrastdf <- transform(contrastdf, contrast=paste0(treatment,"_vs_",control))

if(job=="ALOFT"){
variables <- c("pedu", "pincme", 
                "psesl", 
                "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasg", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq", 
                "genPC1", "genPC2", "genPC3",  "Sex", "cage1", "ceth1", "cwght1", "chght1", "csex1",
                "cgpd5", "cgpd", "cbpd"
                )
 
variable_names <- c("Parental Education", "Parental Income",
                "Subjective SES", 
                "Neighborhood Stress", "Self-disclosure", 
                "Perceived responsiveness", "Parental Warmth", "YR Parent Child Conflict",  "Risky Family", "Youth Depression", 
                "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Stimulated IL-5 cortisol trt", "Stimulated IL-13 cortisol trt", "Stimulated IFNy cortisol trt",    
                "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", 
                "DD Asthma symptoms", "Nightly Asthma", "Asthma Frequency", "Asthma Severity", 
                "Nightly Awakenings",  "Sleep Duration", "Sleep quality",
                "New Genotype PC1", "New Genotype PC2", "New Genotype PC3","Sex_alph", "Age", "Ethnicity", "Weight", "Height", "sex",
                "Female Menarche Status", "Female Puberty Score", "Male puberty score"
                )
variables_df <- data.frame(variable=variables, description=variable_names)
allvarsdf <- data.frame(variable=allvars)
variables_df <- unique(merge(variables_df,allvarsdf,by="variable",all=T))

} else if(job=="CZI"){
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
}

if(withCOMBAT){
    baserun="treatvarint_withCOMBAT_deseq"
} else {
    baserun="treatvarint_noCOMBAT_deseq"
}

runPFAS=FALSE
zingeRrun=FALSE
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

#testing - works if have to turn off parallel for deseq for issues
#library("BiocParallel")
# register(MulticoreParam(8))

#for CZI want to try binary high vs low stress/isel
if(job=="CZI"){
metadata_ls <- lapply(metadata_ls,function(m){
    #cluster="C0"
    cluster_metadata <- data.frame(m)
    cat("running ", unique(m$letter_clusters), "\n")
    cluster_metadata_split <- cluster_metadata %>%
      mutate(
        # Create the quantile-based categories
        PSS_bin = cut(PSS_all_mean,
                       breaks = quantile(PSS_all_mean, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                       labels = c("Low", "Medium", "High"),
                       include.lowest = TRUE),
        ISEL_bin = cut(ISEL_Mean,
                       breaks = quantile(ISEL_Mean, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                       labels = c("Low", "Medium", "High"),
                       include.lowest = TRUE)
      ) 
    return(cluster_metadata_split)
})
}


for (cluster in names(counts_ls)){
    #cluster="C0"
        cat("running ", cluster, "\n")
    fdr=0.1
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    #cluster_metadata <- subset(cluster_metadata, !treats=="LPS-DEX") #removing due to low ind counts
    if(withCOMBAT){
        opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
        load(opfn)
    } else{
        cluster_counts_sce <- counts_ls[[cluster]]
        cluster_counts <- assay(cluster_counts_sce, "counts")
        adjusted_counts <- cluster_counts   
    }

    #run covariates separately for each treatment condition
    if(job=="ALOFT"){
        adjusted_counts <- adjusted
            cluster_metadata <- transform(cluster_metadata, Sex=as.factor(Sex),treats=as.factor(treats))
        if(all(c("0","1") %in% levels(cluster_metadata$Sex))){
            cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "0"))
        }
    } else if(job=="CZI"){
        cluster_metadata <- transform(cluster_metadata, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r),treats=as.factor(treats))
        cluster_metadata <- within(cluster_metadata, sex_alph <- relevel(sex_alph, ref = "Male"))
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
        #keep <- rowSums(adjusted_countsP,na.rm=T) >= (ncol(adjusted_countsP) / 2) #filter to keep genes expressed in 25% of samples
        #adjusted_counts <- unlist(adjusted_countsP[keep, ])
        }
    }
        lapply(c("PSS_bin","ISEL_bin"),function(var){
        #var="pnsi" var="DSES_07" var="PSS_bin"
            if(var=="factor_HS_CRP"){
            cluster_metadata <- subset(cluster_metadata, HS_CRP<10) #advised to remove as likely an infection
            cluster_metadata <- transform(cluster_metadata, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
            cluster_metadata <- within(cluster_metadata, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
            }
            if(var=="PSS_bin"){
                  # Filter to keep only the "Low" and "High" groups
                cluster_metadata <- subset(cluster_metadata, PSS_bin != "Medium") #advised to remove as likely an infection
                cluster_metadata <- within(cluster_metadata, PSS_bin <- relevel(PSS_bin, ref = "Low"))
                cluster_metadata$PSS_bin <- factor(cluster_metadata$PSS_bin)
            }
            if(var=="ISEL_bin"){
                  # Filter to keep only the "Low" and "High" groups
                cluster_metadata <- subset(cluster_metadata, ISEL_bin != "Medium") #advised to remove as likely an infection
                cluster_metadata <- within(cluster_metadata, ISEL_bin <- relevel(ISEL_bin, ref = "Low"))
                cluster_metadata$ISEL_bin <- factor(cluster_metadata$ISEL_bin)
            }

            lapply(1:nrow(contrastdf),function(x) {
            #for (x in 1:nrow(contrastdf)){
            con=contrastdf[x,]
            contrast=paste0(con$treatment,"_vs_",con$control)
                if(!isTRUE(file.size(paste0(outFolder,"stats/",cluster,".",contrast,".",var,".stats_all_cell_types",".",run,".txt")) > 0)){
                cat("running deseq",var,cluster,contrast," \n")
                cluster_metadata_t <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
                cluster_metadata_t <- within(cluster_metadata_t, treats <- factor(relevel(treats, ref = con$control)))
                if(job=="ALOFT"){
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]
                } else if (job=="CZI"){
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age","treats",var)]
                }
                if(var=="age" | var=="cage1"){
                if(job=="ALOFT"){
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","treats",var)]
                } else if (job=="CZI"){
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","treats",var)]
                }}
                if(var=="CVDRISK"){
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","BATCH","PC1","PC2","treats",var)]
                }
                cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
                cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
                cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
                all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
                genestoremove <- data.frame(genes=unique(rownames(which(cluster_counts_t >= .Machine$integer.max, arr.ind = TRUE))))
                cluster_counts_t <- cluster_counts_t[!rownames(cluster_counts_t) %in% genestoremove$genes,]

                if(job=="ALOFT" & var=="cage1"){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + ", var, " + treats + treats:",var)
                } else if (job=="ALOFT" & !var=="cage1"){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + cage1 + ", var, " + treats + treats:",var)
                } else if (job=="CZI" & var=="age"){
                design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + ", var, " + treats + treats:",var)
                } else if(job=="CZI" & var=="CVDRISK"){
                design <-  paste0("~ BATCH + PC1 + PC2 + ",var, " + treats + treats:",var)
                } else if (job=="CZI" & !var=="age"){
                design <-  paste0("~ BATCH + PC1 + PC2 + sex_alph + age + ", var, " + treats + treats:",var)
                }
                dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                              colData = cluster_metadata_var, 
                                              design = as.formula(design))
                if(zingeRrun){
                designb = model.matrix( as.formula(design) , cluster_metadata_var)
                weights <- zeroWeightsLS(counts=cluster_counts_t, design=designb, maxit=200, normalization="DESeq2_poscounts", colData=cluster_metadata_var, designFormula=as.formula(design))
                assays(dds)[["weights"]]=weights
                dds = DESeq2::estimateSizeFactors(dds, type="poscounts")
                dds = estimateDispersions(dds)
                dds = nbinomWaldTest(dds, betaPrior=FALSE, useT=TRUE, df=rowSums(weights)-2) #betaPrior=FALSE should be used for designs with interactions >> note: orignial script is betaPrior=TRUE
                } else {
                dds <- DESeq(dds,parallel=TRUE) #was parallel=TRUE causing issue?
                }
                opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",contrast,".",cluster,".",run,".RDS")
                saveRDS(dds, file=opfn)
    
                con <- transform(con, control=gsub("-",".",con$control),treatment=gsub("-",".",con$treatment),contrast=gsub("-",".",con$contrast))
                res <- results(dds, name =paste0(var,".treats",con$treatment))
                sub.table <- data.frame(res@rownames, res$'baseMean', res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
                names(sub.table) <- c('identifier','baseMean','padj','pvalue', 'logFC','SE')
                sub.table <- sub.table[!is.na(sub.table$padj), ]
                cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
                sub.table$var=var
                sub.table$cluster=cluster
                sub.table$contrast=contrast
                sigDEGs <- subset(sub.table,padj<0.1)
                varmarginal_res <- results(dds, name =var)
                varmarginal_resdf <- data.frame(varmarginal_res@rownames, varmarginal_res$'baseMean', varmarginal_res$'padj', varmarginal_res$'pvalue', varmarginal_res$'log2FoldChange', varmarginal_res$'lfcSE',stringsAsFactors=FALSE)
                names(varmarginal_resdf) <- c('identifier','baseMean','padj','pvalue', 'logFC','SE')
                treatmarginal_res <- results(dds, name =paste0("treats_",con$contrast))
                treatmarginal_resdf <- data.frame(treatmarginal_res@rownames, treatmarginal_res$'baseMean', treatmarginal_res$'padj', treatmarginal_res$'pvalue', treatmarginal_res$'log2FoldChange', treatmarginal_res$'lfcSE',stringsAsFactors=FALSE)
                names(treatmarginal_resdf) <- c('identifier','baseMean','padj','pvalue', 'logFC','SE')
                table <- data.frame(symb=var, description=variables_df[variables_df$variable==var,]$description, cluster=cluster,contrast=contrast)
                table$number_samples <- paste(nrow(cluster_metadata_var))
                table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
                table$gene_number <- paste(nrow(cluster_counts_t))
                table$tested_genes <- paste(nrow(sub.table))
                table$intDEGs <- paste(nrow(subset(sub.table,padj<0.1)))
                table$varDEGs <- paste(nrow(subset(varmarginal_resdf,padj<0.1)))
                table$treatDEGs <- paste(nrow(subset(treatmarginal_resdf,padj<0.1)))
                fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",cluster,".",contrast,".",var,".deseq_interaction.",run,".txt"))
                fwrite(treatmarginal_resdf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",cluster,".",contrast,".",var,".deseq_treat.",run,".txt"))
                fwrite(varmarginal_resdf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",cluster,".",contrast,".",var,".deseq_variable.",run,".txt"))
                fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",cluster,".",contrast,".",var,".stats_all_cell_types",".",run,".txt"))
                }
                })
        })
}



#for (var in c(psychvarstorun,"cage1")[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]){
for (var in updated_vars){
cat("running",var,"\n")
myDir <- paste0(outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0("-",var,".treats"),filenames)]
data_names <- gsub(paste0(".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("stats_all_cell_types-", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
stats$contrast <- gsub(paste0(var,".treats"),"",stats$contrast)
fwrite(stats, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".stats_all_cell_types-",var,".",run,".treatinteraction.txt"))

myDir <- paste0(outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0("_",var,".treats"),filenames)]
data_names <- gsub(paste0(".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
deseqres$contrast <- gsub(paste0(var,".treats"),"",deseqres$contrast)
fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".deseqres_",var,".",run,".treatinteraction.txt"))

myDir <- paste0(outFolder,"sigDEGs/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0("_",var,".treats"),filenames)]
data_names <- gsub(paste0(".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("sigDEGs_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
df <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
df$contrast <- gsub(paste0(var,".treats"),"",df$contrast)
fwrite(df, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".sigDEGs_",var,".",run,".treatinteraction.txt"))
}

myDir <- paste0(outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".treatinteraction.txt",filenames)]
data_names <- gsub(paste0(".",run,".treatinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".deseqres_"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#colnames(deseqres) <- c("identifier","padj","pvalue","logFC","SE","var","cluster","contrast")

myDir <- paste0(outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".treatinteraction.txt",filenames)]
data_names <- gsub(paste0(".",run,".treatinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".stats_all_cell_types-"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
colnames(stats) <- c("symb","description","cluster","contrast","number_samples","number_individuals","gene_number","tested_genes","DEGs_FDR","sigDEGs")
stats <- transform(stats, symb_contrast=paste(symb,contrast,sep=":"))

var50 <- subset(stats, sigDEGs>50)
stats50 <- subset(stats, symb_contrast %in% var50$symb_contrast)
wstats50 <- reshape(stats50[,-c(7,9,11)], idvar = c("symb","description","contrast"), timevar = "cluster", v.names=c("number_individuals","tested_genes","sigDEGs"), direction = "wide",sep=":")

 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

degcols <- grep("sigDEGs",colnames(wstats50))
stats50degcols <- degcols[ apply(wstats50[,degcols],MARGIN=2,FUN=my.max)<50]
#subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit( colnames(wstats50)[stats50degcols],"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(wstats50))] #old version
subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit( colnames(wstats50)[stats50degcols],"[.]"),function(y) y[1]),"NA"),collapse="|"),colnames(wstats50))]
dft <- subsubvars[,-4] %>% flextable() %>% span_header(sep=":")
bodycol=ncol(subsubvars)-3
cols <- seq(3,bodycol,by=3)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1:2, j = c(4:length(subsubvars[,-4])), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars[,-4])), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(cols), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".treatinteraction.",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("number_individuals|tested_genes",colnames(subsubvars))]
dft <- subsubvars_degonly[,-4] %>% flextable() %>% split_header(sep=":")
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = c(3), border = border, part = "all")
dft <- align(dft, i = 2, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars_degonly[,-4])), align = "center", part = "body")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".treatinteraction.",run,".ndegonly.png"))



    set_caption(caption = "Table 8.1")  




