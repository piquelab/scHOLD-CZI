library(DESeq2)
library(qvalue)
library(annotables)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(parallel)
library(ggseurat)
library(sva)
library(flextable)
library(ftExtra)
library(rlist)
library(officer)
library(zingeR)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
args <- c("ALOFT","/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_uniq-n265_psesl-a2_fixed_12-19-2024.txt","ALL","demux",0.2, 50, "/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt") 
args <- c("CZI","/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_ISEL_Mean_03_26_2024.txt","ALL","fastdemux",0.2,11) #for testing

job <- args[1]
base <- args[2]
cov_file=fread(args[3]) #this is the psych cov file
project=args[4]
method=args[5]
resset <- args[6]
dimset=args[7]
sample_batch <- args[8]
zingeRrun=FALSE

if(job=="ALOFT"){
outFolder=paste0(base,method,"_pseudobulk_ctrl/")
#for resolution 0.2 and including V2 chem
outFolder=paste0(outFolder,"lessfilt/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

exp <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt")
eigenvec2 <- merge(exp,cov_file,by="dbgap.ID",all.x=T)
#Ali shared that not all waves were run for singlecell and to make sure everything is grabbed correctly this filter is needed
eigenvec2 <- eigenvec2 %>% dplyr::filter(SCAIP7_18 == 1 | SCAIP1_6 == 1)
eigenvec2$Sample_ID <- eigenvec2$dbgap.ID
psychvarstorun <- colnames(eigenvec2[,24:59])
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

opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.RData")
load(opfn)

combatrun="income_PCs_sex_age_and_treats_adjusted"
run="income_PCs_sex_age_and_treats_adjusted_withWave_ageint"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
} else if(job=="CZI"){
    outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)
#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o[,-c("sex","sex_alph","age")],cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip","SNI_NoP","age","sex","sex_alph","isel") #SNI_NoP is the only variable that should be excluded based on the observed issues in score distributions that don’t make sense (negative values and extreme outliers)
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
focus_vars <-c("factor_HS_CRP","HS_CRP","PSQI_total","HVS_mean","EDS_mean","SES","LogCRP","cytocomp","SLS","PSS_all_mean","DED_all_mean","Chol_HDL","Chol_LDL","BPs_avg","BPd_avg")
# Variables to focus on
variablesL <- c("chronic_sum","cytocomp", "SES","PSS_all_mean", "Chol_HDL", "BPd_avg", "ISEL_Mean", "pr_comp") 
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
variables_df2 <- subset(variables_df, !variable %in% c("SNI_NoP"))
additionalvars <- fread("/rs/rs_grp_schold/covariates/other_covariates/chronic_chond_lables_attributes_varnames.txt")
variables_df <- merge(variables_df2,additionalvars[!additionalvars$Variable %in% c("HS_CRP"),-1],by.x=c("variable","description"),by.y=c("Variable","Column_Label"),all=T)
cat_cov <- fread(paste0(base,"categories_cov.txt"))
cat_cov_m <- merge(variables_df,cat_cov[,-2],by.y="Column_Name",by.x="variable",all.y=T)

opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.RData")
load(opfn)
combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
if(zingeRrun){
    run="SES_PCs_sex_age_and_treats_generem_zingeR_ageint"
}else{
    run="SES_PCs_sex_age_and_treats_generem_ageint"
}
}

combatfolder=paste0(baseoutFolder,combatrun,"/")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleidwave/"))) dir.create(paste0(outFolder,"sampleidwave/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid_all/"))) dir.create(paste0(outFolder,"sampleid_all/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

if(job=="ALOFT"){
withWave=TRUE
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
allvars <- psychvarstorun[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]
firstrunvars <- allvars[c(1:10)]
secondrunvars <- c(allvars[c(1:20)])
} else if(job=="CZI"){
psychvarstorun <- psychvarstorun[!psychvarstorun %in% c("CVDRISK")]
firstrunvars=c("SES","pr_comp","ISEL_Mean","PSS_all_mean","BPd_avg")
allvars=c(psychvarstorun[-c(45:66,72:94)])
secondrunvars=c(psychvarstorun[c(1:13,24:44)])
}
#firstrunvars <- "psesl"
#testing
#library("BiocParallel")
# register(MulticoreParam(8))

for (cluster in names(counts_ls)){
    #cluster="C0"
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    #cluster_metadata <- subset(cluster_metadata, !treats=="LPS-DEX") #removing due to low ind counts
    #opfn <- paste0(combatfolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    #run covariates separately for each treatment condition
    if(job=="ALOFT"){
    cluster_metadata <- transform(cluster_metadata, Sex=as.factor(Sex))
    if(all(c("0","1") %in% levels(cluster_metadata$Sex))){
        cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "0"))
    }
    } else if(job=="CZI"){
    adjusted <- adjusted_counts
    cluster_metadata <- transform(cluster_metadata, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
    cluster_metadata <- within(cluster_metadata, sex_alph <- relevel(sex_alph, ref = "Male"))
    }
    #for(i in unique(cluster_metadata_sce$treats)){
    for(i in c("RNA-CTRL","RNA-LPS")){
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        lapply(secondrunvars,function(var){
        #var="pnsi" var="DSES_07"
            if(var=="factor_HS_CRP"){
            cluster_metadata_t <- subset(cluster_metadata_t, HS_CRP<10) #advised to remove as likely an infection
            cluster_metadata_t <- transform(cluster_metadata_t, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
            cluster_metadata_t <- within(cluster_metadata_t, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
            }
            if (job=="ALOFT"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1",var)]
            resultscontrast=paste0(var,".cage1")
            } else if(job=="CZI"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
            resultscontrast=paste0(var,".age")
            }
            contrast=paste0(var,":age")
                #if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",i,".",run,".txt")) > 0)){
                cat("running deseq",var,cluster,i,contrast," \n")
                cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
                cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
                cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
                all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
                genestoremove <- data.frame(genes=unique(rownames(which(cluster_counts_t >= .Machine$integer.max, arr.ind = TRUE))))
                cluster_counts_t <- cluster_counts_t[!rownames(cluster_counts_t) %in% genestoremove$genes,]

                if(job=="ALOFT"){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + ", var, " + cage1 + ", var, ":cage1")
                } else if (job=="CZI"){
                design <-  paste0("~ PC1 + PC2 + sex_alph + ", var, " + age + " ,var, ":age")
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
                opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",contrast,".",cluster,".",i,".",run,".RDS")
                saveRDS(dds, file=opfn)

                res <- results(dds, name =resultscontrast)
                sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
                names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
                sub.table <- sub.table[!is.na(sub.table$padj), ]
                cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
                sub.table$var=var
                sub.table$cluster=cluster
                sub.table$treat=i
                sub.table$contrast=contrast
                fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",contrast,".",i,".",run,".txt"))
                sigDEGs <- subset(sub.table,padj<fdr)
                sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
                fwrite(sigDEGs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",contrast,".",i,".",run,".txt"))
                table <- data.frame(symb=var, description=variables_df[variables_df$variable==var,]$description, cluster=cluster,treat=i,contrast=contrast)
                table$number_samples <- paste(nrow(cluster_metadata_var))
                table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
                table$gene_number <- paste(nrow(cluster_counts_t))
                table$tested_genes <- paste(nrow(sub.table))
                table$DEGs_FDR <- paste(nrow(sigDEGs))
                table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
                fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",i,".",run,".txt"))
                #}
        })
    }
}


#for (var in c(psychvarstorun,"cage1")[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]){
for (var in secondrunvars){
    for (treat in c("RNA-CTRL","RNA-LPS")){
cat("running",var,treat,"\n")
myDir <- paste0(outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory 
filenames <- filenames[grep(paste0("-",var,":age.",treat),filenames)]
data_names <- gsub(paste0(".",treat,".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("stats_all_cell_types-", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(stats, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".stats_all_cell_types-",var,".",run,".",treat,".ageinteraction.txt"))

myDir <- paste0(outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0("_",var,":age.",treat),filenames)]
data_names <- gsub(paste0(".",treat,".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".deseqres_",var,".",run,".",treat,".ageinteraction.txt"))

myDir <- paste0(outFolder,"sigDEGs/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0("_",var,":age.",treat),filenames)]
data_names <- gsub(paste0(".",treat,".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("sigDEGs_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
df <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(df, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".sigDEGs_",var,".",run,".",treat,".ageinteraction.txt"))
}
}

for (treat in c("RNA-CTRL","RNA-LPS")){
myDir <- paste0(outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".ageinteraction.txt",filenames)]
filenames <- filenames[grep(paste0(treat),filenames)]
data_names <- gsub(paste0(".",run,".",treat,".ageinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".deseqres_"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#colnames(deseqres) <- c("identifier","padj","pvalue","logFC","SE","var","cluster","contrast")

myDir <- paste0(outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".ageinteraction.txt",filenames)]
filenames <- filenames[grep(paste0(treat),filenames)]
data_names <- gsub(paste0(".",run,".",treat,".ageinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".stats_all_cell_types-"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does

var50 <- subset(stats, DEGs_FDR_10>50)
stats50 <- subset(stats, contrast %in% var50$contrast)
wstats50 <- reshape(stats50[,-c(6,8,10)], idvar = c("symb","description","contrast"), timevar = "cluster", v.names=c("number_individuals","tested_genes","DEGs_FDR_10"), direction = "wide",sep=":")

 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

degcols <- grep("DEGs_FDR_10",colnames(wstats50))
stats50degcols <- degcols[ apply(wstats50[,degcols],MARGIN=2,FUN=my.max)<50]
#subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit( colnames(wstats50)[stats50degcols],"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(wstats50))] #old version
subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit( colnames(wstats50)[stats50degcols],"[.]"),function(y) y[1]),"NA"),collapse="|"),colnames(wstats50))]
dft <- subsubvars %>% flextable() %>% span_header(sep=":")
bodycol=ncol(subsubvars)-3
cols <- seq(4,bodycol,by=3)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1:2, j = c(4:length(subsubvars)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(cols), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".ageinteraction.",run,".",treat,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("number_individuals|tested_genes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() %>% split_header(sep=":")
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = c(4), border = border, part = "all")
dft <- align(dft, i = 2, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars_degonly)), align = "center", part = "body")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".ageinteraction.",run,".",treat,".ndegonly.png"))
}


    set_caption(caption = "Table 8.1")  




