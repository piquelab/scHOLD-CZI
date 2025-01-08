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

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_n521_uniq-n265_04-03-2024.txt","ALL","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt",0.2) 

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
run="income_PCs_sex_age_and_treats_adjusted_withWave"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
combatfolder=paste0(baseoutFolder,combatrun,"/")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleidwave/"))) dir.create(paste0(outFolder,"sampleidwave/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid_all/"))) dir.create(paste0(outFolder,"sampleid_all/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

runage=FALSE
withWave=TRUE
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
firstrunvars <- psychvarstorun[c(21:length(psychvarstorun))]
for (cluster in names(counts_ls)){
    #cluster="C0"
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))
    cluster_metadata <- subset(cluster_metadata, !treats=="LPS-DEX") #removing due to low ind counts
    #opfn <- paste0(combatfolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    #run covariates separately for each treatment condition
    cluster_metadata <- transform(cluster_metadata, Sex=as.factor(Sex))
    if(all(c("0","1") %in% levels(cluster_metadata$Sex))){
        cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "0"))
    }
    finalcontrast="cage1.treatsPHA.DEX_vs_cage1.treatsPHA"
    if(runage==TRUE){
        if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",finalcontrast,".",run,".txt")) > 0)){
            var="age"
            cat("running deseq",var,cluster," \n")
            contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA.DEX"))
            contrastdf <- transform(contrastdf, contrast1=paste0("cage1",".","treats",control),contrast2=paste0("cage1",".","treats",treatment))
            contrastdf <- transform(contrastdf,control=gsub("[.]","-",control),treatment=gsub("[.]","-",treatment))
            cluster_metadata_var <- cluster_metadata[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats")]

            mclapply(1:nrow(contrastdf),function(x) {
            con=contrastdf[x,]
            contrast=paste0(con$contrast2,"_vs_",con$contrast1)
            cat("running",contrast,"\n")
            cluster_metadata_var <- subset(cluster_metadata_var, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
            indcount <- plyr::count(cluster_metadata_var,"Sample_ID")
            cluster_metadata_var <- subset(cluster_metadata_var, Sample_ID %in% subset(indcount, freq>1)$Sample_ID)
            cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
            #fwrite(cluster_metadata_var[,c("Sample_ID","Wave")], sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sampleidwave_",var,"-",i,".",run,".txt"))
            cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
            cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
            all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
            if(withWave==TRUE){
            design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + cage1 + treats + treats:","cage1")
            } else {
            design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + cage1 + treats + treats:","cage1")
            }
            #may try this design:
            #design <-  paste0("~ Sample_ID + treats + treats:","cage1")
            cat("running deseq",var,cluster," \n")
            dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                          colData = cluster_metadata_var, 
                                          design = as.formula(design))
            dds <- DESeq(dds,parallel=TRUE)
            opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",contrast,".",cluster,".",run,".RDS")
            saveRDS(dds, file=opfn)

           #if(con$control=="CTRL"){
           #     res <- results(dds, name =con$contrast2)
           #     contrast=paste0(con$contrast2,"_vs_",con$contrast1)
           # } else{
           #     res <- results(dds, contrast = list(c(con$contrast2, con$contrast1)))
           #     contrast=paste0(con$contrast2,"_vs_",con$contrast1)
           # }
            res <- results(dds, name =con$contrast2)
            sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
            names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
            sub.table <- sub.table[!is.na(sub.table$padj), ]
            cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
            sub.table$var=var
            sub.table$cluster=cluster
            fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",contrast,".",run,".txt"))
            sigDEGs <- subset(sub.table,padj<fdr)
            sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
            fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",contrast,".",run,".txt"))
            table <- data.frame(symb="cage1", description=variables_df[variables_df$variable=="cage1",]$description, cluster=cluster)
            table$number_samples <- paste(nrow(cluster_metadata_var))
            table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
            table$gene_number <- paste(nrow(cluster_counts_t))
            table$tested_genes <- paste(nrow(sub.table))
            table$DEGs_FDR <- paste(nrow(sigDEGs))
            table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",run,".txt"))
            })
            }
            }
        lapply(firstrunvars[c(1:length(firstrunvars))][!firstrunvars[c(1:length(firstrunvars))] %in% puberty],function(var){
        #var="pnsi"
            finalcontrast=paste0(var,".treatsPHA.DEX_vs_",var,".treatsPHA")
            contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA.DEX"))
            contrastdf <- transform(contrastdf, contrast1=paste0(var,".","treats",control),contrast2=paste0(var,".","treats",treatment))
            contrastdf <- transform(contrastdf,control=gsub("[.]","-",control),treatment=gsub("[.]","-",treatment))
            cluster_metadata_var <- cluster_metadata[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]

            mclapply(1:nrow(contrastdf),function(x) {
                con=contrastdf[x,]
                contrast=paste0(con$contrast2,"_vs_",con$contrast1)
                if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",run,".txt")) > 0)){
                cat("running deseq",var,cluster,contrast," \n")
                cluster_metadata_var <- subset(cluster_metadata_var, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
                cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
                fwrite(cluster_metadata_var[,c("Sample_ID","Wave")], sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sampleidwave_",contrast,".",run,".txt"))
                cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
                cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
                all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

                if(withWave==TRUE){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + cage1 + ", var, " + treats + treats:",var)
                } else {
                design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + cage1 + ", var, " + treats + treats:",var)
                }
                cat("running deseq",var,cluster," \n")
                dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                              colData = cluster_metadata_var, 
                                              design = as.formula(design))
                dds <- DESeq(dds,parallel=TRUE)
                opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",contrast,".",cluster,".",run,".RDS")
                saveRDS(dds, file=opfn)

                res <- results(dds, name =con$contrast2)
                sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
                names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
                sub.table <- sub.table[!is.na(sub.table$padj), ]
                cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
                sub.table$var=var
                sub.table$cluster=cluster
                sub.table$contrast=contrast
                fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",contrast,".",run,".txt"))
                sigDEGs <- subset(sub.table,padj<fdr)
                sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
                fwrite(sigDEGs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",contrast,".",run,".txt"))
                table <- data.frame(symb=var, description=variables_df[variables_df$variable==var,]$description, cluster=cluster,contrast=contrast)
                table$number_samples <- paste(nrow(cluster_metadata_var))
                table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
                table$gene_number <- paste(nrow(cluster_counts_t))
                table$tested_genes <- paste(nrow(sub.table))
                table$DEGs_FDR <- paste(nrow(sigDEGs))
                table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
                fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",run,".txt"))
                }
                })
        })
}



for (var in c(psychvarstorun,"cage1")[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]){
cat("running",var,"\n")
myDir <- paste0(outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0(var,".treats"),filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.", "", data_names)
data_names <- gsub("stats_all_cell_types-", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
stats$contrast <- gsub(paste0(var,".treats"),"",stats$contrast)
fwrite(stats, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".stats_all_cell_types-",var,".",run,".treatinteraction.txt"))

myDir <- paste0(outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0(var,".treats"),filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.", "", data_names)
data_names <- gsub("deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
deseqres$contrast <- gsub(paste0(var,".treats"),"",deseqres$contrast)
fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".deseqres_",var,".",run,".treatinteraction.txt"))

myDir <- paste0(outFolder,"sigDEGs/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0(var,".treats"),filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.", "", data_names)
data_names <- gsub("sigDEGs_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
df <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
df$contrast <- gsub(paste0(var,".treats"),"",df$contrast)
fwrite(df, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".sigDEGs_",var,".",run,".treatinteraction.txt"))
}

#Work in progress

myDir <- paste0(outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".treatinteraction.txt",filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.treatinteraction.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#colnames(deseqres) <- c("identifier","padj","pvalue","logFC","SE","var","cluster","contrast")

myDir <- paste0(outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".treatinteraction.txt",filenames)]
data_names <- gsub(".income_PCs_sex_age_and_treats_adjusted_withWave.treatinteraction.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.50.stats_all_cell_types-", "", data_names)
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
subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit(colnames(wstats50[,degcols[ apply(wstats50[,degcols],MARGIN=2,FUN=my.max)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(wstats50))]
dft <- subsubvars[,-4] %>% flextable() %>% span_header(sep=":")
cols <- seq(3,36,by=3)
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




