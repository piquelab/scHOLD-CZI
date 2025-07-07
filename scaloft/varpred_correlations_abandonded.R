library(data.table)
library(pheatmap)
library(ggplot2)
library(tidyverse)
future::plan(strategy = 'multicore', workers = 10) #had an issue: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers
options(future.globals.maxSize = 30 * 1024 ^ 3)

treat="CTRL"
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
resmethod="voom"
combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
glmnetfolder=paste0(baseoutFolder,"glmnet/")
outFolder=paste0(glmnetfolder,resmethod,"/")
figuredir=paste0(outFolder,"figures/")


#get the predictable trait and merge (once and for all):
#big table with first column as sample and others as variables

myDir <- outFolder
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep("pred.txt",filenames)]
filenames <- filenames[grep(paste0(treat),filenames)]
data_names <- gsub("pred.txt", "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub(paste0(treat,"."), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], setnames(fread(file.path(myDir, filenames[i]),header = FALSE, sep='\t'),c("Sample_ID",paste0(data_names[i])))) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
df <- as.data.frame(Reduce(function(x, y) merge(x, y,by="Sample_ID",all=T),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(df, paste0(outFolder,"all-predictions.txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

samples <- as.character(df[,1])
rownames(preds) <- as.character(df[,1])
preds <- df[,-1]

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

preds1 <- preds[,!colnames(preds) %in% c("csex1")]

corr <- cor(preds1, use="pairwise.complete.obs")
library("psych")
corr.psych <- corr.test(preds1, adjust="none")
# blank out non-significant correlations:
corr <- corr*(corr.psych$p<0.05)
pairs <- reshape2::melt(corr)
pairs <- pairs[!pairs$value %in% c(0,1),]
write.table(pairs[,1:2], paste0(outFolder,"correlated-metagenes-list.txt"), sep="\t", quote=F, col.names=F, row.names=F)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks.corr <- c(seq(min(corr, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(corr, na.rm=TRUE)/paletteLength, max(corr, na.rm=TRUE), length.out=floor(paletteLength/2)))

corr <- data.frame(corr)
corrsub <- corr[!is.na(corr$Sex),]
corrnocol0 <- corrsub %>% keep(~!all(is.na(.x)))

png(width = 18, height = 15, file=paste0(figuredir,"corr_predicted_pheatmap_clust.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
pheatmap(corrnocol0, cluster_row = TRUE, cluster_col = TRUE, na_col = "grey90")
dev.off()
#annotation_col = annotation_col, annotation_colors = ann_colors

#now for just sig vars
allvarscorr <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))
names(allvarscorr)[5] <- "var_explained"
allvarscorr_1 <- subset(allvarscorr, var_explained>=0.01)

preds1 <- preds[,colnames(preds) %in% allvarscorr_1$variable]
cytokinevars <- c("IL5_co","IL13_co","IFNG_co","IL5_hc","IL13_hc","IFNG_hc")
preds1 <- preds1[,!colnames(preds1) %in% c("csex1","cage1","Sex","genPC1","genPC2","genPC3","cwght1","chght1","cgpd5","cgpd","cbpd",cytokinevars)]
corr <- cor(preds1, use="pairwise.complete.obs")
library("psych")
corr.psych <- corr.test(preds1, adjust="none")
# blank out non-significant correlations:
corr <- corr*(corr.psych$p<0.05)
pairs <- melt(corr)
pairs <- pairs[!pairs$value %in% c(0,1),]
corr <- data.frame(corr)
#corrsub <- corr[!is.na(corr$Sex),]
corrnocol0 <- corr %>% keep(~!all(is.na(.x)))

genedfcols <- merge(data.frame(variable=colnames(corrnocol0)),variables_df,by="variable")
genedfrows <- merge(data.frame(variable=rownames(corrnocol0)),variables_df,by="variable")
rownames(corrnocol0) <- genedfrows$description
colnames(corrnocol0) <- genedfcols$description


png(width = 18, height = 15, file=paste0(figuredir,"corr_predicted_pheatmap_clust_sig.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
pheatmap(corrnocol0, cluster_row = TRUE, cluster_col = TRUE, na_col = "grey90",breaks=myBreaks.corr, color=myColor)
dev.off()






################testing

/nfs/rprdata/ALOFT/AL/GLMnet/alpha0.1-LOO_119_rmXY


corrs <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))

