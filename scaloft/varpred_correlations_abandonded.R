library(data.table)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library("psych")

library(ComplexHeatmap)
cell_fun = function(j, i, x, y, w, h, fill){
    if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
    }
    if( abs(mat[i, j])>0.05){
    if (cor_p[i, j]  < 0.05 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
      grid.text(paste0(sprintf("%.2f", mat[i, j]),"**"), x, y, gp = gpar(fontsize = 12))
    } else if (cor_p[i, j]  <= 0.1 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
      grid.text(paste0(sprintf("%.2f", mat[i, j]),"*"), x, y, gp = gpar(fontsize = 12))
    }
      }
}
lgd_list = list(
    Legend( labels = c("<0.05", "<0.1"), title = "padj",
            graphics = list(
              function(x, y, w, h) grid.text("**", x = x, y = y,
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.text("*", x = x, y = y,
                                               gp = gpar(fill = "black")))
            ))

future::plan(strategy = 'multicore', workers = 10) #had an issue: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers
options(future.globals.maxSize = 30 * 1024 ^ 3)

treat="CTRL"
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
#resset=0.2
resset=0.1
dimset=50
filter <- "CTRLonly" #ALOFT used
resmethod="voom"
#combatrun="income_PCs_sex_age_and_treats_adjusted"
combatrun="income_PCs_sex_age_adjusted"
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
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
#aloft categories
cats <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/aloft_variables_categories.txt")
catsrm <- subset(cats, category %in% c("other","puberty"))
catsrmother <- subset(cats, category %in% c("other"))

varexp_thres <- 0.05

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

allvarscorr <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))
#preds <- fread(paste0(outFolder,"all-predictions.txt"))

colkeep <- !colnames(preds) %in% c("csex1")
preds1 <- preds[,..colkeep]

corr <- cor(preds1, use="pairwise.complete.obs")
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

#complexheat
colkeep <- !colnames(preds) %in% c("Sample_ID",catsrm$variable)
preds1 <- preds[,..colkeep]

corr.psych <- corr.test(preds1, adjust="none")
corr.psych.corr <- corr.psych$r

genedfcols <- merge(data.frame(variable=colnames(corr.psych.corr)),variables_df,by="variable")
genedfrows <- merge(data.frame(variable=rownames(corr.psych.corr)),variables_df,by="variable")
rownames(corr.psych.corr) <- genedfrows$description
colnames(corr.psych.corr) <- genedfcols$description

giveNAs = which(is.na(as.matrix(dist(corr.psych.corr))),arr.ind=TRUE)
head(giveNAs)
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
sum(is.na(as.matrix(dist(corr.psych.corr[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]
if(!is.null(rmv)){
    colkeep <- colnames(corr.psych.corr) %in% rownames(corr.psych.corr[-as.numeric(rmv),])
    mat = corr.psych.corr[-as.numeric(rmv),colkeep]
    cor_p = corr.psych$p.adj[-as.numeric(rmv),colkeep]
} else{
    mat = corr.psych.corr
    colkeep <- colnames(corr.psych$p.adj) %in% genedfcols$variable
    cor_p = corr.psych$p.adj[,colkeep]    
}

hp<- ComplexHeatmap::Heatmap(mat,
                        rect_gp = gpar(type = "none"),
                        column_dend_side = "bottom",
                        #column_title = "NK cells",
                        name = "correlation", col = myColor,
                        cell_fun = cell_fun,
                        cluster_rows = T, cluster_columns = T,
                        row_names_side = "left")

png(width = 15, height = 12, file=paste0(figuredir,"corr_predicted_complexheat_clust.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
draw(hp, annotation_legend_list = lgd_list, ht_gap = unit(1, "cm") )
dev.off()

#now for just sig vars
allvarscorr_s <- subset(allvarscorr, !variable %in% c(catsrm$variable,"pedu") & !variable=="") #currently removing education
names(allvarscorr_s)[5] <- "var_explained"
allvarscorr_1 <- subset(allvarscorr_s, var_explained>=varexp_thres)

colkeep <- colnames(preds1) %in% allvarscorr_1$variable
preds1 <- preds1[,..colkeep]
corr <- cor(preds1, use="pairwise.complete.obs")
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

#got this from https://stackoverflow.com/questions/61469201/pheatmap-won-t-cluster-rows-na-nan-inf-in-foreign-function-call-arg-10
#was getting an error (between some rows, it's not possible to calculate euclidean distances. You need to the euclidean distance matrix to have no NAs to do clustering) and this fixes it
giveNAs = which(is.na(as.matrix(dist(corrnocol0))),arr.ind=TRUE)
head(giveNAs)
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
sum(is.na(as.matrix(dist(corrnocol0[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]
mat = corrnocol0[-as.numeric(rmv),]
cor_p = corr.psych$p.adj[-as.numeric(rmv),colnames(corr.psych$p.adj) %in% genedfcols$variable]
cor_mat = corr.psych$r[,colnames(corr.psych$r) %in% genedfcols$variable]
rownames(cor_mat) <- genedfrows$description
colnames(cor_mat) <- genedfcols$description
cor_mat = cor_mat[-as.numeric(rmv),]
myBreaks.corr <- c(seq(min(corr, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(corr, na.rm=TRUE)/paletteLength, max(corr, na.rm=TRUE), length.out=floor(paletteLength/2)))

png(width = 18, height = 15, file=paste0(figuredir,"corr_predicted_pheatmap_clust_sig.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
pheatmap(mat, cluster_row = TRUE, cluster_col = TRUE, na_col = "grey90",breaks=myBreaks.corr, color=myColor, fontsize = 16)
dev.off()

corr.psych.corr <- corr.psych$r

genedfcols <- merge(data.frame(variable=colnames(corr.psych.corr)),variables_df,by="variable")
genedfrows <- merge(data.frame(variable=rownames(corr.psych.corr)),variables_df,by="variable")
rownames(corr.psych.corr) <- genedfrows$description
colnames(corr.psych.corr) <- genedfcols$description

giveNAs = which(is.na(as.matrix(dist(corr.psych.corr))),arr.ind=TRUE)
head(giveNAs)
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
sum(is.na(as.matrix(dist(corr.psych.corr[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]
if(!is.null(rmv)){
    colkeep <- colnames(corr.psych.corr) %in% rownames(corr.psych.corr[-as.numeric(rmv),])
    mat = corr.psych.corr[-as.numeric(rmv),colkeep]
    cor_p = corr.psych$p.adj[-as.numeric(rmv),colkeep]
} else{
    mat = corr.psych.corr
    colkeep <- colnames(corr.psych$p.adj) %in% genedfcols$variable
    cor_p = corr.psych$p.adj[,colkeep]    
}

hp<- ComplexHeatmap::Heatmap(mat,
                        rect_gp = gpar(type = "none"),
                        row_names_gp = gpar(fontsize = 12), 
                        column_names_gp = gpar(fontsize = 12),
                        column_dend_side = "bottom",
                        #column_title = "NK cells",
                        name = "correlation", col = myColor,
                        cell_fun = cell_fun,
                        cluster_rows = T, cluster_columns = T,
                        row_names_side = "left")

png(width = 12, height = 10, file=paste0(figuredir,"corr_predicted_complexheat_clust_sig.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 640)
draw(hp, annotation_legend_list = lgd_list, ht_gap = unit(1, "cm") )
dev.off()

############################
############################

#now to do correlation of the variables
metaoutFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
scmetadata <- fread(paste0(metaoutFolder,"scmetadata_allind.txt"))

#now for just sig vars
glmnetfolder=paste0(baseoutFolder,"glmnet/")
outFolder=paste0(glmnetfolder,resmethod,"/")
allvarscorr <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))
allvarscorr_s <- subset(allvarscorr, !variable %in% c(catsrm$variable) & !variable=="")
names(allvarscorr_s)[5] <- "var_explained"
allvarscorr_1 <- subset(allvarscorr_s, var_explained>=varexp_thres)

colkeep <- colnames(scmetadata) %in% allvarscorr_1$variable
scmetadata1 <- scmetadata[,..colkeep]
library("psych")
corr.psych <- corr.test(scmetadata1, adjust="none")
corr.psych.corr <- corr.psych$r

genedfcols <- merge(data.frame(variable=colnames(corr.psych.corr)),variables_df,by="variable")
genedfrows <- merge(data.frame(variable=rownames(corr.psych.corr)),variables_df,by="variable")
rownames(corr.psych.corr) <- genedfrows$description
colnames(corr.psych.corr) <- genedfcols$description

giveNAs = which(is.na(as.matrix(dist(corr.psych.corr))),arr.ind=TRUE)
head(giveNAs)
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
sum(is.na(as.matrix(dist(corr.psych.corr[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]
if(!is.null(rmv)){
mat = corr.psych.corr[-as.numeric(rmv),]
colkeep <- colnames(corr.psych$p.adj) %in% genedfcols$variable
cor_p = corr.psych$p.adj[-as.numeric(rmv),colkeep]
} else{
mat = corr.psych.corr
colkeep <- colnames(corr.psych$p.adj) %in% genedfcols$variable
cor_p = corr.psych$p.adj[,colkeep]    
}

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks.corr <- c(seq(min(corr, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(corr, na.rm=TRUE)/paletteLength, max(corr, na.rm=TRUE), length.out=floor(paletteLength/2)))

hp<- ComplexHeatmap::Heatmap(mat,
                        rect_gp = gpar(type = "none"),
                        column_dend_side = "bottom",
                        #column_title = "NK cells",
                        name = "correlation", col = myColor,
                        cell_fun = cell_fun,
                        cluster_rows = T, cluster_columns = T,
                        row_names_side = "left")

png(width = 15, height = 12, file=paste0(figuredir,"corr_basevariables_complexheat_clust_sig.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
draw(hp, annotation_legend_list = lgd_list, ht_gap = unit(1, "cm") )
dev.off()


#all not just sig
colkeep <- colnames(scmetadata) %in% allvarscorr_s$variable
scmetadata1 <- scmetadata[,..colkeep]
corr.psych <- corr.test(scmetadata1, adjust="none")
corr.psych.corr <- corr.psych$r

genedfcols <- merge(data.frame(variable=colnames(corr.psych.corr)),variables_df,by="variable")
genedfrows <- merge(data.frame(variable=rownames(corr.psych.corr)),variables_df,by="variable")
rownames(corr.psych.corr) <- genedfrows$description
colnames(corr.psych.corr) <- genedfcols$description

giveNAs = which(is.na(as.matrix(dist(corr.psych.corr))),arr.ind=TRUE)
head(giveNAs)
tab = sort(table(c(giveNAs)),decreasing=TRUE)
checkNA = sapply(1:length(tab),function(i){
sum(is.na(as.matrix(dist(corr.psych.corr[-as.numeric(names(tab[1:i])),]))))
})
rmv = names(tab)[1:min(which(checkNA==0))]
if(!is.null(rmv)){
mat = corr.psych.corr[-as.numeric(rmv),]
colkeep <- colnames(corr.psych$p.adj) %in% genedfcols$variable
cor_p = corr.psych$p.adj[-as.numeric(rmv),colkeep]
} else{
mat = corr.psych.corr
colkeep <- colnames(corr.psych$p.adj) %in% genedfcols$variable
cor_p = corr.psych$p.adj[,colkeep]    
}

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks.corr <- c(seq(min(mat, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1), seq(max(mat, na.rm=TRUE)/paletteLength, max(mat, na.rm=TRUE), length.out=floor(paletteLength/2)))

hp<- ComplexHeatmap::Heatmap(mat,
                        rect_gp = gpar(type = "none"),
                        column_dend_side = "bottom",
                        #column_title = "NK cells",
                        name = "correlation", col = myColor,
                        cell_fun = cell_fun,
                        cluster_rows = T, cluster_columns = T,
                        row_names_side = "left")

png(width = 15, height = 12, file=paste0(figuredir,"corr_basevariables_complexheat_clust.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
draw(hp, annotation_legend_list = lgd_list, ht_gap = unit(1, "cm") )
dev.off()
