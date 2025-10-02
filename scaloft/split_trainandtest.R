#need to split gene expression data into training and testing data sets
#1. Import your gene expression matrix (rows as genes_cluster, columns as samples) and corresponding clinical data
#2. split data

#require(caTools)
#set.seed(101) 
#sample = sample.split(data[,1], SplitRatio = .7)
#train = subset(data, sample == TRUE)
#test  = subset(data, sample == FALSE)


library(tidyr)
library(tidyverse) #gave me issues 031725 -- hax node
library(glmnet)
library(data.table)
library(methods)
library(plyr)
library(ggplot2)
library(ggpubr)

future::plan(strategy = 'multicore', workers = 10) #had an issue: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers
options(future.globals.maxSize = 15 * 1024 ^ 3)

myvar="pedu"
treat="CTRL"

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
variables_dftorun <- subset(variables_df, !variable %in% c("csex1"))
#44 vars
firstrunvars <- variables_dftorun$variable[c(1:15)]
secondrunvars <- variables_dftorun$variable[c(16:30)]
thirdrunvars <- variables_dftorun$variable[c(31:length(variables_dftorun$variable))]

#aloft categories
cats <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/aloft_variables_categories.txt")


base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.1
dimset=50
filter <- "CTRLonly" #ALOFT used
resmethod="voom"
varexp_thres <- 0.05

#combatrun="income_PCs_sex_age_and_treats_adjusted"
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
combatrun="income_PCs_sex_age_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
load(opfn)
glmnetfolder=paste0(baseoutFolder,"glmnet/")
if (!file.exists(glmnetfolder)) dir.create(glmnetfolder, showWarnings=F)
outFolder=paste0(glmnetfolder,resmethod,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#u_eigenvec2 <- unique(fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt"))
u_eigenvec2 <- unique(fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/covariates/",filter,".eigenvec2_u.txt")))
#to match names need to change dash to .
u_eigenvec2 <- transform(u_eigenvec2, Sample_ID=gsub("-",".",Sample_ID))

#was cinfirming NA was due to gene-cluster combo
#new_DF <- data[rowSums(is.na(data)) > 0,]
#new_DF[new_DF$V1=="ENSG00000115977",] %>% summarise_all(~all(is.na(.)))
#data_na <- ldply(lapply(names(counts_ls), function(cluster) {
#  Res1 <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/",cluster,".",treat,".residuals_",resmethod,".txt"))
#  Res1$ensg_cluster <- paste0(Res1$V1,"_",cluster)
#  new_DF <- Res1[rowSums(is.na(Res1)) > 0,]  
#  return(new_DF)
#  }),data.frame)

#for (treat in c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")){ #only have ctrl so far
  treat="CTRL"
	cat("running",treat,"\n")
data <- lapply(names(counts_ls), function(cluster) {
	Res1 <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals_ctrlonly/",cluster,".",treat,".residuals_",resmethod,".txt"))
  Res1$mean_gene=rowMeans(Res1[,-1], na.rm = TRUE)
  Res1$ensg_cluster <- paste0(Res1$V1,"_",cluster)
	Res1 <- Res1 %>% dplyr::select(V1, ensg_cluster, mean_gene, everything())
	return(Res1)
	})
data_rlist <- rbind.fill(data)

# Function to replace NAs in a row with the mean of non-NA values in that row

# Apply the function to each row of the data frame
df_imputed <- ldply(lapply(1:nrow(data_rlist), function(i) {
  #cat("running",i,"\t")
  row=data_rlist[i,]
    na_indices <- is.na(row)
  if (any(na_indices)) {
    row[na_indices] <- row$mean_gene
  }
  return(row)
  }),data.frame)

> dim(df_imputed)
[1] 110592    211

df_imputed_sub <- ldply(lapply(1:nrow(data_rlist), function(i) {
  #cat("running",i,"\t")
  row=data_rlist[i,]
    na_indices <- is.na(row)
if(sum(na_indices)/length(na_indices)>0.5){ #tested multiple percents, all the same amount removed so discarding this
      row[na_indices] <- row$mean_gene
} else{
  row <- NULL
}
  return(row)
  }),data.frame)
dim(df_imputed_sub)
[1] 36685   211


#df_imputed <- apply(data_rlist, 1, replace_na_with_row_mean)

u_eigenvec2_resind <- subset(u_eigenvec2, Sample_ID %in% colnames(df_imputed))
cv_original <- u_eigenvec2_resind[ order(match(u_eigenvec2_resind$Sample_ID, colnames(df_imputed[,-c(1:3)]))), ]

identical(colnames(df_imputed[,-c(1:3)]), cv_original$Sample_ID)

# keep the full covariate and GE data before subsetting:
cvfull <- cv_original
datafull <- df_imputed[,-c(1:3)]
if (!file.exists(paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))) {
    cat("making file \n")
table <- data.frame(variable=NA, description=NA, cor=NA,pvalue=NA,improve=NA,N=NA)
fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))
}

#allvars <- ldply(lapply(firstrunvars, function(myvar){
lapply(c(variables_dftorun$variable), function(myvar){
  # subset data to only samples that have y:
  cv <- na.omit(cv_original,cols=myvar)
  # subset data make sure they're in the same order:
  data_cv <- df_imputed[,colnames(df_imputed) %in% c("V1","ensg_cluster",cv$Sample_ID)]
  identical(colnames(data_cv[,-c(1:2)]),cv$Sample_ID)
  #data_sub <- data_cv[complete.cases(data_cv),] #removes na rows -> left me with only C0, so trying back to setting to 0

  df <- data.frame(starting_genecluster=length(df_imputed$ensg_cluster),complete_genecluster=length(data_cv$ensg_cluster),starting_uniquegene=length(unique(data_cv$V1)),complete_uniquegene=length(unique(data_sub$V1)))
  fwrite(df, file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".", myvar, "testedgenenum.txt"), col.names=TRUE, quote=FALSE, sep="\t")

  # train glmnet model:
  cat(myvar,"training \n")
  y <- as.numeric(unlist(cv[,..myvar]))
  N <- length(y)
  x <- t(data_cv[,-c(1:2)])
  colnames(x) <- data_cv$ensg_cluster
  set.seed(12)
  mymodel <- cv.glmnet(x,y,family="gaussian", alpha=0.1, type.measure="mse", nfold=length(y),trace.it = TRUE)

  #Get lambda from model
  ind = which(mymodel$lambda==mymodel$lambda.min)
  coefs <- mymodel$glmnet.fit$beta[,ind]
  nzero <- coefs[which(coefs!=0)]
  lambda <- mymodel$lambda[ind]

  ensts <- names(nzero)
  genes <- as.character(names(nzero))

  intercept <- coef(mymodel)[1]
  enstsweights <- coefs[ensts]

  cvm <- mymodel$cvm[ind]
  cvm_null <- mymodel$cvm[1]
  improve <- 1-cvm/cvm_null

  # save the glmnet model: (enstsweights etc was not originally saved but is used for plotting)
  cat("saving \n")
  save(mymodel,enstsweights,genes,file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".", myvar, ".RData"))

  # make predictions on all, including dropouts and check correlation:
  rownames(datafull) <- df_imputed$ensg_cluster
  vardata <- datafull[ensts,]
  varpredictions <- enstsweights %*% as.matrix(vardata) + intercept
  varpredictions <- t(data.frame(varpredictions))
  #rownames(varpredictions) <- gsub("[.]", "-", rownames(varpredictions))
  identical(rownames(varpredictions), cvfull$Sample_ID)
  corr <- cor.test(varpredictions[,1], as.numeric(unlist(cvfull[,..myvar])))

  # save the names of genes that went into the model:
  if (length(genes)>0){
  	gene_df=data.frame(ensg_cluster=genes)
  fwrite(gene_df, sep="\t",file=paste0(outFolder,project,".",resset,".",dimset,".",treat, ".", myvar, "-pred-genes.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE)
  }
  # and the predictions:
  rownames(varpredictions) <- gsub("[.]", "-", rownames(varpredictions))
  fwrite(varpredictions, file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".", myvar, "pred.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")

  # plot the correlation between actual and predicted:
  cat("plotting \n")

  dfplot <- cbind(x=as.numeric(unlist(cvfull[,..myvar])), y=varpredictions[,1])
  p <- ggplot(dfplot, aes(x=x, y=y)) +
    theme_bw()+
    xlab(paste0(variables_df[variables_df$variable==myvar,"description"]))+
    ylab(paste0("GLMnet-predicted ", paste0(variables_df[variables_df$variable==myvar,"description"])))+
    ggtitle(paste0("Correlation: ", myvar, " vs. GLMnet-predicted ", myvar, " on ", length(enstsweights), " gene_clusters"))+
    geom_point()+ #aes(color=sig)
    #geom_vline(xintercept = 0)+
    #geom_hline(yintercept = 0)+
    geom_smooth(method = "lm", se = FALSE)+
  #  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
    stat_cor(color="blue",method="pearson",cor.coef.name = "R", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
      theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
      axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
      axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
      legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
  png(width = 10, height = 10, file=paste0(figuredir,project,".",resset,".",dimset,".",treat,".", myvar,"_expanded_baseline_res_corr.png"), pointsize=12, 
        bg = "transparent", canvas = "white", units = "in", res = 1200)
  print(p)
  dev.off()

  # simple plot:
  png(width = 10, height = 10, file=paste0(figuredir,project,".",resset,".",dimset,".",treat, ".GLMnet-model-",myvar,".png"), pointsize=12, 
        bg = "transparent", canvas = "white", units = "in", res = 1200)
  plot(mymodel, main=paste0(length(genes), " genes best predict ", myvar))
  dev.off()

  # add official variable name and save correlation:
  tab <- data.frame(variables_df[variables_df$variable==myvar,],cor=corr$estimate,pvalue=corr$p.value,improve, N)
  fwrite(tab, file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"), sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE) #if rerunning, make sure to delete first or doubles variables results

  #return(tab)
})#,data.frame)

#fwrite(allvars, file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#} #treatment loop

#going back in and plotting splitting scaip1 and 2 as colors and equations
#run script through making cvfull then load:
lapply(variables_dftorun$variable, function(myvar){
varpredictions <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".", myvar, "pred.txt"))
load(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".", myvar, ".RData"))

dfplot <- cbind(as.numeric(unlist(cvfull[,..myvar])), varpredictions[,2],wave=as.factor(ifelse(cvfull$SCAIP7_18==0,"SCAIP1","SCAIP2")))
colnames(dfplot) <- c("x","y","wave")
p <- ggplot(dfplot, aes(x=x, y=y,color=wave, group=wave)) +
  theme_bw()+
  xlab(paste0(variables_df[variables_df$variable==myvar,"description"]))+
  ylab(paste0("GLMnet-predicted ", paste0(variables_df[variables_df$variable==myvar,"description"])))+
  ggtitle(paste0("Correlation: ", myvar, " vs. GLMnet-predicted ", myvar, " on ", length(enstsweights), " gene_clusters"))+
  geom_point()+ #aes(color=sig)
  #geom_vline(xintercept = 0)+
  #geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", se = FALSE)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(aes(color=wave),method="pearson",cor.coef.name = "R", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 10, height = 10, file=paste0(figuredir,project,".",resset,".",dimset,".",treat,".", myvar,"_expanded_baseline_res_corr_scaip1v2color.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
})

#larger table containing gene nums, plus corr from glmnet, data from Justyna? , #degs
filter <- "CTRLonly" #ALOFT used
deseqrun="income_PCs_sex_age_adjusted_withWave"
#deseqrun="income_PCs_sex_age_and_treats_adjusted_withWave" #old
#deseqoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/",deseqrun,"/") #old
deseqoutFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/",deseqrun,"/")

genenum <- ldply(lapply(c(variables_dftorun$variable,"age"), function(myvar){
  if(file.exists(paste0(outFolder,project,".",resset,".",dimset,".",treat, ".", myvar, "-pred-genes.txt"))){
  df <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".", myvar, "testedgenenum.txt"))
  #fixed script for future runs ignore next 2 lines
  #df <- df[,-1]
  colnames(df) <- c("starting_genecluster","complete_genecluster","starting_uniquegene","complete_uniquegene")
  modelgenes <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat, ".", myvar, "-pred-genes.txt"), header=F)
  df_n <- cbind(variable=myvar,df, modelgenenum=dim(modelgenes)[1])
  return(df_n)
  }
}),data.frame)
degnum <- ldply(lapply(c(variables_df$variable,"age"), function(var){
    cat("running",treat,var,"\n")
    if(file.exists(paste0(deseqoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",treat,".",deseqrun,".txt"))){
    stats <- fread(paste0(deseqoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",treat,".",deseqrun,".txt"))
    stats_s <- data.frame(variable=var,totalDEGs=sum(stats$DEGs_FDR_10,na.rm=T))
    return(stats_s)
    }
}),data.frame)
allvarscorr <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))
names(allvarscorr)[5] <- "var_explained"
justyna <- fread(paste0(base,"justyna_glmnet_1ctable.txt"))
colnames(justyna) <- gsub("full","justyna",colnames(justyna))
names(justyna)[c(5)] <- c("N.training.justyna")
justyna$varexp.training.justyna <- justyna$varexp.training/100
justyna_lab <- fread(paste0(base,"justyna_variable_symb.txt"))
justyna_symb <- merge(justyna[,-4],justyna_lab[,-c(1)],by.y="Short_form_label",by.x="variable")

df <- merge(genenum,allvarscorr,by="variable")
df2 <- merge(df,degnum,by="variable",all=T)
df3 <- merge(df2,justyna_symb[,c("symb","corr.justyna","pvalue.justyna","varexp.training.justyna","N.training.justyna")],by.x="variable",by.y="symb",all=T)

df4 <- df3 %>% dplyr::select(c(variable,description,complete_genecluster,N,modelgenenum,cor,pvalue,var_explained,totalDEGs,corr.justyna,pvalue.justyna,varexp.training.justyna,N.training.justyna))
fwrite(df4, file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".statstable_glmnet.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

df4 <- fread(paste0(outFolder,project,".",resset,".",dimset,".",treat,".statstable_glmnet.txt"))
sigus <- df4[df4$var_explained>=varexp_thres,]
sigjust <- df4[df4$varexp.training.justyna>=0.01,]
justna <- df4[is.na(df4$varexp.training.justyna),]
sigusonly <- sigus[!sigus$variable %in% sigjust$variable,]
sigusonlynotjustna <- sigusonly[!sigusonly$variable %in% justna$variable,]

library(ggrepel)
df4sig <- subset(df4,var_explained>=varexp_thres & !is.na(varexp.training.justyna))
  p <- ggplot(df4sig, aes(x=varexp.training.justyna, y=var_explained,color=variable)) +
    theme_bw()+
    geom_point(size=5)+ #aes(color=sig)
    geom_abline()+
    geom_hline(yintercept=varexp_thres,linetype="dashed", color = "red")+
    geom_vline(xintercept=0.01,linetype="dashed", color = "red")+
    geom_label_repel(aes(label = description),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
    #geom_smooth(method = "lm", se = FALSE)+
    stat_cor(color="blue",method="pearson",cor.coef.name = "R", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
      theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
      axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
      axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
      legend.position="none",strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
  png(width = 10, height = 10, file=paste0(figuredir,project,".",resset,".",dimset,".",treat,".","_var_explained_vsJustyna.png"), pointsize=12, 
        bg = "transparent", canvas = "white", units = "in", res = 1200)
  print(p)
  dev.off()

#calculating cell proportion for each model

#file=paste0(outFolder,project,".",resset,".",dimset,".",treat, ".", myvar, "-pred-genes.txt")

genedf <- ldply(lapply(c(variables_dftorun$variable,"age"), function(myvar){
  if(file.exists(paste0(outFolder,project,".",resset,".",dimset,".",treat, ".", myvar, "-pred-genes.txt"))){
  modelgenes <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat, ".", myvar, "-pred-genes.txt"), header=F)
  
  modelgenes <- transform(modelgenes, ensg=sapply(strsplit(modelgenes$V1, "_"),function(y) y[1]),cluster=sapply(strsplit(modelgenes$V1, "_"),function(y) y[2]))
  proplist <- ldply(lapply(names(counts_ls),function(c){
    modelc <- subset(modelgenes, cluster==c)
    modelnotc <- subset(modelgenes, !cluster==c)
    prop=length(modelc$V1)/length(modelgenes$V1)
    df <- data.frame(variable=myvar, cluster=c, geneincluster=length(modelc$V1),modelgenenum=length(modelgenes$V1),cellprop=prop)
    return(df)
    }),data.frame)
  #return(proplist)
  }
}),data.frame)

#subset for variables in 1%explained
allvarscorr <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))
allvarscorr <- merge(allvarscorr,cats,by=c("variable","description"),all=T)
allvarscorr_s <- subset(allvarscorr, !variable %in% c("csex1") & !variable=="")
names(allvarscorr_s)[5] <- "var_explained"
allvarscorr_1 <- subset(allvarscorr_s, var_explained>=varexp_thres)

#number of sig variables (not empty and not duping sex)
length(unique(allvarscorr_1$variable))
length(unique(allvarscorr_s$variable)) #total variables
length(unique(allvarscorr_s[!is.na(allvarscorr_s$cor),]$variable)) #total variables tested 
paste(unique(allvarscorr_1$description),collapse=", ")
paste(unique(allvarscorr_s[is.na(allvarscorr_s$cor),]$description),collapse=", ")
length(unique(allvarscorr_s[!is.na(allvarscorr_s$cor) & !allvarscorr_s$category=="other",]$variable))#total vars not other and not NA
length(unique(allvarscorr_1[!allvarscorr_1$category=="other",]$variable))#total sig vars not other 

summary(allvarscorr_s$cor)
summary(allvarscorr_1$cor)

allvarscorr_1_sub <- subset(allvarscorr_1, !category %in% c("other","glucocorticoid"))
genedf_1 <- subset(genedf, variable %in% allvarscorr_1_sub$variable)
#cytokinevars <- c("IL5_co","IL13_co","IFNG_co","IL5_hc","IL13_hc","IFNG_hc")
#genedf_2 <- subset(genedf_1, !variable %in% c("csex1","cage1","Sex","genPC1","genPC2","genPC3","cwght1","chght1","cgpd5","cgpd","cbpd",cytokinevars) & !variable=="")
genedfcols <- merge(genedf_1,cats,by="variable",all.x=T)

top10m <- reshape2::dcast(cluster ~ description, data=genedfcols, value.var="cellprop",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$cluster
top10ms <- top10m[apply(top10m[,-1], 1, function(x) !all(x==0)),]
top10msa <- top10ms[, colSums(top10ms != 0, na.rm = TRUE) > 0]
myMat <- top10msa[,-1]

library(ggcorrplot)
png(width = 10, height = 10, file=paste0(figuredir,"cellprop_heatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2, digits = 1) +
scale_fill_gradient2(low = "white", high = "red", breaks=c(0, 0.75), limit=c(0, 0.75)) + labs(fill = "Cell prop")
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
dev.off()

#high cell prop
calculate_range <- function(x) {
  c(median(x,na.rm=T), max(x,na.rm=T), ifelse((max(x,na.rm=T)-median(x,na.rm=T))>0.2, 1, 0))
}
colrange <- as.data.frame(apply(myMat, 2, calculate_range))
rownames(colrange) <- c("median","max","specific")
keep <- colrange[3,] ==1
specific <- myMat[,keep]
calculate_importance <- function(x) {
  c(sum(x >0.3,na.rm=T),length(x),sum(x >0.3,na.rm=T)/length(x))
}
rowimportance <- as.data.frame(apply(myMat, 1, calculate_importance))
rownames(rowimportance) <- c("num_important_vars","num_vars","prop_important")


#variables not sig
allvarscorr_L1 <- subset(allvarscorr_s, var_explained<varexp_thres)
genedf_1 <- subset(genedf, variable %in% allvarscorr_L1$variable)
cytokinevars <- c("IL5_co","IL13_co","IFNG_co","IL5_hc","IL13_hc","IFNG_hc")
genedf_2 <- subset(genedf_1, !variable %in% c("csex1","cage1","Sex","genPC1","genPC2","genPC3","cwght1","chght1","cgpd5","cgpd","cbpd",cytokinevars))
genedfcols <- merge(genedf_2,variables_dftorun,by="variable")

top10m <- reshape2::dcast(cluster ~ description, data=genedfcols, value.var="cellprop",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$cluster
top10ms <- top10m[apply(top10m[,-1], 1, function(x) !all(x==0)),]
top10msa <- top10ms[, colSums(top10ms != 0, na.rm = TRUE) > 0]
myMat <- top10msa[,-1]

png(width = 10, height = 10, file=paste0(figuredir,"cellprop_heatmap_nonsig.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2, digits = 1) +
scale_fill_gradient2(low = "white", high = "red", breaks=c(0, max(genedfcols$cellprop)), limit=c(0, max(genedfcols$cellprop))) + labs(fill = "Cell prop")
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
dev.off()

#compare cellcount and voom vs last run
varexp_thres <- 0.05
previousoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/glmnet/voom/")
oldresset=0.2
previousallvarscorr <- fread(file=paste0(previousoutFolder,project,".",oldresset,".",dimset,".",treat,".GLMnet-correlations.txt"))
names(previousallvarscorr)[5] <- "var_explained"

merged_treat <- merge(allvarscorr_s,previousallvarscorr,by=c("variable","description"))
merged_treat_all <- merged_treat %>% mutate(sig = case_when(
var_explained.x>=varexp_thres & var_explained.y>=varexp_thres ~ "4BOTH_sig",
var_explained.x >= varexp_thres ~ "3current_sig",
var_explained.y >= varexp_thres ~ "2previous_sig",    
))
merged_treat_all$sig[is.na(merged_treat_all$sig)] <- "1Not_Sig"

#forplotting <- subset(merged_treat_all, variable %in% variables_dftorun[c(1:10),"variable"])
forplotting <- subset(merged_treat_all, !sig=="1Not_Sig" & variable %in% variables_dftorun[,"variable"])

cols <- c("4BOTH_sig" = "green","1Not_Sig" = "grey","3current_sig" = "red", "2previous_sig" = "blue")
p <- ggplot(forplotting, aes(x=var_explained.x, y=var_explained.y)) +
  #facet_wrap(.~variable)+
  theme_bw()+
  xlab("current sig")+
  ylab("previous sig")+
  geom_point(size=5,aes(color=sig))+ #aes(color=sig)     
  geom_vline(xintercept = 0.1, linetype="dotted",color="red")+
  geom_hline(yintercept = 0.1, linetype="dotted",color="red")+
  geom_abline()+
  geom_label_repel(aes(label = description),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50')+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/current_vs_previous.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

# corr matrix signatures
library(psych)
library(rlist)
library(pheatmap)
#1%var explained varaivles
sigvars <- subset(variables_dftorun, variable %in% allvarscorr_1$variable)
predlist <- lapply(variables_dftorun$variable, function(v) {
  varpredictions <- fread(paste0(outFolder,project,".",resset,".",dimset,".",treat,".", v, "pred.txt"))
  pred=varpredictions[,2]
  names(pred)[1] <- v
  rownames(pred) <- varpredictions$V1
  return(pred)
  })
preddf <- list.cbind(predlist)
corrmat <- corr.test(preddf,adjust="none",ci=F)
myMat <- corrmat$r
pMat <- corrmat$p
myMat[is.na(pMat)] <- NA
pMat[is.na(myMat)] <- NA

png(width = 17, height = 8, file=paste0(figuredir,treat,".signaturecor_heatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2, insig = "blank",p.mat = pMat, sig.level=0.1,digits = 1) +
scale_fill_gradient2(low = "blue", high =  "red", mid = "white", midpoint = 0)
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
dev.off()

png(width = 17, height = 8, file=paste0(figuredir,treat,".signaturecor_pheatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
pheatmap(myMat, cluster_row = TRUE, cluster_col = TRUE, na_col = "grey90")
dev.off()



