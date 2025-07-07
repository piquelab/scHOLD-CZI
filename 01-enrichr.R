library(plyr);library(dplyr)
library(msigdbr)
library(clusterProfiler)
  require(splitstackshape)
  library(gridExtra)
library(condformat)
library(data.table)
library(ggplot2)
library(cowplot)

m_df <- bind_rows(msigdbr(species = "Homo sapiens", category = "H"),
                  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS"),
                  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"))
msigdbr_t2g_all = m_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
msigdbr_t2g_BP <- m_df %>% filter(gs_cat=="C5") %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

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
opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
load(opfn)
combatrun="income_PCs_sex_age_and_treats_adjusted"

puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
allvars <- psychvarstorun[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]
firstrunvars <- allvars[c(1:10)]
secondrunvars <- c(allvars[c(1:20)],"cage1")
} else if(job=="CZI"){
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",11,0.2) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
contrastdf <- data.frame(control=c("RNA-CTRL","RNA-LPS"),treatment=c("RNA-LPS","RNA-LPS-DEX"))
outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
figuredir=paste0(outFolder,"figures/")
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

opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.RData")
load(opfn)
combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
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
allPFAS <- psychvarstorun[c(33:53)]
combPFAS <- allPFAS[seq(1,21,3)]
PFASvars <- c(allPFAS,"sumPFAS","log_sumPFAS")
allvars <- c(reordered_psychvarstorun[c(1:82)],"age","Lead","sumPFAS","log_sumPFAS")
clusters=c("C0", "C1", "C2", "C3", "C4", "C5", "C6")
}

baserun="limma_var"
log2cpm=FALSE
voom=TRUE
cpmqqnorm=FALSE
iselcov=FALSE
WHRcov=FALSE
runPFAS=TRUE
if(voom ){
    run=paste0(baserun,"_voom")
} else if (cpmqqnorm){
    run=paste0(baserun,"_cpmqqnorm")
} else{
    run=baserun
}
if(iselcov){
    run=paste0(run,"_iselcov")
}
if(WHRcov){
    run=paste0(run,"_WHRcov")
}
outFolder=paste0(baseoutFolder,run,"/")

enricheroutFolder=paste0(outFolder,"enrichr/")
if (!file.exists(enricheroutFolder)) dir.create(enricheroutFolder, showWarnings=F)

figuredir=paste0(enricheroutFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

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
variables_df <- merge(variables_df,allvarsdf,by="variable",all=T)

} else if(job=="CZI"){
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
variables_df <- unique(merge(variables_df,allvarsdf,by="variable",all=T))
}

varsrun = subset(variables_df, variable %in% c(allPFAS))
statsfull <- unique(fread(paste0(outFolder,"summary.",run,".txt")))

zoom="allPFAS"
lapply(treatmentsfirst, function(i) {
  lapply(varsrun$variable, function(var){
    cat("running",i,var,"\n")
    stats <- subset(statsfull, treat==i & variable==var)
    stats_sig <- subset(stats, varDEGs>100)
    #ressub <- subset(res, treat==i & variable==var)
    dge <- fread(paste0(outFolder,"limmares.",i,".",var,".",run,".txt"))
    names(dge)[4] <- "gene_symbol"
    #first subset variables that have at least one cluster with 50DEGs
    if(dim(subset(stats, varDEGs>100))[1]<1){ 
        cat("not enough DEGs","\n")
        } else{
dge <- subset(dge, cluster %in% stats_sig$cluster)
ddf <- split(dge,dge$cluster)

list_em <- lapply(ddf,function(f){
#f <- ddf[[1]]
  n <- unique(as.character(f$cluster))
message(paste0("running ",n))
f_df <- subset(f,adj.P.Val<0.1) #& abs(logFC)>0.5)
f_df_up <- subset(f_df, logFC>=0)
f_df_up <- transform(f_df_up, direction="up")
f_df_down <- subset(f_df, logFC<0)
f_df_down <- transform(f_df_down, direction="down")
me <- lapply(list(down=f_df_down,up=f_df_up),function(u){
  #u <-f_df_down
  nn <- ifelse(mean(u$logFC,na.rm=T)<0,"downregulated","upregulated")
  if(!dim(u)[1]==0){
em <- enricher(u$gene_symbol, TERM2GENE=msigdbr_t2g_BP, universe=f$gene_symbol, pvalueCutoff =0.1)
  if(!is.null(em) & dim(em)[1]>0){
res <- dplyr::filter(em@result, p.adjust < 0.05)
res <- res[order(res$p.adjust),]
   fwrite(res, sep='\t', quote=F, row.names=F, file=paste0(enricheroutFolder,n,"_",nn,"_gse_enrichr_fullres.",i,".",var,".",run,".txt"))
  if(!dim(res)[1]==0 ){
   s <- cSplit(as.data.table(res)[c(1:10),], "geneID", "/",direction = "long")
   s_dge <- merge(s, dge[dge$cluster==n,c("gene_symbol","logFC","adj.P.Val")], by.x="geneID",by.y="gene_symbol")
   s_dge <- s_dge[order(s_dge$adj.P.Val,-abs(s_dge$logFC)),]
   names(s_dge)[c(10:11)] <- c("dge_log2FoldChange","dge_adj.P.Val")
   #s <- transform(s, Z=qnorm(pvalue)) #remember direction is nonsensical
   fwrite(s_dge, sep='\t', quote=F, row.names=F, file=paste0(outFolder,n,"_",nn,"_gse_enrichr.txt"))
ss <- ldply(lapply(split(s_dge,s_dge$Description),function(i){
      ds <- paste(i$geneID[c(1:5)],collapse=';')
      df <- cbind(i[1,c("Description","qvalue","pvalue","p.adjust")],data.frame(geneID=ds))
      return(df)
      }),data.frame)[,-1]
   fwrite(ss, sep='\t', quote=F, row.names=F, file=paste0(outFolder,n,"_",nn,"_gse_enrichrtop.",i,".",var,".",run,".txt"))
  if(nn=="upregulated"){ 
      hc <- "indianred1"
      lc <- "darkred"
      title <- "Upregulated DEGs"
  } else {
      hc <- "lightblue"
      lc <- "darkblue"
    title <- "Downregulated DEGs"
  }
  em <- dplyr::mutate(em, direction=nn)
 #png(width = 8, height = 10, file=paste0("figure/",n,"_",nn,"_dotplot_enrichr.png"), pointsize=12, 
 #       bg = "transparent", units = "in", res = 1200)
  #p <-dotplot(em, showCategory=5) + scale_color_gradient(low = lc, high = hc, na.value = "darkgrey") + ggtitle(title)
 # print(p)
 # dev.off()
}}}
if(exists("em")){
  if(dim(em)[1]>0){
  return(em)
}
}
  })
#fig <-plot_grid(plotlist=me, nrow=1, ncol=2, align="h")
#figfn <- paste(figuredir, n,".dotplot.png", sep="")
#png(figfn, width=3000, height=2000, res=180)
#print(fig)
#dev.off()
}) #cluster loop
#subvars <- list_em[sapply(list_em, nrow)>0]
list_em <- list_em[sapply(list_em, is.null)] <- NULL
if(!is.null(list_em)){
combined_gsea<-merge_result(unlist(list_em,recursive=F)) #must be a single level named list
figfn <- paste(figuredir,"me.dotplot.",i,".",var,".",run,".png", sep="")
png(figfn, width=3000, height=2000, res=180)
p <- dotplot(combined_gsea, showCategory=5)
print(p)
dev.off()
} else { cat("em is null \n")}
} #if loop vardeg>100
})
})


#dont sep up and down
zoom="allPFAS"
lapply(treatmentsfirst, function(i) {
  lapply(varsrun$variable, function(var){
    cat("running",i,var,"\n")
    stats <- subset(statsfull, treat==i & variable==var)
    stats_sig <- subset(stats, varDEGs>100)
    #ressub <- subset(res, treat==i & variable==var)
    dge <- fread(paste0(outFolder,"limmares.",i,".",var,".",run,".txt"))
    names(dge)[4] <- "gene_symbol"
    #first subset variables that have at least one cluster with 50DEGs
    if(dim(subset(stats, varDEGs>100))[1]<1){ 
        cat("not enough DEGs","\n")
        } else{
dge <- subset(dge, cluster %in% stats_sig$cluster)
ddf <- split(dge,dge$cluster)

list_em <- lapply(ddf,function(f){
#f <- ddf[[1]]
  n <- unique(as.character(f$cluster))
message(paste0("running ",n))
f_df <- subset(f,adj.P.Val<0.1) #& abs(logFC)>0.5)
em <- enricher(f_df$gene_symbol, TERM2GENE=msigdbr_t2g_BP, universe=f$gene_symbol, pvalueCutoff =0.1)
  if(!is.null(em) & dim(em)[1]>0){
res <- dplyr::filter(em@result, p.adjust < 0.05)
res <- res[order(res$p.adjust),]
   fwrite(res, sep='\t', quote=F, row.names=F, file=paste0(enricheroutFolder,n,"_gse_enrichr_fullres.",i,".",var,".",run,".txt"))
  if(!dim(res)[1]==0 ){
   s <- cSplit(as.data.table(res)[c(1:10),], "geneID", "/",direction = "long")
   s_dge <- merge(s, dge[dge$cluster==n,c("gene_symbol","logFC","adj.P.Val")], by.x="geneID",by.y="gene_symbol")
   s_dge <- s_dge[order(s_dge$adj.P.Val,-abs(s_dge$logFC)),]
   names(s_dge)[c(10:11)] <- c("dge_log2FoldChange","dge_adj.P.Val")
   #s <- transform(s, Z=qnorm(pvalue)) #remember direction is nonsensical
   fwrite(s_dge, sep='\t', quote=F, row.names=F, file=paste0(outFolder,n,"_gse_enrichr.txt"))
ss <- ldply(lapply(split(s_dge,s_dge$Description),function(i){
      ds <- paste(i$geneID[c(1:5)],collapse=';')
      df <- cbind(i[1,c("Description","qvalue","pvalue","p.adjust")],data.frame(geneID=ds))
      return(df)
      }),data.frame)[,-1]
   fwrite(ss, sep='\t', quote=F, row.names=F, file=paste0(outFolder,n,"_gse_enrichrtop.",i,".",var,".",run,".txt"))
 #png(width = 8, height = 10, file=paste0("figure/",n,"_",nn,"_dotplot_enrichr.png"), pointsize=12, 
 #       bg = "transparent", units = "in", res = 1200)
  #p <-dotplot(em, showCategory=5) + scale_color_gradient(low = lc, high = hc, na.value = "darkgrey") + ggtitle(title)
 # print(p)
 # dev.off()
}}
if(exists("em")){
  if(dim(em)[1]>0){
  return(em)
}
}
#fig <-plot_grid(plotlist=me, nrow=1, ncol=2, align="h")
#figfn <- paste(figuredir, n,".dotplot.png", sep="")
#png(figfn, width=3000, height=2000, res=180)
#print(fig)
#dev.off()
}) #cluster loop
#subvars <- list_em[sapply(list_em, nrow)>0]
list_em <- list_em[sapply(list_em, is.null)] <- NULL
if(!is.null(list_em)){
combined_gsea<-merge_result(unlist(list_em,recursive=F)) #must be a single level named list
figfn <- paste(figuredir,"me.dotplot.",i,".",var,".",run,".png", sep="")
png(figfn, width=3000, height=2000, res=180)
p <- dotplot(combined_gsea, showCategory=5)
print(p)
dev.off()
} else { cat("em is null \n")}
} #if loop vardeg>100
})
})





myDir <- enricheroutFolder #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames1 <- filenames[grep("upregulated_gse_enrichrtop.txt", filenames)] #pick specific files from list
filenames2 <- filenames[grep("downregulated_gse_enrichrtop.txt", filenames)] #pick specific files from list
filenames <- c(filenames1,filenames2)
data_names <- gsub("_gse_enrichrtop.txt", "", filenames) #remove file ending
shortnames <- data_names
for(i in 1:length(filenames)) assign(shortnames[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')[,analysis:=shortnames[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
#combining data
ddf <- lapply(shortnames, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(shortnames)
rm(list=ls(pattern="C[0-9]"))
diff_analysis <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
diff_analysis <- transform(diff_analysis, cluster=gsub('_[^_]*$',"",analysis),direction=gsub(".*_","",analysis),p.adjust=round(p.adjust,3))
diff_analysis <- subset(diff_analysis, !cluster=="")
diff_analysis <- diff_analysis[,c("cluster","Description","geneID","p.adjust","direction")]

lapply(split(diff_analysis,diff_analysis$cluster), function(i){
  #i <- split(diff_analysis,diff_analysis$cluster)[[1]]
  n <- unique(i$cluster)
ss <- condformat(i) %>%
  rule_text_color(geneID, ifelse(direction=="downregulated","darkblue", ifelse(direction=="upregulated", "red4", "black"))) %>%
  condformat2grob()
png(width = 16, height = 8, file=paste0(figuredir,"enrichrtop_table.",n,".png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
grid.arrange(ss)
dev.off()
})



