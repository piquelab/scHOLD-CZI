#################### INTERACTION -- THIS VERSION WITHOUT SUBSETTING TO EQTLs
#The interaction term is a tab-delimited text file or dataframe mapping sample ID to interaction value(s) 
#(if multiple interactions are used, the file must include a header with variable names).
#each row is sample id, each col is variable


#make each int file in R
R
library(data.table)
library(plyr);library(dplyr)

residual_path="residuals_ctrlonly"

data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/"
treat="CTRL"
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.1
dimset=50
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
variables_dftorun <- subset(variables_df, !variable %in% c("csex1","Sex","cage1","Wave","genPC1","genPC2","genPC3"))
fwrite(variables_dftorun, file=paste0(data_path,"variables_dftorun.txt"),sep="\t",col.names=F,row.names=F, quote=F)
#running only signature sig variables
cats <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/aloft_variables_categories.txt")
catsrm <- subset(cats, category %in% c("other","puberty"))
filter="CTRLonly" #ALOFT
glmnetfolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/glmnet/")
normmethod="voom" 
glmnetvoomfolder=paste0(glmnetfolder,normmethod,"/")
allvarscorr <- fread(file=paste0(glmnetvoomfolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))
allvarscorr_s <- subset(allvarscorr, !variable %in% c(catsrm$variable) & !variable=="")
names(allvarscorr_s)[5] <- "var_explained"
allvarscorr_1 <- subset(allvarscorr_s, var_explained>=0.05)
variables_dftorun <- subset(variables_df, variable %in% allvarscorr_1$variable)
fwrite(variables_dftorun, file=paste0(data_path,"variables_dftorun_glmnet.txt"),sep="\t",col.names=F,row.names=F, quote=F)

preds <- fread(paste0(glmnetvoomfolder,"all-predictions.txt"))

outFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
clusters <- fread(paste0(outFolder,"ALL.0.1.50.cluster_celltype.txt"))
opfn <- paste0(outFolder,"ALL.",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
load(opfn)

outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/"
filenames_res <- list.files(paste0(data_path,"residuals_ctrlonly/")) #file list from directory
filenames_cov <- list.files(paste0(data_path,"residuals_ctrlonly/covariates/")) #file list from directory
best_df <- fread(file=paste0(outFolder,"bestPCs_table.txt"))
#
for(c in clusters$cluster[!clusters$cluster %in% c("C9","C10")]){
  for(v in variables_dftorun$variable){
    best.PCs <- subset(best_df, cluster==c)$PCs
    pheno <- fread(paste0(data_path,"residuals_ctrlonly/",filenames_res[grep(paste0("phenotypes.",c,".",treat,".residuals_voom.sort.bed.gz$"),filenames_res)]))
    covar <- fread(paste0(data_path,"residuals_ctrlonly/covariates/",filenames_cov[grep(paste0(c,".",treat,".PC1-",best.PCs,".covariates_voom.txt"),filenames_cov)]))

    phenoind <- colnames(pheno[,-c(1:4)])
    covarind <- colnames(covar)[-1]
    cluster_metadata_sce <- metadata_ls[[c]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata_t <- subset(cluster_metadata, treats==treat )
    covvar <- as.data.frame(cluster_metadata_t[,v])
    notna <- complete.cases(covvar)
    cvinds <- cluster_metadata_t$Sample_ID[notna]
    common_samples <- Reduce(intersect, list(phenoind,covarind,cvinds))

    colkeep <- colnames(preds) %in% c("Sample_ID",v)
    predsv <- na.omit(preds[,..colkeep])
    predsind <- predsv$Sample_ID
    common_samplesp <- Reduce(intersect, list(phenoind,covarind,predsind))
    predsub <-na.omit(subset(predsv, Sample_ID %in% common_samplesp))

    cluster_metadata_var <- na.omit(subset(cluster_metadata_t[,c("Sample_ID",v)], Sample_ID %in% common_samples))
    covvarv <- cluster_metadata_var[,-1, drop=F]
    rownames(covvarv) <- cluster_metadata_var$Sample_ID #NEED ROWNAMES
    names(covvarv)[1] <- v
    fwrite(covvarv,file=paste0(data_path,"residuals_ctrlonly/covariates/",c,".",treat,".",v,".for_tensorqtl_int.allsnps.txt"),sep="\t",col.names=T,row.names=T, quote=F)
    phenosubcov <- subset(pheno, select = c(colnames(pheno)[1:4],cluster_metadata_var$Sample_ID))
    fwrite(phenosubcov,paste0(data_path,"residuals_ctrlonly/","phenotypes.",c,".",treat,".",v,".residuals_voom.allsnps.sort.bed.gz"),sep="\t",col.names=T,row.names=F, quote=F)
    covarsubcov <- subset(covar, select = c(colnames(covar)[1],cluster_metadata_var$Sample_ID))
    fwrite(covarsubcov,paste0(data_path,"residuals_ctrlonly/covariates/",c,".",treat,".",v,".PC1-",best.PCs,".covariates_voom.allsnps.txt"),sep="\t",col.names=T,row.names=F, quote=F)

    #using signatures
    covvarv <- predsub[,-1, drop=F]
    rownames(covvarv) <- predsub$Sample_ID #NEED ROWNAMES
    names(covvarv)[1] <- v
    fwrite(covvarv,file=paste0(data_path,"residuals_ctrlonly/covariates/",c,".",treat,".",v,".for_tensorqtl_int.signature.allsnps.txt"),sep="\t",col.names=T,row.names=T, quote=F)
    phenosubcov <- subset(pheno, select = c(colnames(pheno)[1:4],predsub$Sample_ID))
    fwrite(phenosubcov,paste0(data_path,"residuals_ctrlonly/","phenotypes.",c,".",treat,".",v,".residuals_voom.signature.allsnps.sort.bed.gz"),sep="\t",col.names=T,row.names=F, quote=F)
    covarsubcov <- subset(covar, select = c(colnames(covar)[1],predsub$Sample_ID))
    fwrite(covarsubcov,paste0(data_path,"residuals_ctrlonly/covariates/",c,".",treat,".",v,".PC1-",best.PCs,".covariates_voom.signature.allsnps.txt"),sep="\t",col.names=T,row.names=F, quote=F)

  } #var
} #cluster

#HPC test
data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly"

out_path="${outFolder}/interaction_allsnps"
mkdir -p ${out_path}

#cluster="C0"
treat="CTRL"

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/CTRLonly/ALL.0.1.50.cluster_celltype.txt | egrep -v "C9|C10"`;do
PC=`grep -w $cluster ${outFolder}/bestPCs_table.txt | cut -f2`
for var in `cut -f1 ${data_path}/variables_dftorun_glmnet.txt`; do
#for chr in {1..22};do
  chr=NA #can run all chr at once
  njobs=`squeue -u fh8591 -r| wc -l`
  maxjobs=500
  while [ "$njobs" -gt "$maxjobs" ];do 
  echo waiting for jobspace, sleeping ...
  sleep 300 
  njobs=`squeue -u fh8591 -r| wc -l`
  done #end while
      if [ -s "${out_path}/${cluster}_${treat}_${var}.tensorqtlint.cis_qtl_top_assoc.txt.gz" ] ; then
      echo " Output already exists. Skipping..."
    else
    echo submitting cluster $cluster variable $var 
    sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",var="$var",PC="${PC}",data_path="$data_path",out_path="$out_path" ${data_path}/tensorQTL/src/run_tensorqtl_int.sh 
    sleep 1
  fi
#done #chr end
done #var end
done #cluster end
mkdir ${out_path}/logs
mv ${out_path}/*.log ${out_path}/logs/

#for signatures
out_path="${outFolder}/interaction_signature_allsnps"
mkdir -p ${out_path}

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/CTRLonly/ALL.0.1.50.cluster_celltype.txt | egrep -v "C9|C10"`;do
PC=`grep -w $cluster ${outFolder}/bestPCs_table.txt | cut -f2`
for var in `cut -f1 ${data_path}/variables_dftorun_glmnet.txt`; do
#for chr in {1..22};do
  chr=NA #can run all chr at once
  njobs=`squeue -u fh8591 -r| wc -l`
  maxjobs=500
  while [ "$njobs" -gt "$maxjobs" ];do 
  echo waiting for jobspace, sleeping ...
  sleep 300 
  njobs=`squeue -u fh8591 -r| wc -l`
  done #end while
      if [ -s "${out_path}/${cluster}_${treat}_${var}.tensorqtlint.cis_qtl_top_assoc.txt.gz" ] ; then
      echo " Output already exists. Skipping..."
    else
    echo submitting cluster $cluster variable $var 
    sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",var="$var",PC="${PC}",data_path="$data_path",out_path="$out_path" ${data_path}/tensorQTL/src/run_tensorqtl_int_signature_allsnps.sh 
    sleep 1
  fi
#done #chr end
done #var end
done #cluster end


#combine
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/interaction_signature_allsnps/"
#outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/interaction_signature/"
#outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/interaction/"
if (!file.exists(paste0(outFolder,"figures"))) dir.create(paste0(outFolder,"figures"), showWarnings=F)

clusters <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/CTRLonly/ALL.0.1.50.cluster_celltype.txt")

myDir <- outFolder #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep("tensorqtlint.cis_qtl_top_assoc.txt.gz", filenames)] #pick specific files from list
int <- ldply(lapply(clusters$cluster,function(cluster){
  vl <- ldply(lapply(variables_dftorun$variable,function(v){
    cat("running",cluster,v,"\n")
    filenamesc <- filenames[grepl(paste0(cluster,".",treat,".",v), filenames)] #pick specific files from list
    if(length(filenamesc)>0){
    data_names <- gsub(".tensorqtlint.cis_qtl_top_assoc.txt.gz", "", filenamesc) #remove file ending
    #shortnames <- gsub("[.].*", "", data_names)
    for(i in 1:length(filenamesc)) assign(data_names[i], fread(file.path(myDir, filenamesc[i]),header = T)[,analysis:=data_names[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
    #combining data
    all_chrs <- lapply(data_names, function(x) get(x)) #grab data from list of df names
    names(all_chrs) <- data_names
    singelchrrun <- ldply(all_chrs, data.frame)[,-1]
    #singelchrrun <- all_chrsdf[grep(".NA",all_chrsdf$.id),] #after testing decided to run all chr together for interaction test
    #table(singelchrrun$analysis,singelchrrun$pval_adj_bh<0.1)
    #pval_emt = pval_gi * tests_emt and pval_adj_bh=p.adjust(pval_emt,method="BH")
    singelchrrun <- transform(singelchrrun, genotype_pval_emt=pval_g*tests_emt, variable_pval_emt=pval_i*tests_emt)
    singelchrrun <- transform(singelchrrun, genotype_padj=p.adjust(genotype_pval_emt,method="BH"), variable_padj=p.adjust(variable_pval_emt,method="BH"),
      cluster=cluster,treat=treat,variable=v)
    singelchrrun <- singelchrrun %>% relocate(c(cluster,treat,variable))  # have to use ensgene as there was multi gene symbols
    #table(singelchrrun$analysis,singelchrrun$genotype_padj<0.1)
    #table(singelchrrun$analysis,singelchrrun$variable_padj<0.1)
    return(singelchrrun)
    } 
  }),data.frame)
}),data.frame)

fwrite(int, paste0(outFolder,treat,".GxE_abundance_perclus.txt"), sep='\t', quote=F, row.names=F)

tested <- as.data.frame(table(int$cluster,int$variable))
intsig <- subset(int,pval_adj_bh<0.1 )
tablesig <- as.data.frame(table(intsig$cluster,intsig$variable))
outtable <- merge(tested,tablesig,by=c("Var1","Var2"))
colnames(outtable) <- c("cluster","variable","tested","int_FDR10")
outtabledf <- merge(outtable,variables_dftorun)
fwrite(outtabledf[,c("cluster","variable","description","tested","int_FDR10")], paste0(outFolder,treat,".GxE_summary_perclus.txt"), sep='\t', quote=F, row.names=F)

#table with categories (variable | colorbar | GxE eGenes)
intsig <- merge(intsig,cluster_celltype,by="cluster")
int_sigc <- merge(intsig,cats,by="variable")
intsigdf <- ddply(int_sigc, c("description","category"), plyr::summarize,
    'GxE eGenes'=length(unique(phenotype_id)))

intsigdf <- ddply(int_sigc, c("description","category","cell_type"), plyr::summarize,
    'GxE eGenes'=length(unique(phenotype_id)))
intsigdfc <- reshape2::dcast(intsigdf,description + category ~ cell_type,value.var="GxE eGenes")
fwrite(intsigdfc[order(intsigdfc$category),],sep='\t', quote=F, row.names=F, col.names=T,paste0(outFolder,"GxE_table.",resset,".",dimset,".txt"))


library(ggrastr)
i="CTRL"
lapply(split(int,int$variable),function(v){
    
    v <- transform(v, cluster=as.factor(cluster))
    v <- v %>%
   group_by(cluster)%>%
   arrange(pval_emt) %>%
   mutate(observed=-log10(pval_emt), expected=-log10(ppoints(length(pval_emt))))

p0 <- ggplot(v, aes(x=expected, y=observed, color=cluster))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    #facet_grid(variable~treat, scales="free_y")+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    ggtitle(paste0(unique(v$variable)))+
    theme_bw()
    figfn <- paste0(outFolder,"figures/GxE_abundance_pvalues_",unique(v$variable),".pcl_qqplot.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()
})