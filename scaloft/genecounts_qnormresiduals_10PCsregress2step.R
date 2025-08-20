
#this script was adapted to regress out the top 10 PCs for a QC check
R
library(tidyverse)
library(edgeR)
library(limma)
library(annotables)
library(data.table)
library(plyr);library(dplyr)

job="ALOFT"
#load in original counts info for clusters and treats used
if(job=="ALOFT"){
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
load(opfn)
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

# load annotation
u_eigenvec2 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt")
treatments=c("CTRL","LPS","PHA","PHA-DEX")

} else if(job=="CZI"){
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.15) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
outFolder=paste0(base,"fastQTL/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
figuredir=paste0(outFolder,"figures/")
combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
load(opfn)

leadvar <- fread("/rs/rs_grp_schold/covariates/other_covariates/HOLD LEAD 5.27.25.csv")
cov_pluslead <- merge(cov_file,leadvar,by="pID",all=T)
cov_pluslead <- transform(cov_pluslead, Lead=ifelse(Lead==-99,NA,Lead))
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o[,-c("sex","sex_alph","age")],cov_pluslead,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip","SNI_NoP","age","sex","sex_alph","isel","pID") #SNI_NoP is the only variable that should be excluded based on the observed issues in score distributions that don’t make sense (negative values and extreme outliers)
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
old_cytokines <- psychvarstorun[c(8,14:24)] #old cytokines
keep=!colnames(eigenvec2) %in% c(notrun_var,old_cytokines)
u_eigenvec2 <- eigenvec2[,..keep]
treatments=c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")

}
if (!file.exists(paste0(outFolder,"vcf/"))) dir.create(paste0(outFolder,"vcf/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"counts/"))) dir.create(paste0(outFolder,"counts/"), showWarnings=F)

cv <- unique(u_eigenvec2)
##check for individuals with more than one sample
cv %>% count(Sample_ID) %>% dplyr::filter(n>1)

library("AnnotationHub")
ah <- AnnotationHub()
if(length(ah["AH98047"]) == 0) {
  edb <- ah[["AH75011"]]
} else {
  edb <- ah[["AH98047"]]
}
geneIDs <- genes(edb) %>%
  as.data.frame() %>% 
  setDT(keep.rownames = "ensembl_gene_id") %>%
  .[, c("ensembl_gene_id","entrezid","symbol","seqnames","start","end","strand","gene_biotype", "description")]
names(geneIDs)[c(1,2,4,8)] <- c("ensgene","entrez","chr","biotype")
geneIDs <- subset(geneIDs, biotype=="protein_coding")
geneIDs <- transform(geneIDs, chr=as.character(chr),strand=as.character(strand))
#this 1bp positioning for bed file is based on examples from qtltools and tensorqtl
geneIDs <- transform(geneIDs, strand_start=ifelse(strand=="+",start,end),strand_end=ifelse(strand=="+",start+1,end+1))
#geneIDs <- transform(geneIDs, strand_start=start,strand_end=start+1)
geneIDs <- subset(geneIDs, chr %in% c(1:22))

lapply(names(counts_ls),function(clus){
  for (i in treatments){
    cluster_metadata_sce <- metadata_ls[[clus]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata_t <- subset(cluster_metadata, treats==i)
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",clus,".",combatrun,".RData")
    load(opfn)
    if(job=="ALOFT"){
    adjusted_counts <- adjusted
    }
    cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_t))]
    cluster_metadata_t <- cluster_metadata_t[which(rownames(cluster_metadata_t) %in% colnames(cluster_counts_t)),]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_t))
    genes <- rownames(cluster_counts_t)
    colnames(cluster_counts_t) <- cluster_metadata_t$Sample_ID
    cluster_counts_tg <- cbind(genes, as.data.frame(cluster_counts_t))
    all_counts_bed <- merge(geneIDs[,c("chr","strand_start","strand_end","symbol","ensgene")],cluster_counts_tg,by.x="symbol",by.y="genes")
    all_counts_bed <- all_counts_bed %>% relocate(ensgene,.after =strand_end) %>% dplyr::select(-symbol) # have to use ensgene as there was multi gene symbols
    #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
    names(all_counts_bed)[c(2:4)] <- c("start","end","gene_id")
    fwrite(all_counts_bed, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"counts/phenotypes.",clus,".",i,".bed"))
    samples <- data.frame(samples=colnames(all_counts_bed[,-c(1:4)]))
    fwrite(samples, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(outFolder,"vcf/sample_list_fastqtl.",clus,".",i,".txt"))
    #  cat(clus,i, length(samples$sample),"\n")
  }
})
ctrl-z #background while running vcfs 

#subset vcf to samples in data
module load bcftools/1.19
#ALOFT
bcftools view --header-only /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.vcf.gz > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.txt
sed 's/RI\_//g' /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.txt > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.gsubRI.txt
bcftools reheader -h /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.gsubRI.txt /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.vcf.gz > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf
#takes awhile:
bcftools view -S /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/sample_list_fastqtl.txt /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf
bgzip /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf && tabix -p vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf.gz

#vcf per treatment and cluster
treat="CTRL"
for cluster in `awk 'NR>2 && NR<8{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
#for cluster in `awk 'NR>7 {print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
for cluster in `awk 'NR>11{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
echo running $cluster
bcftools view -S /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/sample_list_fastqtl.$cluster.$treat.txt /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf
echo zipping
bgzip --threads 10 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf && tabix -p vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz
echo making header
bcftools view --header-only /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz | sed '/^##/d' | cut -f10- > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.vcfheader.txt
done

#not functional yet
#data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
#treat="CTRL"
#for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt | grep -v "C4"`;do
#  echo submitting cluster $cluster
#  sbatch --export=cluster="${cluster}",treat="${treat}",data_path="$data_path" ${data_path}/src/vcf_subset.sh 
#  sleep 1
#done



#CZI
mkdir /rs/rs_grp_schold/CZI/RNA/analysis/vcf
bcftools view --header-only /rs/rs_grp_schold/CZI/RNA/analysis/ref.ac1.removekin.vcf.gz > /rs/rs_grp_schold/CZI/RNA/analysis/vcf/ref.ac1.headeronly.txt
treat="RNA-CTRL"
for cluster in C0 C1 C2 C3 C4 C5 C6 ;do
echo running $cluster
bcftools view -S /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/sample_list_fastqtl.$cluster.$treat.txt /rs/rs_grp_schold/CZI/RNA/analysis/ref.ac1.removekin.vcf.gz > /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/vcf/ref.ac1.removekin.$cluster.$treat.filtered.vcf
echo zipping
bgzip /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/vcf/ref.ac1.removekin.$cluster.$treat.filtered.vcf && tabix -p vcf /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/vcf/ref.ac1.removekin.$cluster.$treat.filtered.vcf.gz
echo making header
bcftools view --header-only /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/vcf/ref.ac1.removekin.$cluster.$treat.filtered.vcf.gz | sed '/^##/d' | cut -f10- > /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/vcf/ref.ac1.removekin.$cluster.$treat.filtered.vcfheader.txt
done

fg 

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

job="ALOFT"
voom=TRUE
cpmqqnorm=FALSE
cluster="C6"
if (!file.exists(paste0(outFolder,"covariates/"))) dir.create(paste0(outFolder,"covariates/"), showWarnings=F)

cv <- unique(u_eigenvec2)

clusters <- names(counts_ls)
treat=treatments[1]

lapply(clusters,function(cluster){
#  lapply(treatments,function(treat){
    cat("running ",cluster,treat,"\n")

countFile <- fread(paste0(outFolder,"counts/phenotypes.",cluster,".",treat,".bed"))
data1 <- countFile[,-c(1:3)]
data <- data1[data1$gene_id %in% geneIDs$ensgene,]
data[is.na(data)] <- 0

## Normalization of data
# make edgeR object
dge <- DGEList(counts=data)
#sum(colnames(dge$counts)==cv$Sample_ID) #does the order need to match or is this just a check for how many samples have all info
sum(colnames(dge$counts) %in% cv$Sample_ID)

#Transform counts to counts per million
dge <- calcNormFactors(dge)
cpm0 <- cpm(dge)
rownames(cpm0) <- data$gene_id
samples <- dim(data)[2]
#Remove genes that are lowly expressed
table(rowSums(dge$counts==0)==samples) #Shows how many transcripts have 0 count across all samples
# as did GTEx:
keep.exprs <- rowSums(cpm0>=0.1)>=(0.2*samples) 
#& rowSums(data>=6)>=(0.2*samples) #Only keep transcripts that have cpm>0.1 in at least 20% of the samples
#Julong does not use the 6 count min
# subset samples in cv to match samples in express data
cv_d <- subset(cv, Sample_ID %in% colnames(data))

## # Normalize data
if(voom){
method="voom"
dim(dge) #before filter
#dge <- dge[keep.exprs,]
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE] #this was the old script line, must have the wrong set of libraries loaded for it to run
dim(dge) #genes that pass filter
fwrite(data.frame(genes=dge$genes$gene_id), sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,cluster,".",treat,".residualsgenes_",method,".txt"))

if(job=="ALOFT"){
design <- model.matrix(~0+ as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
} else if (job=="CZI"){
design <- model.matrix(~0+ as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2) )  
}
v_e <- voom(dge, design, plot=FALSE, normalize.method="quantile")
genes_normed_baseline <- data.frame(v_e$E)
##
X <- design
H <- X %*% solve(t(X) %*% X) %*% t(X)
dim(H)
He <- (diag(rep(1,ncol(H)))-H)
Res <- as.matrix(genes_normed_baseline) %*% He
sum(abs(t(He)-He)) #Should be almost 0
#1.476147e-12 
#save the residuals:
Res <- data.frame(Res)
colnames(Res) <- colnames(dge)
rownames(Res) <- dge$genes$gene_id
write.table(Res, paste0(outFolder,cluster,".",treat,".residuals_",method,".txt"), sep="\t", row.names=TRUE, quote=FALSE)
}
if(cpmqqnorm){
dim(cpm0) #before filter
cpm <- cpm0[keep.exprs,]
#cpm <- cpm[keep.exprs,, keep.lib.sizes=FALSE] #this was the old script line, must have the wrong set of libraries loaded for it to run
dim(cpm) #genes that pass filter
method="qnorm"
#quantile normalize across rows - aka across samples for each gene
mat_qnorm <- apply(cpm,1,function(x){qqnorm(rank(x, ties.method = "random"), plot = F)$x})
#model 
design_expanded <- model.matrix(~0+ as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
# Fit the linear model
model <- lm(mat_qnorm ~ design_expanded, data = data.frame(mat_qnorm, cv_d))
# Extract residuals
residuals <- t(residuals(model))
Res1 <- data.frame(residuals)
colnames(Res1) <- colnames(cpm)
#rownames(Res1) <- cpm0$genes$gene_id
write.table(Res1, paste0(outFolder,cluster,".",treat,".residuals_",method,".txt"), sep="\t", row.names=TRUE, quote=FALSE)
}
})

## this script calculates PCs for the bed file used for FastQTL analysis
#IBD eQTL Project - Rectum samples - protein coding genes only, the residuals

#library(tidyverse)
library(irlba)
nPCs=20
if(voom){
    method="voom"
}else{
    method="qnorm"
}
lapply(clusters,function(cluster){
#  lapply(treatments,function(treat){
    cat("running ",cluster,treat,"\n")
## load normalized data / or residuals:
Res1 <- fread(paste0(outFolder,cluster,".",treat,".residuals_",method,".txt"))
if(job=="ALOFT"){
vcfind <- colnames(fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/","vcf/ref.ac1.",cluster,".",treat,".vcfheader.txt")))
} else if (job=="CZI"){
  vcfind <- colnames(fread(paste0(outFolder,"vcf/ref.ac1.removekin.",cluster,".",treat,".filtered.vcfheader.txt")))
}
Res <- Res1[,-1]

# run PCA on residuals
PCs <- try(prcomp_irlba(t(Res), n=nPCs))
if(inherits(PCs, "try-error")){
  print("An error occurred")
} else {

summary(PCs)
mypcs <- as.data.frame(PCs$x)
rownames(mypcs) <-colnames(Res)
samples <- colnames(Res)

covs <- as_tibble(t(mypcs))
#order to vcf
covsord <- covs[,match(vcfind, colnames(covs)) ]
cvs <- cbind(id=colnames(mypcs),covsord)

fwrite(cvs, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"covariates/PCcovariates_",method,".",cluster,".",treat,".txt"))
fwrite(cvs, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(outFolder,"covariates/PCcovariates_",method,"nohead.",cluster,".",treat,".txt"))

PCvar <- data.frame(summary(PCs)$importance)[2,]
PCvar <- t(PCvar)

# save the proportion of variance expplained by PCs:
fwrite(PCvar, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"covariates/var_explained_by_GE_PCs_",method,".",cluster,".",treat,".txt"))
}
})
ctrl-z #background while grabbing PCs

#had 1-30 PCs but unnecessary to do that many
treat="CTRL"
method="voom"
#ALOFT
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
  for i in $(seq 1 20); do 
  if [ -f /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/covariates/PCcovariates_${method}.$cluster.$treat.txt ]; then 
  head -n $(($i+1)) /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/covariates/PCcovariates_${method}.$cluster.$treat.txt > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/covariates/$cluster.$treat.PC1-$i.covariates_${method}.txt
  fi
  done
done

#CZI
treat="RNA-CTRL"
for cluster in C0 C1 C2 C3 C4 C5 ;do
  for i in $(seq 1 20); do 
  if [ -f /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/covariates/PCcovariates_${method}-FastQTL.$cluster.$treat.txt ]; then 
  head -n $(($i+1)) /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/covariates/PCcovariates_${method}-FastQTL.$cluster.$treat.txt > /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/covariates/$cluster.$treat.PC1-$i.covariates_${method}-FastQTL.txt
  fi
  done
done


#step 3 regress out PCs
fg 

lapply(names(counts_ls),function(cluster){
#  lapply(treatments,function(treat){
for (pcnum in seq(1:20)){
    if(isTRUE(file.size(paste0(outFolder,"covariates/",cluster,".",treat,".PC1-",pcnum,".covariates_",method,".txt")) > 0)){
      cat("running ",cluster,treat,pcnum,"\n")
    PCs <- fread(paste0(outFolder,"covariates/",cluster,".",treat,".PC1-",pcnum,".covariates_",method,".txt"))
    design_expanded <- model.matrix(~0+ t(PCs[,-1]))
    Res1 <- fread(paste0(outFolder,cluster,".",treat,".residuals_",method,".txt"))
    Resnoid <- Res1[,-1]
    cv_d <- subset(cv, Sample_ID %in% colnames(Resnoid))
    resgenes <- fread(paste0(outFolder,cluster,".",treat,".residualsgenes_",method,".txt"))$genes
    model <- lm(t(Resnoid) ~ design_expanded, data = data.frame(t(Resnoid), cv_d))
    # Extract residuals
    residuals <- t(residuals(model))
    Res <- data.frame(residuals)
    colnames(Res) <- colnames(Resnoid)
    rownames(Res) <- resgenes
    # save the residuals
    write.table(Res, paste0(outFolder,cluster,".",treat,".PC1-",pcnum,".residuals_",method,"_PCregress2step.txt"), sep="\t", row.names=TRUE, quote=FALSE)

    #add bed file info
    Res$ensgene <- resgenes
    all_counts_bed <- merge(geneIDs[,c("chr","strand_start","strand_end","ensgene")],Res,by="ensgene")
    all_counts_bed <- all_counts_bed %>% relocate(ensgene,.after =strand_end)  # have to use ensgene as there was multi gene symbols
    #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
    names(all_counts_bed)[c(2:4)] <- c("start","end","gene_id")
    #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
    fwrite(all_counts_bed, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"phenotypes.",cluster,".",treat,".PC1-",pcnum,".residuals_",method,"_PCregress2step.bed"))
  }
}
})
q
()

#testing the results
#Sample IDs are specified in the header line. This line needs to start with a hash key (i.e. #).
#rm /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.C*.CTRL.PC1-*.residuals_qnorm_10PCregress2step.sort.bed.gz
#rm /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.C*.CTRL.PC1-*.residuals_qnorm_10PCregress2step.sort.bed.gz*
module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0
method="voom"
#ALOFT
treat="CTRL"
for cluster in C0 C1 C10 C11 C2 C3 C4 C5 C6 C7 C8 C9;do
for pcnum in $(seq 1 20); do
i=`ls -1 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/phenotypes.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.bed | grep -v 'sort'`
  echo "running " $i
  #less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sortBed -i"}' > ${i%.*}.sort.bed
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sort -k1,1 -k2,2n "}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
done
done

#CZI
treat="RNA-CTRL"
for cluster in C0 C1 C2 C3 C4 C5 ;do
for pcnum in $(seq 1 20); do
i=`ls -1 /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.bed | grep -v 'sort'`
  echo "running " $i
  #less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sortBed -i"}' > ${i%.*}.sort.bed
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sort -k1,1 -k2,2n "}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
done
done

module swap gnu7/7.3.0 gnu9

mkdir /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/permutations
mkdir /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/logs
mkdir /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/results
mkdir /rs/rs_grp_schold/CZI/RNA/analysis/fastQTL/logs/failedruns

######################################################################
###########################################################################

treat="CTRL"
method="voom"
#cluster="C0"
#i=10
for cluster in `awk 'NR>10{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
for i in $(seq 1 20); do
    njobs=`squeue -u fh8591 -r | wc -l`
    maxjobs=500
    if [ "$njobs" -gt "$maxjobs" ]; then
    echo waiting for jobspace
    sleep 1000
    else
for j in $(seq 1 30); do
#j=1
FILENAME=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations_10PCregress2step.chunk$j.txt.gz
#if [[ $(wc -l <$FILENAME) -ge 2 ]];then
#echo already run $cluster.$treat.PC$i.chunk$j
#else 
echo running $cluster
sbatch -q primary -N1-1 -n 2 --mem=15G -t 15000 --job-name=$cluster.$treat.PC$i.chunk$j -o /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/%j.out\
    --wrap "module load misc2; \
    fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$i.residuals_${method}_PCregress2step.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 30"
sleep 1
#fi
done
fi
done
echo finished submitting $cluster
#sleep 2000 #30ish min wait to try and not hit the max jobs limit -- may still be an issue
done

#identify jobs that did not finish due to node error
for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/*.out | awk '$6 == "Jun" && $7 >= 6 {print $9}'`; do 
j=`grep "Chunk processed" $i | sed 's/^[^0-9]*//' | sed 's/ .*//'1`
cluster=`grep "Scanning phenotype data" $i | cut -d . -f2`
treat=`grep "Scanning phenotype data" $i | cut -d . -f3`
cluster=`grep "Scanning phenotype data" $i | cut -d . -f2`
pcnum=`grep "Scanning phenotype data" $i | cut -d . -f4 | cut -d - -f2`
#   jobid=`echo $i | sed 's/.*\///' | sed 's/.out//'`
#    s=`sacct | grep "$jobid" | awk -v OFS='\t' '{print $2, $1}' | grep -v "batch" | grep -v "extern" | sort | uniq | tr '.' '\t' `
#    cluster=$(echo "$s" | awk '{print $1}')
#    treat=$(echo "$s" | awk '{print $2}')
#    pcnum=$(echo "$s" | awk '{print $3}' | sed 's/PC//')
#    j=$(echo "$s" | awk '{print $4}'| sed 's/chunk//')
    #echo $cluster $treat $pcnum $j
    if grep -q "Running time:" $i; then
#        #echo "file ran to completion"
        echo $cluster $treat $pcnum $j >> slurm/filesran.txt
        rm $i
    else
    #echo $cluster $treat PC$pcnum chunk$j within $i didnt finish
    echo $cluster $treat $pcnum $j >> slurm/needtorun_conditions.txt
#    mv $i slurm/oldrerun/
    fi 
done 
less slurm/filesran.txt | sort | uniq > slurm/filesran_unique.txt
less slurm/needtorun_conditions.txt | sort | uniq > slurm/needtorun_unique.txt

R
library(data.table)
f <- fread("slurm/needtorun_unique.txt",select=c(1:4))
c <- fread("slurm/filesran_unique.txt",select=c( 1:4))
notXY1 <- unique(merge(f,c,all.x = TRUE)[!merge(f,c)])
#f <- transform(f, combo=paste(V1,V2,V3,V4,sep="_"))
#c <- transform(c, combo=paste(V1,V2,V3,V4,sep="_"))
#fnotc <- subset(f, !combo %in% c$combo)
fwrite(notXY1, file="slurm/failedruns10_nocommpleted.txt", quote=F,sep='\t',row.names=F,col.names=F)
dim(notXY1)
q()
n
cat slurm/failedruns10_nocommpleted.txt | wc -l

#identify jobs that failed between time points
start="2025-04-28T14:00:00"
end="2025-04-28T23:59:59"
sacct | awk -v start=$start -v end=$end '$3 >= start && $3 <= end { print $0 }' | grep 'FAILED' | awk -v OFS='\t' '{print $2, $1}' | grep -v "batch" | grep -v "extern" | sort | uniq | tr '.' '\t' > slurm/failedruns.txt
sacct | awk -v start=$start -v end=$end '$3 >= start && $3 <= end { print $0 }' | grep 'COMPLETED' | awk -v OFS='\t' '{print $2, $1}' | grep -v "batch" | grep -v "extern" | sort | uniq | tr '.' '\t' > slurm/completedruns.txt

sacct --format="JobID,JobName%30,State"| awk -v OFS='\t' '{print $2, $1, $3}' | grep 'FAILED' | grep -v "batch" | grep -v "extern" | sort | uniq | tr '.' '\t' > slurm/failedruns.txt
sacct --format="JobID,JobName%30,State"| awk -v OFS='\t' '{print $2, $1, $3}' | grep 'COMPLETED' | grep -v "batch" | grep -v "extern" | sort | uniq | tr '.' '\t' > slurm/completedruns.txt

library(data.table)
f <- fread("slurm/failedruns.txt",select=c(1:4))
c <- fread("slurm/completedruns.txt",select=c( 1:4))
notXY1 <- unique(merge(f,c,all.x = TRUE)[!merge(f,c)])
notXY1 <- transform(notXY1, V4=gsub("chunk","",V4))
fwrite(notXY1, file="slurm/failedruns_nocommpleted.txt", quote=F,sep='\t',row.names=F,col.names=F)


#this reads in the failed runs and reruns fastqtl. I'm increasing time and memory here in case that was a problem
cat slurm/failedruns10_nocommpleted.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
#    oldjob=$(echo "$i" | awk '{print $4}')
    j=$(echo "$i" | awk '{print $4}')
#    j=$(echo "$i" | awk '{print $4}'| sed 's/chunk//')
    echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
    njobs=`squeue -u fh8591 | wc -l`
    maxjobs=200
    if [ "$njobs" -gt "$maxjobs" ]; then
    echo waiting for jobspace
    sleep 1000
    else
    #for j in $(seq 1 30); do --exclude=node\[117-118\]
    sbatch -q primary -N1-1 -n 5 --mem=25G -t 20000 --job-name=$cluster.$treat.PC$pcnum.chunk$j -o /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/%j.out\
        --wrap "module load misc2; \
        fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
        --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
        --permute 1000 10000 --window 1e6 \
        --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
        --chunk $j 30"
    #echo $i >> slurm/failedruns_ran.txt
    sleep 1
    #done
    fi
done

##############################3
#can run this while other jobs still going -- not fixed!!!
c="C0"
for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/*.out | awk '$6 == "Jun" && $7 >= 6 {print $9}'`; do 
    if grep -q "Default GSL error handler invoked" ; then #only setting jobs from previous submission batch that had an error
      slurmj=`grep "Chunk processed" $i | sed 's/^[^0-9]*//' | sed 's/ .*//'1`
      slurmcluster=`grep "Scanning phenotype data" $i | cut -d . -f2`
      slurmtreat=`grep "Scanning phenotype data" $i | cut -d . -f3`
      slurmpcnum=`grep "Scanning phenotype data" $i | cut -d . -f4 | cut -d - -f2`
grep $c slurm/failedruns10_nocommpleted.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
#    oldjob=$(echo "$i" | awk '{print $4}')
    j=$(echo "$i" | awk '{print $4}')
#    j=$(echo "$i" | awk '{print $4}'| sed 's/chunk//')
if [[ ${slurmj}${slurmcluster}${slurmtreat}${slurmpcnum} == ${j}${cluster}${treat}${pcnum} ]]; then

    echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
    njobs=`squeue -u fh8591 | wc -l`
    maxjobs=500
    if [ "$njobs" -gt "$maxjobs" ]; then
    echo waiting for jobspace
    sleep 1000
    else
    #for j in $(seq 1 30); do --exclude=node\[117-118\]
    sbatch -q primary -N1-1 -n 5 --mem=25G -t 20000 --job-name=$cluster.$treat.PC$pcnum.chunk$j -o /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/%j.out\
        --wrap "module load misc2; \
        fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
        --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
        --permute 1000 10000 --window 1e6 \
        --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
        --chunk $j 30"
    #echo $i >> slurm/failedruns_ran.txt
    sleep 1
    #done
    fi #end of max jobs
    fi #end of error grep 
    fi #end of match jobinfo
done
####################################


#rm fastQTL/output.txt
#exec > >(tee -a fastQTL/output.txt)
method="voom"
runnum=1
cat > fastQTL/endoffileerror${runnum}.txt

cat slurm/failedruns10_nocommpleted.txt | while read i || [[ -n $i ]];do #used for runnum=1
#cat fastQTL/endoffileerror.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
    j=$(echo "$i" | awk '{print $4}')
if ! zcat  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz
then
    echo $cluster $treat $pcnum $j >> fastQTL/endoffileerror${runnum}.txt
fi
done | gzip -c >  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.eQTL.txt.gz

#testing interactive
#/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/C0.CTRL.PC1-2.permutations_PCregress2step.chunk22.txt.gz

module load misc2
#cat slurm/failedruns3_nocommpleted.txt | while read i || [[ -n $i ]];do
#test C0
runnum=2
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
    grep "$cluster" /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/endoffileerror${runnum}.txt > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/endoffileerror${runnum}_${cluster}.txt
done
#for iline in `sed -n "${iline}p" fastQTL/endoffileerror_C1.txt`; do

#exec > >(tee -a fastQTL/endoffileerror${runnum}_${cluster}_output.txt)
exec > >(tee -a fastQTL/run${runnum}_output.txt)

cat fastQTL/endoffileerror${runnum}.txt | while read i || [[ -n $i ]];do
#cat fastQTL/endoffileerror${runnum}_${cluster}.txt | while read i || [[ -n $i ]];do
#cat fastQTL/endoffileerror.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
#    oldjob=$(echo "$i" | awk '{print $4}')
    j=$(echo "$i" | awk '{print $4}')
#    j=$(echo "$i" | awk '{print $4}'| sed 's/chunk//')
    echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
     fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
        --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
        --permute 1000 10000 --window 1e6 \
        --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
        --chunk $j 30
done

less fastQTL/run${runnum}_output.txt | grep "cluster C" -B5 > fastQTL/Cgenesfailed.txt
tail fastQTL/run${runnum}_output.txt >> fastQTL/Cgenesfailed.txt

#add ensg to line info to remove and rerun fastqtl
#ensg to remove:
#first check no 100% completed - C0 was 100% done one 1
ensgremovec0=`less fastQTL/C0genesfailed.txt | grep "ENSG" | tail -n1 | sed 's/.*\[//' | sed 's/\]//'`
ensgremovec2=`less fastQTL/C2genesfailed.txt | grep "ENSG" | sed 's/.*\[//' | sed 's/\]//'`
ensgremovec3=`less fastQTL/C3genesfailed.txt | grep "ENSG" | sed 's/.*\[//' | sed 's/\]//'`
ensgremovec4=`less fastQTL/C4genesfailed.txt | grep "ENSG" | sed 's/.*\[//' | sed 's/\]//'`
ensgremovec7=`less fastQTL/C7genesfailed.txt | grep "ENSG" | sed 's/.*\[//' | sed 's/\]//'`
printf "%s\n" "$ensgremovec0" "$ensgremovec2" "$ensgremovec3" "$ensgremovec4" "$ensgremovec7" > fastQTL/e5_generem.txt
paste fastQTL/endoffileerror5.txt fastQTL/e5_generem.txt > fastQTL/endoffileerror5generem.txt 

ensgremove=`less fastQTL/Cgenesfailed.txt | grep "ENSG" | sed 's/.*\[//' | sed 's/\]//'`


#less fastQTL/endoffileerror5generem.txt | egrep "C0|C2|C3|C4" > fastQTL/endoffileerror5generem_C0_c4.txt
exec > >(tee -a fastQTL/endoffileerror5C0_4_output.txt) #actually running all here 

cat fastQTL/endoffileerror5generem.txt  | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
    j=$(echo "$i" | awk '{print $4}')
    ensg=$(echo "$i" | awk '{print $5}')
    echo cluster $cluster treat $treat PCnum $pcnum chunk $j running without $ensg
    less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step.sort.bed.gz | grep -v $ensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step_rem$ensg.sort.bed.gz
    tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step_rem$ensg.sort.bed.gz 
     fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
        --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step_rem$ensg.sort.bed.gz  \
        --permute 1000 10000 --window 1e6 \
        --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_10PCregress2step.chunk$j.txt.gz \
        --chunk $j 30
done

#repeat this section as needed
newrun="3"
oldrun="2"
treat="CTRL"
less fastQTL/endoffileerror${oldrun}_output.txt | grep "cluster C" -B6 > fastQTL/Cgenesfailed.txt
tail fastQTL/endoffileerror${oldrun}_output.txt >> fastQTL/Cgenesfailed.txt
#open fastQTL/Cgenesfailed.txt and delete sections that completed

ensgremove=`less fastQTL/Cgenesfailed.txt | grep "\[ENSG" | sed 's/.*\[//' | sed 's/\]//'`
printf "%s\n" "$ensgremove" > fastQTL/e_generem.txt

rm fastQTL/endoffileerror${newrun}.txt
cat fastQTL/endoffileerror${oldrun}.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
    j=$(echo "$i" | awk '{print $4}')
if ! zcat  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_10PCregress2step.chunk$j.txt.gz
then
    echo $cluster $treat $pcnum $j >> fastQTL/endoffileerror${newrun}.txt
fi
done | gzip -c >  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_10PCregress2step.eQTL.txt.gz

#the orders don't match up always, need to rearrange fastQTL/endoffileerror${newrun}.txt to match fastQTL/Cgenesfailed.txt
paste fastQTL/endoffileerror${newrun}.txt fastQTL/e_generem.txt > fastQTL/endoffileerror${newrun}generem.txt 

R 
library(data.table)
newrun="7"
#oldrun=paste0(c("1"),collapse="")
oldrunlast=6
firstrun=FALSE
if(newrun==1){
    cat("no need to make a new file")
    } else if(newrun==2){
    ori <- fread(paste0("fastQTL/endoffileerror",oldrunlast,"generem.txt"),header=F)
    e6 <- fread(paste0("fastQTL/endoffileerror",newrun,"generem.txt"),header=F)
    allrem <- merge(ori, e6, by=c("V1"))
    fwrite(allrem, quote=F,row.names=F,col.names=F,sep='\t', file=paste0("fastQTL/endoffileerror",oldrunlast,newrun,"generem.txt"))
}else{
    ori <- fread(paste0("fastQTL/endoffileerror",oldrunlast-1,oldrunlast,"generem.txt"),header=F)
    e6 <- fread(paste0("fastQTL/endoffileerror",newrun,"generem.txt"),header=F)
    allrem <- merge(ori, e6, by=c("V1"))
    fwrite(allrem, quote=F,row.names=F,col.names=F,sep='\t', file=paste0("fastQTL/endoffileerror",oldrunlast,newrun,"generem.txt"))
}
q()
n

rm fastQTL/endoffileerror${newrun}_output.txt
exec > >(tee -a fastQTL/endoffileerror${newrun}_output.txt) #actually running all here 

cat fastQTL/endoffileerror${oldrun}${newrun}generem.txt  | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
    j=$(echo "$i" | awk '{print $4}')
    ensgs=$(echo "$i" | awk '{for(i=5;i<=NF;++i)print $i}')
    ensg=$(echo $ensgs | tr " " "|")
    ensgname=$(echo $ensgs | tr " " "_")
    echo cluster $cluster treat $treat PCnum $pcnum chunk $j running without $ensg
    less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step.sort.bed.gz | egrep -v $ensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step_rem$ensgname.sort.bed.gz
    tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step_rem$ensgname.sort.bed.gz 
     fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
        --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step_rem$ensgname.sort.bed.gz  \
        --permute 1000 10000 --window 1e6 \
        --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_10PCregress2step.chunk$j.txt.gz \
        --chunk $j 30
done

#mkdir /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/endoffileerror
#cat fastQTL/endoffileerror2.txt | while read i || [[ -n $i ]];do
#    cluster=$(echo "$i" | awk '{print $1}')
#    treat=$(echo "$i" | awk '{print $2}')
#    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
#    j=$(echo "$i" | awk '{print $4}')
#    mv /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_10PCregress2step.chunk$j.txt.gz /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/endoffileerror/$cluster.$treat.PC1-$pcnum.permutations_10PCregress2step.chunk$j.txt.gz
#done

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
    for i in $(seq 1 20); do
    echo combining cluster $cluster pcnum $i
        for j in $(seq 1 30); do
        zcat  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations_10PCregress2step.chunk*.txt.gz | gzip -c >  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations_10PCregress2step.eQTL.txt.gz
        done
    done
done

#low gene tested numbers, where is it coming from:
cluster="C0"
treat="CTRL"
clusters <- names(counts_ls)
treatments <- unique(metadata_ls[[1]]$treats)
best_df <- fread(paste0(outFolder,"bestPCs_table_PCregress2step.txt"))

df <- ldply(lapply(1:nrow(best_df),function(x){
info <- best_df[x,]
cluster=info$cluster
pcnum=info$PCs
    cat("running ",cluster,pcnum,"\n")
#just not gonna worry about the handful of removed genes here:
bed <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.",cluster,".",treat,".PC1-",pcnum,".residuals_qnorm_10PCregress2step.sort.bed.gz"))
permutations <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/",cluster,".",treat,".PC1-",pcnum,".permutations_10PCregress2step.eQTL.txt.gz"))
inout <- data.frame(inputgenes=dim(bed)[1],outputgenes=dim(permutations)[1])
return(inout)
}),data.frame)

cluster="C6"
countFile <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.",cluster,".",treat,".bed"))




library(data.table)
library(qqman)
library(qvalue)
library(ggplot2)
library(plyr)
library(tidyverse)

#cluster="C0"
treat="CTRL"
FDR <- 0.1
outFolder <- "/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/"
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.RData")
load(opfn)

bestPCtable <- lapply(names(counts_ls),function(cluster){
#  lapply(treatments,function(treat){
    cat("running ",cluster,treat,"\n")

myDir <- "/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/old_qnorm/" #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep("permutations_10PCregress2step.eQTL.txt.gz", filenames)] #pick specific files from list
filenames <- filenames[grepl(paste0(cluster,".",treat), filenames)] #pick specific files from list
data_names <- gsub(".permutations_10PCregress2step.eQTL.txt.gz", "", filenames) #remove file ending
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = F)[,analysis:=data_names[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
#combining data
all_PCs <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(all_PCs) <- data_names
# add colnames:
for(i in 1:length(all_PCs)){
    #cat("running",i,"\n")
    colnames(all_PCs[[i]]) <- c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval","analysis")
}
# add qvalue to each:
res <- data.frame("PCs"=integer(0),"eGenes"=integer(0))
for(i in 1:length(all_PCs)){
all_PCs[[i]]$bqval <- qvalue(all_PCs[[i]]$bpval)$qvalues
PC=gsub(".*-","",names(all_PCs)[i])
PC=gsub("[.].*","",PC)
res <- rbind(res, as.numeric(c(PC,sum(all_PCs[[i]]$bqval<FDR,na.rm =TRUE)))
)    }
colnames(res) <- c("PCs","eGenes")
best.PCs <- res[res$eGenes==max(res$eGenes),"PCs"]
best.index <- paste0(cluster,".",treat,".PC1-",best.PCs)
res <- res[order(res$PCs),]
if(length(best.PCs)>1){
    best.index <- paste0(cluster,".",treat,".PC1-",1)
    res <- res[order(res$PCs),]
    best.PCs <- best.PCs[1]
}
pc_results <- all_PCs[[1]]
ci=0.95
fwrite(res, file=paste0(outFolder,"results/",cluster,".",treat,".eGenes-per-GEPCs_PCregress2step.txt"), sep='\t', quote=F, row.names=F)

# save the best results:
pc_signif_pairs <- all_PCs[[best.index]]
fwrite(pc_signif_pairs, paste0(outFolder,"results/",cluster,".",treat,".FastQTL_results_best_", best.PCs, ".GEPCs_PCregress2step.txt"), sep='\t', quote=F, row.names=F)

# subset to significant only:
pc_signif_pairs <- pc_signif_pairs[pc_signif_pairs$bqval<FDR,]
pairs <- pc_signif_pairs[,c(1,6),] #geneid and snpid
fwrite(pairs, file=paste0(outFolder,"results/",cluster,".",treat,".PC",best.PCs,"_significant_topeeQTL_pairs_PCregress2step.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
##save SNP IDs
snps <- pc_signif_pairs[,"sid"]
fwrite(snps, file=paste0(outFolder,"results/",cluster,".",treat,".PC",best.PCs,"_significant_topeeQTL_snps_PCregress2step.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

pc_results_bp <- pc_results %>% select(pid, sid, bpval) %>% filter(!is.na(bpval)) %>%
            arrange(bpval) %>%
            mutate(r=rank(bpval, ties.method = "random"),
                   pexp=r/length(bpval),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(bpval)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(bpval)-r)))

pc_results_bp$group <- "Corrected p-value"
pc_results_bp <- pc_results_bp %>% dplyr::rename(pval = bpval) 

png(width = 12, height = 12, file=paste0(outFolder,"results/figures/",cluster,".",treat,".PC1-",best.PCs,"_eGene_qqplotpermuted_pvalue_only_PCregress2step.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
    p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
        ##    facet_grid(Origin ~ Location) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(best.index," eGene QQ Plot")) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()
    return(data.frame(cluster=cluster,PC=best.PCs))
})


#making QC tables
bestPCtableu <- ldply(bestPCtable, data.frame)
best_df <- ldply(lapply(names(counts_ls),function(c){
    cat("running", c, "\n")
    best.PCs <- subset(bestPCtableu, cluster==c)$PC
    pheno <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.",c,".",treat,".PC1-",best.PCs,".residuals_qnorm_10PCregress2step.sort.bed.gz"))
    pc_signif_pairs <- fread(paste0(outFolder,"results/",c,".",treat,".FastQTL_results_best_", best.PCs, ".GEPCs_PCregress2step.txt"))
    df <- data.frame(cluster=c,PCs=best.PCs,numInd=length(colnames(pheno))-4,testedgenes=dim(pc_signif_pairs)[1],eGenes_10=dim(pc_signif_pairs[pc_signif_pairs$bqval<0.1,])[1],eGenes_5=dim(pc_signif_pairs[pc_signif_pairs$bqval<0.05,])[1])
    return(df)
}),data.frame)
fwrite(best_df, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"bestPCs_table_PCregress2step.txt"))

#0 PCs run
module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0
for i in `ls -1 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.*.residuals_qnorm.bed | grep -v 'sort'`; do
  echo "running " $i
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sortBed -i"}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
done
module swap gnu7/7.3.0 gnu9

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
for j in $(seq 1 30); do
sbatch -q primary -N1-1 -n 2 --mem=12G -t 10000 --job-name=$cluster.$treat.PC0.chunk$j -o /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/%j.out \
    --wrap "module load misc2; \
    fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.residuals_qnorm.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.chunk$j.txt.gz \
    --chunk $j 30"
    sleep 1
done
done

for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/*.out | awk '$6 == "May" && $7 >= 16 {print $9}'`; do 
j=`grep "Chunk processed" $i | sed 's/^[^0-9]*//' | sed 's/ .*//'1`
cluster=`grep "Scanning phenotype data" $i | cut -d . -f2`
treat=`grep "Scanning phenotype data" $i | cut -d . -f3`
cluster=`grep "Scanning phenotype data" $i | cut -d . -f2`
    #echo $cluster $treat $pcnum $j
    if grep -q "Running time:" $i; then
        echo $cluster $treat $j >> slurm/filesranPC0.txt
        #rm $i
    else
    echo $cluster $treat $j >> slurm/needtorun_conditionsPC0.txt
    fi 
done 
#everything ran, good to continue -- no?
sacct --format="JobID,JobName%30,State"| awk -v OFS='\t' '{print $2, $1, $3}' | grep 'FAILED' | grep -v "batch" | grep -v "extern" | sort | uniq | tr '.' '\t' > slurm/failedruns.txt
sacct --format="JobID,JobName%30,State"| awk -v OFS='\t' '{print $2, $1, $3}' | grep 'COMPLETED' | grep -v "batch" | grep -v "extern" | sort | uniq | tr '.' '\t' > slurm/completedruns.txt

library(data.table)
f <- fread("slurm/failedruns.txt",select=c(1:4),header=F)
c <- fread("slurm/completedruns.txt",select=c( 1:4),header=F)
notXY1 <- unique(merge(f,c,all.x = TRUE)[!merge(f,c)])
notXY1 <- transform(notXY1, V4=gsub("chunk","",V4))
fwrite(notXY1, file="slurm/failedruns_nocommpleted.txt", quote=F,sep='\t',row.names=F,col.names=F)

cat slurm/failedruns_nocommpleted.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    j=$(echo "$i" | awk '{print $4}')
#    j=$(echo "$i" | awk '{print $4}'| sed 's/chunk//')
    echo cluster $cluster treat $treat chunk $j running
    njobs=`squeue -u fh8591 | wc -l`
    maxjobs=800
    if [ "$njobs" -gt "$maxjobs" ]; then
    echo waiting for jobspace
    sleep 1000
    else
    sbatch -q primary -N1-1 -n 5 --mem=25G -t 20000 --job-name=$cluster.$treat.chunk$j -o /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/%j.out\
        --wrap "module load misc2; \
    fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.residuals_qnorm.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.chunk$j.txt.gz \
    --chunk $j 30"
    sleep 1
    fi
done

rm fastQTL/endoffileerror.txt
cat slurm/failedruns_nocommpleted.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    j=$(echo "$i" | awk '{print $4}')
if ! zcat /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.chunk$j.txt.gz
then
    echo $cluster $treat $j >> fastQTL/endoffileerror.txt
fi
done | gzip -c >  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.eQTL.txt.gz

cat fastQTL/endoffileerror.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    j=$(echo "$i" | awk '{print $3}')
#    j=$(echo "$i" | awk '{print $4}'| sed 's/chunk//')
    echo cluster $cluster treat $treat chunk $j running
    njobs=`squeue -u fh8591 | wc -l`
    maxjobs=800
    if [ "$njobs" -gt "$maxjobs" ]; then
    echo waiting for jobspace
    sleep 1000
    else
    sbatch -q primary -N1-1 -n 5 --mem=25G -t 20000 --job-name=$cluster.$treat.chunk$j -o /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/%j.out\
        --wrap "module load misc2; \
    fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.residuals_qnorm.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.chunk$j.txt.gz \
    --chunk $j 30"
    sleep 1
    fi
done

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
echo running $cluster
for j in $(seq 1 30); do
     zcat  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.chunk$j.txt.gz
done | gzip -c >  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.eQTL.txt.gz;
done


library(data.table)
library(qqman)
library(qvalue)
library(ggplot2)
library(plyr)
library(tidyverse)

#cluster="C0"
treat="CTRL"
FDR <- 0.1
outFolder <- "/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/"
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
load(opfn)

lapply(names(counts_ls),function(cluster){
#  lapply(treatments,function(treat){
    cat("running ",cluster,treat,"\n")

all_PCs <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/",cluster,".",treat,".PC0.permutations.eQTL.txt.gz"))
colnames(all_PCs) <- c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
# add qvalue to each:
all_PCs$bqval <- qvalue(all_PCs$bpval)$qvalues
res <- sum(all_PCs$bqval<FDR,na.rm =TRUE)
ci=0.95

# save the best results:
best.PCs <- "PC0"
pc_signif_pairs <- all_PCs
fwrite(pc_signif_pairs, paste0(outFolder,"results/",cluster,".",treat,".FastQTL_results_best_", best.PCs, ".GEPCs.txt"), sep='\t', quote=F, row.names=F)

# subset to significant only:
pc_signif_pairs <- pc_signif_pairs[pc_signif_pairs$bqval<FDR,]
pairs <- pc_signif_pairs[,c(1,6),] #geneid and snpid
fwrite(pairs, file=paste0(outFolder,"results/",cluster,".",treat,".PC",best.PCs,"_significant_topeeQTL_pairs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
##save SNP IDs
snps <- pc_signif_pairs[,"sid"]
fwrite(snps, file=paste0(outFolder,"results/",cluster,".",treat,".PC",best.PCs,"_significant_topeeQTL_snps.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

pc_results_bp <- all_PCs %>% select(pid, sid, bpval) %>% filter(!is.na(bpval)) %>%
            arrange(bpval) %>%
            mutate(r=rank(bpval, ties.method = "random"),
                   pexp=r/length(bpval),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(bpval)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(bpval)-r)))

pc_results_bp$group <- "Corrected p-value"
pc_results_bp <- pc_results_bp %>% dplyr::rename(pval = bpval) 

png(width = 12, height = 12, file=paste0(outFolder,"results/figures/",cluster,".",treat,".",best.PCs,"_eGene_qqplotpermuted_pvalue_only.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
    p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
        ##    facet_grid(Origin ~ Location) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(best.PCs," eGene QQ Plot")) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()
})





########################################
############################################
cluster="C0"
treat="CTRL"
i=10
pc <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/covariates/",cluster,".",treat,".PC1-",i,".covariates-FastQTL.txt"))
pc0 <- fread("/wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_fastqtl/GEPCA/covariates/IBD_eQTL_Rectum-PC1-0.full.covariates-FastQTL.txt")
pc1 <- fread("/wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_fastqtl/GEPCA/covariates/IBD_eQTL_Rectum-PC1-1.full.covariates-FastQTL.txt")
pc10 <- fread("/wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_fastqtl/GEPCA/covariates/IBD_eQTL_Rectum-PC1-10.full.covariates-FastQTL.txt")
