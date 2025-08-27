#save residuals for built in tensorqtl (treatment)

library(tidyr)
library(tidyverse) #gave me issues 031725 -- hax node
library(glmnet)
library(data.table)
library(methods)
library(plyr)
library(ggplot2)
library(ggpubr)
library(irlba)

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

base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
resmethod="voom"
resFolder <- paste0(base,"residuals/")
if (!file.exists(paste0(resFolder,"counts/"))) dir.create(paste0(resFolder,"counts/"), showWarnings=F)
contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA-DEX")) #no lps vs lps-dex as too few ind
x=2
con=contrastdf[x,]
contrast <- paste0(con$treatment,"_vs_",con$control)
opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)
data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
filenames_cov <- list.files(paste0(data_path,"residuals/covariates/")) #file list from directory
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
best_df <- fread(file=paste0(outFolder,"bestPCs_table.txt"))

allsamples <- ldply(lapply(names(counts_ls),function(clus){
  data_control <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/",clus,".",con$control,".residuals_",resmethod,".txt"))
  data_treat <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/",clus,".",con$treatment,".residuals_",resmethod,".txt"))

  colnames(data_control) <- gsub("[.]","-",colnames(data_control))
  colnames(data_treat) <- gsub("[.]","-",colnames(data_treat))

  names(data_control)[c(2:length(colnames(data_control)))] <- paste0(colnames(data_control)[-c(1)],"-",con$control)
  names(data_treat)[c(2:length(colnames(data_treat)))] <- paste0(colnames(data_treat)[-c(1)],"-",con$treatment)
  inds_ctrl <- colnames(data_control)[-c(1)]
  inds_treat <- colnames(data_treat)[-c(1)]

  datam <- merge(data_control,data_treat, by=c("V1"))
  cv_data <- data.frame(Sample_ID=c(inds_ctrl,inds_treat),treats=c(rep(con$control,length(inds_ctrl)),rep(con$treatment,length(inds_treat))))

  #chr | start | end | gene id | ind1 | ind2 ...
  all_counts_bed <- merge(geneIDs[,c("chr","strand_start","strand_end","symbol","ensgene")],datam,by.x="ensgene",by.y="V1")
  all_counts_bed <- all_counts_bed %>% relocate(ensgene,.after =strand_end) %>% dplyr::select(-symbol) # have to use ensgene as there was multi gene symbols
  #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
  names(all_counts_bed)[c(2:4)] <- c("start","end","gene_id")
  #fwrite(all_counts_bed, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"counts/phenotypes.",clus,".",contrast,".bed"))

  best.PCs <- subset(best_df, cluster==clus)$PCs
  covar <- fread(paste0(data_path,"residuals/covariates/",filenames_cov[grep(paste0(clus,".",con$treatment,".PC1-",best.PCs,".covariates_voom.txt"),filenames_cov)]))
  names(covar)[c(2:length(colnames(covar)))] <- paste0(colnames(covar)[-c(1)],"-",con$control)
  pc_signif_pairs <- subset(fread(paste0(outFolder,"results/",clus,".",treat,".best_", best.PCs, ".GEPCs.txt")),qval<0.1)
  phenosub <- subset(all_counts_bed, gene_id %in% pc_signif_pairs$phenotype_id)
  
  vcfind <- colnames(fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/","vcf/ref.ac1.",clus,".",con$treatment,".vcfheader.txt")))
  # run PCA on residuals
  PCs <- prcomp_irlba(t(data_treat[,-1]), n=best.PCs)
  summary(PCs)
  mypcs <- as.data.frame(PCs$x)
  rownames(mypcs) <-inds_treat

  covs <- as_tibble(t(mypcs))
  #order to vcf
  #covsord <- covs[,match(vcfind, colnames(covs)) ]
  cvs <- cbind(id=colnames(mypcs),covs)
#    fwrite(cvs, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"covariates/PCcovariates_",resmethod,".",clus,".",con$treatment,".txt"))

  bothcovar <- cbind(covar,cvs[,-1])

  phenoind <- colnames(phenosub[,-c(1:4)])
  covarind <- colnames(bothcovar)[-1]
  common_ind <- intersect(phenoind, covarind)

  bothcovarsubcov <- subset(bothcovar, select = c(colnames(covar)[1],common_ind))
  phenoindsubcov <- subset(phenosub, select = c(colnames(phenosub)[c(1:4)],common_ind))
  cv_datasubcov <- subset(cv_data, Sample_ID %in% common_ind)
  cv_datasubcov <- transform(cv_datasubcov, treats=ifelse(treats==con$control,0,1))

  covvarv <- cv_datasubcov[,-1, drop=F]
  rownames(covvarv) <- cv_datasubcov$Sample_ID #NEED ROWNAMES
  names(covvarv)[1] <- contrast
  fwrite(covvarv,file=paste0(data_path,"residuals/covariates/",clus,".",contrast,".for_tensorqtl_int.txt"),sep="\t",col.names=T,row.names=T, quote=F)
  fwrite(bothcovarsubcov,paste0(data_path,"residuals/covariates/",clus,".",contrast,".PC1-",best.PCs,".covariates_voom.txt"),sep="\t",col.names=T,row.names=F, quote=F)
  fwrite(phenoindsubcov,paste0(data_path,"residuals/","phenotypes.",clus,".",contrast,".residuals_voom.eQTLonly.bed"),sep="\t",col.names=T,row.names=F, quote=F)
  samples <- data.frame(cluster=clus,samples=common_ind)
  return(samples)
}), data.frame)

allsamplesu <- unique(allsamples[,-1,drop=F])
fwrite(allsamplesu, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(resFolder,"vcf/sample_list.",contrast,".txt"))

less /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.gsubRI.txt | grep -v "^##" > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.sampleheaderonly.gsubRI.txt

vcfheader <- fread("/rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.sampleheaderonly.gsubRI.txt")
dat <- select(vcfheader,starts_with("AL-"))
ctrlvcfind <- data.frame(samples=paste0(colnames(dat),"-",con$control))
treatvcfind <- data.frame(samples=paste0(colnames(dat),"-",con$treatment))
fwrite(ctrlvcfind, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(resFolder,"vcf/sample_list.",con$control,".txt"))
fwrite(treatvcfind, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(resFolder,"vcf/sample_list.",con$treatment,".txt"))

module load bcftools/1.19
control="CTRL"
treatment="PHA"
bcftools reheader --samples /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/sample_list.$control.txt -o /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$control.vcf.gz /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf
bcftools reheader --samples /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/sample_list.$treatment.txt -o /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.vcf.gz /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf
tabix -p vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$control.vcf.gz
tabix -p vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.vcf.gz
bcftools merge /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$control.vcf.gz /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.vcf.gz -o /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.$control.vcf.gz

bgzip /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.$control.vcf.gz && tabix -p vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.$control.vcf.gz

plink2 --vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink/ref.ac1

#add # to header and sort and gzip
module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0
method="voom"
#ALOFT
var="PHA_vs_CTRL"
for cluster in C0 C1 C10 C11 C2 C3 C4 C5 C6 C7 C8 C9;do
  i=`ls -1 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/phenotypes.$cluster.$var.residuals_${method}.eQTLonly.bed | grep -v 'sort'`
  echo "running " $i
  #less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sortBed -i"}' > ${i%.*}.sort.bed
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sort -k1,1 -k2,2n "}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
done
module swap gnu7/7.3.0 gnu9


conda init bash
source ~/.bashrc
module swap gnu7/7.3.0 gnu9

conda activate tensorqtl_p3.11_env
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
module load R/4.3.2 #old: module load r/4.2.0

data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output"
mkdir -p ${outFolder}/interaction
out_path="${outFolder}/interaction"
var="PHA_vs_CTRL"
cluster="C0"
PC=`grep $cluster ${outFolder}/bestPCs_table.txt | cut -f2`
geno="${data_path}/residuals/vcf_plink/ref.ac1" # prefix for plink triplet files
out_prefix="${cluster}_${var}_tensorqtlint" # prefix for output file
pheno=${data_path}/residuals/phenotypes.$cluster.${var}.residuals_voom.eQTLonly.sort.bed.gz
covar="${data_path}/residuals/covariates/${cluster}.${var}.PC1-$PC.covariates_voom.txt"
int="${data_path}/residuals/covariates/${cluster}.${var}.for_tensorqtl_int.txt"

if [ -s "${out_path}/${out_prefix}.cis_qtl_top_assoc.txt.gz" ] ; then
    echo "[$(date)] Output already exists. Skipping..."
    exit 0
fi

### -------- LOGGING -------- ###
echo "[$(date)] Starting tensorQTL interaction on node: $(hostname)"
echo "Contrast=$var, CLUSTER=$cluster, PC=$PC"
echo "pheno file: $pheno"
echo "geno file: $geno"
echo "covar file: $covar"
echo "interactions file: $int"
echo "Output file: ${out_path}/${out_prefix}"
echo "settings: mode=cis_nominal;best_only;fdr=0.1"

python3 -m tensorqtl ${geno} ${pheno} ${out_prefix} \
    --covariates ${covar} \
    --interaction ${int} \
    --best_only \
    --fdr 0.1 \
    --mode cis_nominal \
    -o ${out_path}

echo end on "[$(date)]"