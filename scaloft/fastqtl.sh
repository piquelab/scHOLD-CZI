for i in $(seq 1 30); do
for j in $(seq 1 30); do
    sbatch -q primary -N1-1 -n 2 --mem=12G -t 10000 --job-name=$i$j \
    --wrap "module load misc; \
    fastQTL --vcf /wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_bi-allelic_SNPs/IBD_eQTL_rectum_DNA_genotypes_filtered_SNPs.vcf.gz \
    --bed /rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/IBD-eQTL_Rectum_covariate_corrected_voom-q-norm.bed.gz \
    --permute 1000 10000 --window 1e5 --out output/PC1-$i.permutations.chunk$j.txt.gz \
    --cov ../GEPCA/covariates/IBD_eQTL_Rectum-PC1-$i.full.covariates-FastQTL.txt --chunk $j 30"
sleep 1
done;
done

run_FastQTL_threaded.py ${genotypes}.vcf.gz ${phenotypes}.bed.gz ${prefix} --covariates ${covariates}.txt.gz 
--permute 1000 10000 --window 1e6 --ma_sample_threshold 10 --maf_threshold 0.01 --chunks 100 --threads 10


/wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_bi-allelic_SNPs/IBD_eQTL_rectum_DNA_genotypes_filtered_SNPs.vcf.gz
#I think this should work
/rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.vcf.gz

#phenotypes file I think is the gene expression
#Chr    start   end     geneID      [ind1] [ind2] [ind3]
#bgzip myPhenotypes.bed && tabix -p bed phenotypes.bed.gz
#FastQTL requires for input VCF and BED files to contain the same samples, so it is necessary to do some filtering steps 
#before running the FastQTL tool. The filtering is done within the FastQTL Preprocessing tool which is based on the custom
#R script which retains only the intersection of samples between VCF and BED and moreover deals with GT values which are 
#different from 0/0, 0/1, 1/0 or 1/1 with .|. (FastQTL treats .|. as missing entries which are internally imputed as mean 
#dosage at the variant site).
#1. make sample list
#2. bcftools view -S sample_file.txt file.vcf > filtered.vcf
library(dplyr)
library(data.table)
library(tidyverse)
u_eigenvec2 <- u_eigenvec2 %>%
  select(Sample_ID, everything())
fwrite(u_eigenvec2, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt")

base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.RData")
load(opfn)
all_counts_ls <- lapply(names(counts_ls),function(cluster){
    #cluster=unique(all_metadata$letter_clusters)[1]
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    genes <- rownames(adjusted)
    adjusted <- cbind(genes, as.data.frame(adjusted))
    return(adjusted)
})
all_counts <- all_counts_ls %>% reduce(left_join, by = "genes")
samples <- data.frame(samples=unique(gsub("[.]","-",sapply(strsplit(colnames(all_counts),"_"),function(y)y[3]))[-1]))
fwrite(samples, sep='\t', quote=F, row.names=F, col.names=F, file="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/sample_list_fastqtl.txt")

module load bcftools/1.19
bcftools view --header-only /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.vcf.gz > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.txt
sed 's/RI\_//g' /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.txt > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.gsubRI.txt
bcftools reheader -h /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.gsubRI.txt /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.vcf.gz > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf

bcftools view -S /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/sample_list_fastqtl.txt /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf
bgzip /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf && tabix -p vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf.gz

#COMBAT counts to bed file
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
#Chr    start   end     ID      Samples...
all_counts_bed <- merge(geneIDs[,c("chr","start","end","symbol","ensgene")],all_counts,by.x="symbol",by.y="genes")
all_counts_bed <- all_counts_bed %>% relocate(ensgene,.after =end) %>% select(-symbol) # have to use ensgene as there was multi gene symbols
all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
fwrite(all_counts_bed, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.bed")

all_counts_bed <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.bed")
clusters <- names(counts_ls)
treatments <- unique(metadata_ls[[1]]$treats)
lapply(clusters,function(clus){
  lapply(gsub("-",".",treatments),function(i){
    cat("running ",clus,i)
    colsclus <- colnames(all_counts_bed)[grepl(clus,colnames(all_counts_bed))]
    colstreat <- colsclus[grepl(i,colsclus)]
    cols <- grepl(paste0(colstreat,collapse="|"),colnames(all_counts_bed))
    df_cols <- all_counts_bed[,..cols]
    colnames(df_cols) <- gsub("[.]","-",colnames(df_cols))
    colnames(df_cols) <- sapply(strsplit(colnames(df_cols),"_"),function(y)y[3])
    df <- unique(cbind(all_counts_bed[,c(1:4)],df_cols))
    df[is.na(df)] <- 0
    fwrite(df, sep='\t', quote=F, row.names=F, col.names=T, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.",clus,".",i,".bed"))
    samples <- data.frame(samples=colnames(df[,-c(1:4)]))
    fwrite(samples, sep='\t', quote=F, row.names=F, col.names=F, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/sample_list_fastqtl.",clus,".",i,".txt"))
  })
})

#Sample IDs are specified in the header line. This line needs to start with a hash key (i.e. #).
module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0
for i in `ls -1 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.*.bed | grep -v 'sort'`; do
  echo "running " $i
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sortBed -i"}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
done

#cov file
The file is TAB delimited
First row gives the sample ID and each additional one corresponds to a single covariate
First column gives the covariate ID and each additional one corresponds to a sample
The file should have S+1 rows and C+1 columns where S and C are the numbers of samples and covariates, respectively.

# transpose
u_eigenvec2 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt")
#now testing subsetting the cov file to essential covariates
u_eigenvec2 <- u_eigenvec2[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","pincme")]
t_eigenvec2 <- transpose(u_eigenvec2)
# get row and colnames in order
colnames(t_eigenvec2) <- u_eigenvec2$Sample_ID
rownames(t_eigenvec2) <- colnames(u_eigenvec2)
t_eigenvec2$variable <- colnames(u_eigenvec2)
t_eigenvec2 <- t_eigenvec2 %>%
  dplyr::select(variable, everything())
#fwrite(t_eigenvec2, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl.txt") #this one contains all columns
fwrite(t_eigenvec2, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl_sub.txt")


### covariates
#cv <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl.txt", header=T,data.table=F)
cv <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl_sub.txt", header=T,data.table=F)
ind1 <- colnames(cv)[-1]

### bed file phenotypes_fastqtl.C9.PHA.fh.sort.bed.gz
phe <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.C0.CTRL.sort.bed.gz",header=T)
ind2 <- colnames(phe)[5:length(phe)]

## vcf file
ind3 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/sample_list_fastqtl.C0.CTRL.txt", header=F)$V1

identical(ind1, ind2)
identical(ind2, ind3)
identical(ind1,ind3)

not_identical <- ind1 != ind2
diff_indices <- which(not_identical)

fixind1 <- cv[c("variable",ind2)]
ind1 <- colnames(fixind1)[-1]
#fwrite(fixind1, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl_fixedindorder.txt")
fwrite(fixind1, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl_sub_fixedindorder.txt")

mkdir fastQTL
#moved files into here since it was more than I th0ought
#for i in $(seq 1 30); do
#for j in $(seq 1 30); do
j=1
cluster="C0"
treat="CTRL"
    sbatch -q primary -N1-1 -n 2 --mem=12G -t 10000 --job-name=$cluster.$treat.chunk$j \
    --wrap "module load misc2; \
    fastQTL --vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl_sub.$cluster.$treat.sort.bed.gz \
    --permute 1000 10000 --window 1e5 --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/$cluster.$treat.permutations.chunk$j.txt.gz \
    --cov /rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl_sub_fixedindorder.txt \
    --chunk $j 300"
#sleep 1
#done

prob <- phe[phe$ensgene=="ENSG00000188976",]
png(width = 12, height = 12, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/ENSG00000188976.hist.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
hist(unlist(prob[,-c(1:4)]))
dev.off()
#not obviously super low expression, just testing removing this 
probgenes <- c("ENSG00000188976","ENSG00000187961","ENSG00000187583")
fixphe <- subset(phe, !ensgene %in% probgenes)
fwrite(fixphe, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl_sub.C0.CTRL.sort.bed")
less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl_sub.C0.CTRL.sort.bed | awk -v OFS='\t' 'NR == 1{print "#"$0;next}; NR > 1 {print $0}' | bgzip >/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl_sub.C0.CTRL.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl_sub.C0.CTRL.sort.bed.gz
