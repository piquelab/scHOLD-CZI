for i in $(seq 1 30); do
for j in $(seq 1 30); do
    sbatch -q primary -N1-1 -n 2 --mem=12G -t 10000 --job-name=$i$j \
    --wrap "module load misc; \
    fastQTL --vcf /wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_bi-allelic_SNPs/IBD_eQTL_rectum_DNA_genotypes_filtered_SNPs.vcf.gz \
    --bed /rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/IBD-eQTL_Rectum_covariate_corrected_voom-q-norm.bed.gz \
    --permute 1000 10000 --window 1e5 --out output/PC1-$i.permutations.chunk$j.txt.gz \
    --cov /wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_fastqtl/GEPCA/covariates/IBD_eQTL_Rectum-PC1-$i.full.covariates-FastQTL.txt --chunk $j 30"
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
for (clus in clusters){
  for (i in gsub("-",".",treatments)){
    cat("running ",clus,i)
    colsclus <- colnames(all_counts_bed)[grepl(paste0("^",clus,"_"),colnames(all_counts_bed))]
    colstreat <- colsclus[grepl(paste0(i,"$"),colsclus)]
    cols <- grepl(paste0(colstreat,collapse="$|"),colnames(all_counts_bed))
    df_cols <- all_counts_bed[,..cols]
    colnames(df_cols) <- gsub("[.]","-",colnames(df_cols))
    colnames(df_cols) <- sapply(strsplit(colnames(df_cols),"_"),function(y)y[3])
    df <- unique(cbind(all_counts_bed[,c(1:4)],df_cols))
    df[is.na(df)] <- 0
    #write.table(df,sep='\t', quote=F, row.names=F, col.names=T, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.",clus,".",i,".bed"))
    fwrite(df, sep='\t', quote=F, row.names=F, col.names=T, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.",clus,".",i,".bed"))
    samples <- data.frame(samples=colnames(df[,-c(1:4)]))
    fwrite(samples, sep='\t', quote=F, row.names=F, col.names=F, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/sample_list_fastqtl.",clus,".",i,".txt"))
  }
}

#test <- table(colnames(df_cols))>1

#Sample IDs are specified in the header line. This line needs to start with a hash key (i.e. #).
module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0
for i in `ls -1 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.*.residuals_qnorm.bed | grep -v 'sort'`; do
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
u_eigenvec2 <- unique(fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt"))
#now testing subsetting the cov file to essential covariates
#u_eigenvec2 <- u_eigenvec2[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","pincme")]
u_eigenvec2 <- u_eigenvec2 %>% dplyr::select(-c(BATCH,Sex,cage1,Wave,genPC1,genPC2,genPC3,SCAIP1_6_genPC1,SCAIP1_6_genPC2,SCAIP1_6_genPC3,genPC1_pub,genPC2_pub,genPC3_pub,genPC1_old,genPC2_old,genPC3_old,csex1,wave_old))
design_expanded <- model.matrix(~0+ as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3) )

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
phe <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.C0.CTRL.residuals_qnorm.sort.bed.gz",header=T)
ind2 <- colnames(phe)[5:length(phe)]

## vcf file
ind3 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/sample_list_fastqtl.C0.CTRL.txt", header=F)$V1

identical(ind1, ind2)
identical(ind2, ind3)
identical(ind1,ind3)

not_identical <- ind1 != ind2
diff_indices <- which(not_identical)

fixind1 <- cv[c("variable",ind2)]
ind1 <- colnames(fixind1)[-1]
#fwrite(fixind1, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl_fixedindorder.txt")
fwrite(fixind1, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl_sub_fixedindorder.txt")

fixind1 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/covfile_fastqtl_fixedindorder.txt", header=T,data.table=F)
ind1 <- colnames(fixind1)[-1]
identical(ind1, ind2)
identical(ind1,ind3)

mkdir fastQTL
#moved files into here since it was more than I th0ought
mkdir fastQTL/covariates
mkdir fastQTL/nominal
mkdir fastQTL/permutations
mkdir fastQTL/results
mkdir fastQTL/results/figures

#had 1-30 PCs but unnecessary to do that many
treat="CTRL"
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
  for i in $(seq 1 20); do head -n $(($i+1)) /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/PCcovariates-FastQTL.$cluster.$treat.txt > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/covariates/$cluster.$treat.PC1-$i.covariates-FastQTL.txt; done
done
#cluster="C0"
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
#for cluster in `awk 'NR>1&&NR<7{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
#for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
echo running $cluster
for i in $(seq 1 20); do
for j in $(seq 1 30); do
#i=1
#j=1
FILENAME=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations.chunk$j.txt.gz
#FILESIZE=$(stat -c%s "$FILENAME")
if [[ $(wc -l <$FILENAME) -ge 2 ]];then
#if (( FILESIZE > 2)); then
echo already run $cluster.$treat.PC$i.chunk$j
else 
sbatch -q primary -N1-1 -n 2 --mem=12G -t 10000 --job-name=$cluster.$treat.PC$i.chunk$j \
    --wrap "module load misc2; \
    fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.$cluster.$treat.residuals_qnorm.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations.chunk$j.txt.gz \
    --cov /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/covariates/$cluster.$treat.PC1-$i.covariates-FastQTL.txt \
    --chunk $j 30"
sleep 1
fi
done
done
echo finished submitting $cluster
sleep 2000 #30ish min wait to try and not hit the max jobs limit -- may still be an issue
done

#combine output
#for cluster in `awk 'NR>1&&NR<7{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
echo running $cluster
for i in $(seq 1 20); do
for j in $(seq 1 30); do
     zcat  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations.chunk$j.txt.gz
done | gzip -c >  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations.eQTL.txt.gz;
done
done
#run fastqtl_chosebestPCs.R

#Also run no PCs (for interaction later)
treat="CTRL"
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
echo running $cluster
for j in $(seq 1 30); do
  sbatch -q primary -N1-1 -n 2 --mem=12G -t 10000 --job-name=$cluster.$treat.PC0.chunk$j \
    --wrap "module load misc2; \
   fastQTL \
      --permute 1000 10000 --window 1e6 \
      --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
      --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.$cluster.$treat.residuals_qnorm.sort.bed.gz \
      --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.chunk$j.txt.gz \
      --chunk ${j} 30"
      sleep 1
done
echo finished submitting $cluster
done
#combine output
#for cluster in `awk 'NR>1&&NR<7{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
echo running $cluster
for j in $(seq 1 30); do
     zcat  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.chunk$j.txt.gz
done | gzip -c >  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC0.permutations.eQTL.txt.gz;
done

#i=`less fastQTL/bestPCs_table.txt | awk -v pat="$cluster\$" '$1~pat{print $2}'` #this grabs the best PC number

#################
###################
#Also run no PCs (for interaction later)
   fastQTL \
      --vcf /wsu/home/groups/piquelab/SCAIP/vcf/SCAIP1-6_filtered_AF.vcf.gz \
      --bed ./1_normalized.data/${cell}_${treat}.bed.gz \
      --out ${outfn}.nominals.chunk${j}.txt.gz --window 1e6 --chunk ${j} 30








find /wsu/home/groups/piquelab/IBD_eQTL/ -name "IBD_Rectum.bed.gz"


#final run
mkdir fastQTL/bestPCout
for i in {1..30}; do
     sbatch -q express -p erprp --mem=10G --time=24:00:00 -N 1-1 -n 1 --job-name=$cluster.$treat.${i} --out=slurm_fastQTL_chunk.$cluster.$treat.${i}.output --wrap "
     module load misc;
     fastQTL --vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf.gz --bed IBD_Rectum.bed.gz --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/bestPCout/$cluster.$treat.nominals_chunk${i}.txt.gz --window 1e6 --chunk ${i} 30 \
     --cov /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/covariates/$cluster.$treat.PC1-8.covariates-FastQTL.txt "
     sleep 1;
done



#############################
##############################
for j in $(seq 1 30); do
     zcat  output/PC1-0.permutations.chunk$j.txt.gz
done | gzip -c >  output/PC1-0.permutations.eQTL.txt.gz




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
