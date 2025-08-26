#jaxqtl from https://github.com/mancusolab/jaxQTL https://www.medrxiv.org/content/10.1101/2025.01.18.25320755v1

git clone https://github.com/mancusolab/jaxqtl.git
#jaxqtl directory : /wsu/home/groups/piquelab/cindy/jaxqtl

conda create -n jaxqtl python=3.10.9
#warning: newer version of conda available: 
#current= 22.9.0
#newest = 25.5.1
conda init bash
source ~/.bashrc
conda activate jaxqtl
cd /wsu/home/groups/piquelab/cindy/jaxqtl
pip install -e .
pip install lineax
pip install qtl

mkdir /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/logs

#example
example_data_path="./tutorial/input"
example_out_path="./tutorial/output"

celltype="CD4_NC"

# genelist to perform cis-eQTL mapping
example_chr=22
example_chunk_file="genelist_10"

# choose test method: score test (recommended) or wald
test_method="score"

# choose cis or nominal scan
mode="cis"
example_window=500000 # default extend 500kb on either side, i.e., [start-window, end+window]

# jaxQTL by default compute expression PCs using the entire data provided in *.bed.gz
# to disable adding expression PCs, set this to 0
num_expression_pc=2

example_pheno="${example_data_path}/${celltype}.N100.bed.gz"
example_geno="${example_data_path}/chr${example_chr}" # prefix for plink triplet files
example_covar="${example_data_path}/donor_features.tsv"

# choose gene list for eQTL mapping
example_genelist="${example_data_path}/${example_chunk_file}"

# choose eQTL model: NB for negative binomial, poisson, gaussian
model="NB"

# if using permutation method to calibrate gene-level p value, set number of permutation
nperm=1000

# choose platform: cpu, gpu, tpu
platform="cpu"

# prefix for output file
example_out="${example_out_path}/${celltype}_chr${example_chr}_${example_chunk_file}_jaxqtl_${model}"
jaxqtl \
 --geno ${example_geno} \
 --covar ${example_covar} \
 --pheno ${example_pheno} \
 --model ${model} \
 --mode ${mode} \
 --window ${example_window} \
 --genelist ${example_genelist} \
 --test-method ${test_method} \
 --nperm ${nperm} \
 --addpc ${num_expression_pc} \
 --standardize \
 -p ${platform} \
 --out ${example_out}
#gives tpu warning, but this is safe to ignore


#Phenotypes are provided in BED format, with a single header line starting with # and the first four columns corresponding to: 
#chr, start, end, phenotype_id, with the remaining columns corresponding to samples (the identifiers must match those in the genotype input). 
#The BED file can specify the center of the cis-window (usually the TSS), with start == end-1, or alternatively, start and end positions, in which case the cis-window is [start-window, end+window]
#actually appear to be straight counts, not residuals
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
load(opfn)
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
treatments=c("CTRL","LPS","PHA","PHA-DEX")
} else if(job=="CZI"){
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",11,0.2) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
figuredir=paste0(outFolder,"figures/")
combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
load(opfn)
counts_ls$C6 <- NULL
metadata_ls$C6 <- NULL
treatments=c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")
}

outFolder=paste0(base,"jaxQTL/")
if (!file.exists(paste0(outFolder,"vcf/"))) dir.create(paste0(outFolder,"vcf/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"counts/"))) dir.create(paste0(outFolder,"counts/"), showWarnings=F)

cv <- unique(u_eigenvec2)
##check for individuals with more than one sample
cv %>% count("Sample_ID") %>% dplyr::filter(freq>1)

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
geneIDs <- subset(geneIDs, chr %in% c(1:22))

i="CTRL"
clus="C4"
lapply(names(counts_ls),function(clus){
  #for (i in treatments){
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
    names(all_counts_bed)[c(1:4)] <- c("Chr","start","end","Geneid")
    #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr)) #when converting vcf to plink, the chr bit is removed so don't want it here
    fwrite(all_counts_bed, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"counts/phenotypes.",clus,".",i,".bed"))
    samples <- data.frame(samples=colnames(all_counts_bed[,-c(1:4)]))
    fwrite(samples, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(outFolder,"vcf/sample_list.",clus,".",i,".txt"))
      cat(clus,i, length(samples$sample),"\n")
  #}
})

#now need to format and sort phenotype file
module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0

for i in /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/counts/phenotypes.*.bed ; do
#for i in /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/counts/phenotypes.C0.CTRL.bed ; do
  if [ ! -f "${i%.*}.sort.bed.gz" ]; then
echo sorting and zipping $i
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sort -k1,1 -k2,2n "}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
fi
done

module swap gnu7/7.3.0 gnu9

#want to create chunks to run (200-300 genes recommended per chunk) in A single-column (no header) file specifying gene identifiers
R
library(data.table)
library(plyr);library(dplyr)
#job="CZI"
job="ALOFT"
#load in original counts info for clusters and treats used
if(job=="ALOFT"){
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
#cov_file=fread("/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt")
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
# load annotation
u_eigenvec2 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt")
u_eigenvec2 <- transform(u_eigenvec2, Wavenumeric=as.numeric(as.factor(Wave))) #poossibly can't have non-numeric?
#u_eigenvec2 <- transform(u_eigenvec2, Wavenonum=gsub('[[:digit:]]+', '', Wave)) #poossibly can't have non-numeric?
opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)
basic_covs <- unique(u_eigenvec2[,c("Sample_ID","Sex","cage1","Wavenumeric","genPC1","genPC2","genPC3")])
basic_covs %>% count(Sample_ID) %>% dplyr::filter(n>1)
} else if(job=="CZI"){
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",11,0.2) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
outFolder=paste0(base,"jaxQTL/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
leadvar <- fread("/rs/rs_grp_schold/covariates/other_covariates/HOLD LEAD 5.27.25.csv")
cov_pluslead <- merge(cov_file,leadvar,by="pID",all=T)
cov_pluslead <- transform(cov_pluslead, Lead=ifelse(Lead==-99,NA,Lead))
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o[,-c("sex","sex_alph","age")],cov_pluslead,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
load(opfn)
basic_covs <- unique(eigenvec2[,c("Sample_ID","sex_alph","age","PC1","PC2")])
basic_covs %>% count(Sample_ID) %>% dplyr::filter(n>1)
}

clus="C0"
treat="CTRL"
n=50 #300 is too big
if (!file.exists(paste0(outFolder,"gene_chunks_50/"))) dir.create(paste0(outFolder,"gene_chunks_50/"), showWarnings=F)
lapply(names(counts_ls),function(clus){
pheno <- fread(paste0(outFolder,"counts/phenotypes.",clus,".",treat,".sort.bed.gz"))
#split into chromosomes
lapply(split(pheno,pheno[,1]),function(chr){
	chrname=unique(unlist(chr[1,1]))
	genes <- chr$Geneid
	splitlist <- split(genes, ceiling(seq_along(genes)/n))
	counter=1
	lapply(splitlist, function(i) {
		cat(counter,"\t")
		genedf <- data.frame(ensgene=i)
		fwrite(genedf, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(outFolder,"gene_chunks_50/chr",chrname,".gene_list_",counter,"_",clus,"_",treat,".txt"))
		counter <<- counter + 1
	})
})
})

#Covariates are provided as a tab-delimited tsv file dataframe (samples x covariates) with column headers.
#my covariates are sex, age, wave, genpcs. expression PCs are computed by jaxqtl

if (!file.exists(paste0(outFolder,"covariates/"))) dir.create(paste0(outFolder,"covariates/"), showWarnings=F)
fwrite(basic_covs, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"covariates/basic_covs.tsv"))

#Genotypes can be provided in PLINK1 bed/bim/fam format. (We will accomodate other formats in the future verssions)
module load plink/2.0
cluster="C0"
treat="CTRL"
plink2 --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz --make-bed --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat
#Start time: Mon Jul  7 16:59:31 2025
#End time: Mon Jul  7 19:26:41 2025
chr=1
#note that this removes "chr" from chr nomenclature
#for chr in {1..22};do
#  echo running chr $chr
#  plink2 --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz --chr chr${chr} --make-bed --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/vcf/ref.ac1.$cluster.$treat.$chr
#done

data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL"
vcf_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf"
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt | grep -v "C4"`;do
for chr in {1..22};do
  echo submitting chr $chr cluster $cluster
  sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",data_path="$data_path",vcf_path="$vcf_path" ${data_path}/src/vcf_to_plink.sh 
  sleep 1
done
done

mkdir -p /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output

data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output"

cluster="C0"

# genelist to perform cis-eQTL mapping
chr=1
n_chunks=`ls ${data_path}/gene_chunks/${chr}.* | wc -l`
for j in $(seq 1 $n_chunks);do
chunk_file="${chr}.gene_list_${j}"

# choose test method: score test (recommended) or wald
test_method="score"

# choose cis or nominal scan
mode="cis"
window=1000000 # default extend 500kb on either side, i.e., [start-window, end+window] ##CANT USE 1e6

# jaxQTL by default compute expression PCs using the entire data provided in *.bed.gz
# to disable adding expression PCs, set this to 0
num_expression_pc=2

pheno="${data_path}/counts/phenotypes.${cluster}.${treat}.sort.bed.gz"
geno="${data_path}/vcf/ref.ac1.${cluster}.${treat}.${chr}" # prefix for plink triplet files
covar="${data_path}/covariates/basic_covs.tsv"

# choose gene list for eQTL mapping
genelist="${data_path}/gene_chunks/chr${chunk_file}_${cluster}_${treat}.txt"

# choose eQTL model: NB for negative binomial, poisson, gaussian
model="NB"

# if using permutation method to calibrate gene-level p value, set number of permutation
nperm=1000

# choose platform: cpu, gpu, tpu
platform="cpu"

# prefix for output file
out="${out_path}/${cluster}_${treat}_${chunk_file}_jaxqtl_${model}"

jaxqtl \
 --geno ${geno} \
 --covar ${covar} \
 --pheno ${pheno} \
 --model ${model} \
 --mode ${mode} \
 --window ${window} \
 --genelist ${genelist} \
 --test-method ${test_method} \
 --nperm ${nperm} \
 --addpc ${num_expression_pc} \
 --standardize \
 -p ${platform} \
 --out ${out}

done

cluster="C0"
PC=2
normmethod="voom"
chr=1
treat="CTRL"

data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output"
mkdir -p /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output

ls ${data_path}/gene_chunks_50/chr${chr}.*_${cluster}_${treat}.txt


for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt | grep -v "C4"`;do
for chr in {1..22};do
  n_chunks=`ls ${data_path}/gene_chunks_50/chr${chr}.*_${cluster}_${treat}.txt | wc -l`
  for PC in {2..10}; do
  njobs=`squeue -u fh8591 -r| wc -l`
  maxjobs=500
  while [ "$njobs" -gt "$maxjobs" ];do 
  echo waiting for jobspace, sleeping ...
  sleep 300 
  njobs=`squeue -u fh8591 -r| wc -l`
  done
  for j in $(seq 1 $n_chunks);do
    echo submitting cluster $cluster chr $chr chunk $j PC $PC
    sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",PC="$PC",j="$j",data_path="$data_path",out_path="$out_path" ${data_path}/src/run_jaxqtl.sh 
    sleep 1
  done
  done
done
done

treat="CTRL"
for i in ${data_path}/gene_chunks_50/chr*.gene_list*_$treat.txt; do
  chr=`echo "$i" |cut -d"/" -f8 | cut -d"." -f1 | sed 's/chr//g'`
  j=`echo "$i" |cut -d"/" -f8 | cut -d"_" -f3 `
  cluster=`echo "$i" |cut -d"/" -f8 | cut -d"_" -f4 `
  njobs=`squeue -u fh8591 -r| wc -l`
  maxjobs=500
  while [ "$njobs" -gt "$maxjobs" ];do 
  echo waiting for jobspace, sleeping ...
  sleep 300 
  njobs=`squeue -u fh8591 -r| wc -l`
  done #while end
  for PC in {2..10}; do
    chunk_file="${chr}.gene_list_${j}"
    out="${out_path}/${cluster}_${treat}_${chr}.gene_list_${j}_jaxqtl_${model}_PC${PC}.cis_score.tsv.gz"
    if [ -f "$out" ] ; then
      echo " Output already exists. Skipping..."
    else
      echo submitting cluster $cluster chr $chr chunk $j PC $PC
      sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",PC="$PC",j="$j",data_path="$data_path",out_path="$out_path" ${data_path}/src/run_jaxqtl.sh 
      sleep 1
    fi #if not file exists
  done #PCnum end
done #gene list combo end

#check what has run
treat="CTRL"
model="NB"
#2868
rm ${out_path}/not_completed.txt ${out_path}/completed.txt
rm ${out_path}/not_completed.txt ${out_path}/not_completed.txt
for i in ${data_path}/gene_chunks_50/chr*.gene_list*_$treat.txt; do
  chr=`echo "$i" |cut -d"/" -f8 | cut -d"." -f1 | sed 's/chr//g'`
  j=`echo "$i" |cut -d"/" -f8 | cut -d"_" -f3 `
  cluster=`echo "$i" |cut -d"/" -f8 | cut -d"_" -f4 `
  for PC in {2..10}; do
    chunk_file="${chr}.gene_list_${j}"
    out="${out_path}/${cluster}_${treat}_${chr}.gene_list_${j}_jaxqtl_${model}_PC${PC}.cis_score.tsv.gz"
    #echo $out
    if [ -f "$out" ] ; then
      #echo " Output already exists."
      echo "${cluster}_${treat}_${chr}_${j}_${PC}" >> ${out_path}/completed.txt
    else
      echo $out not run
      echo "${cluster}_${treat}_${chr}_${j}_${PC}" >> ${out_path}/not_completed.txt
    fi
done
    #sleep 1
done
#25812 - 10305
less ${out_path}/completed.txt | wc -l #16403 23900
less ${out_path}/not_completed.txt | wc -l #9409 1912

cluster="C6"
chr=8
j=1
PC=2
sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",PC="$PC",j="$j",data_path="$data_path",out_path="$out_path" ${data_path}/src/run_jaxqtl.sh 
chr=9
j=4
PC=7
sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",PC="$PC",j="$j",data_path="$data_path",out_path="$out_path" ${data_path}/src/run_jaxqtl.sh 



for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt| grep -v "C4"`;do
    for i in $(seq 2 10); do
    echo combining cluster $cluster PC $i $treat
        zcat $out_path/${cluster}_${treat}_*.gene_list_*_jaxqtl_NB_PC${i}.cis_score.tsv.gz | gzip -c > $out_path/${cluster}_${treat}.jaxqtl_NB_PC${i}.cis_score.tsv.gz
    done
done

for i in /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/logs/*.out; do 
if grep -q "genelist: /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/gene_chunks_50/chr5.gene_list_3_C1_CTRL.txt" $i; then
echo $i
fi
done

i=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/gene_chunks_50/chr5.gene_list_3_C1_CTRL.txt

/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/C0_CTRL_1.gene_list_1_jaxqtl_NB.cis_score.tsv.gz

mkdir -p $out_path/logs
mv $out_path/*.log $out_path/logs

R
library(data.table)
library(qvalue)
library(ggplot2)
library(plyr)
library(tidyverse)

cluster="C0"
treat="CTRL"
chr="1"
PC=2
ci=0.95

if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

#phenotype_id    chrom   num_var variant_id      pos     tss_distance    ma_count        af      beta_shape1     beta_shape2     beta_converged  opt_status      true_nc pval_nominal    slope   slope_se        pval_beta       alpha_cov       model_converged
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/"
#jaxqtlouttest <- fread(paste0(outFolder,cluster,"_",treat,"_",chr,".gene_list_",j,"_jaxqtl_NB.cis_score.tsv.gz"))
jaxqtlout <- fread(paste0(outFolder,cluster,"_",treat,"_",chr,".jaxqtl_NB.cis_score.tsv.gz"))
jaxqtlout <- transform(jaxqtlout, model_converged=as.numeric(model_converged),beta_converged=as.numeric(beta_converged),
	pval_beta=as.numeric(pval_beta), slope=as.numeric(slope),slope_se=as.numeric(slope_se), pval_nominal=as.numeric(pval_nominal))
converged <- subset(jaxqtlout,model_converged > 0 & beta_converged > 0)
converged <- transform(converged, padj=p.adjust(pval_beta, method="BH"), bqval=qvalue(pval_beta)$qvalues)
table(converged$padj<0.1)
#87
table(converged$bqval<0.1)

png(width = 12, height = 12, file=paste0(outFolder,"figures/",cluster,".",treat,".",chr,"_eGene_qqplotpermuted_pvalue.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
    p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(chr," eGene QQ Plot")) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()

converged$pval_beta[converged$pval_beta<1e-20] <- 1e-20
pc_results_bp <- converged %>% select(phenotype_id, variant_id, pval_beta) %>% filter(!is.na(pval_beta)) %>%
            arrange(pval_beta) %>%
            mutate(r=rank(pval_beta, ties.method = "random"),
                   pexp=r/length(pval_beta),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)))

pc_results_bp$group <- "Corrected p-value"
pc_results_bp <- pc_results_bp %>% dplyr::rename(pval = pval_beta) 

png(width = 12, height = 12, file=paste0(outFolder,"figures/",cluster,".",treat,".",chr,"_eGene_qqplotpermuted_pvalue_clipped.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
    p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(chr," eGene QQ Plot")) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()



