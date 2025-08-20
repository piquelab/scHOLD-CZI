module load plink/2.0
mkdir -p /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink
cluster="C3"
treat="CTRL"
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt `;do
for chr in {1..22};do
  if [ ! -f /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink/ref.ac1.$cluster.$treat.${chr}.psam ];then
echo running cluster $cluster chr $chr
plink2 --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz --chr chr${chr} --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink/ref.ac1.$cluster.$treat.${chr}
fi
done
done

R
library(data.table)
library(plyr)
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
combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
load(opfn)
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
treat="CTRL"
method="voom"

lapply(names(counts_ls),function(cluster){
    cat("running ",cluster,treat,"\n")
## load normalized data / or residuals:
  Res <- fread(paste0(outFolder,cluster,".",treat,".residuals_",method,".txt"))
  all_counts_bed <- merge(geneIDs[,c("chr","strand_start","strand_end","ensgene")],Res,by.x="ensgene",by.y="V1")
  all_counts_bed <- all_counts_bed %>% relocate(ensgene,.after =strand_end)  # have to use ensgene as there was multi gene symbols
  #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
  names(all_counts_bed)[c(2:4)] <- c("start","end","gene_id")
  #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
  fwrite(all_counts_bed, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"phenotypes.",cluster,".",treat,".residuals_",method,".bed"))
})

#treat="CTRL"
#method="voom"
#fwrite(cvs, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"covariates/PCcovariates_",method,".",cluster,".",treat,".txt"))
#for i in $(seq 1 20); do 
#  head -n $(($i+1)) /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/covariates/PCcovariates_${method}.$cluster.$treat.txt > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/covariates/$cluster.$treat.PC1-$i.covariates_${method}.txt
#done

module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0
method="voom"
#ALOFT
treat="CTRL"
for cluster in C0 C1 C10 C11 C2 C3 C4 C5 C6 C7 C8 C9;do
i=`ls -1 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/phenotypes.$cluster.$treat.residuals_${method}.bed | grep -v 'sort'`
  echo "running " $i
  #less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sortBed -i"}' > ${i%.*}.sort.bed
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sort -k1,1 -k2,2n "}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
done
module swap gnu7/7.3.0 gnu9



cd /wsu/home/groups/piquelab/cindy/jaxqtl
pip3 install tensorqtl

  WARNING: The scripts f2py, f2py3 and f2py3.7 are installed in '/wsu/home/fh/fh85/fh8591/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
mlxtend 0.18.0 requires scipy>=1.2.1, but you have scipy 1.1.0 which is incompatible.
multiqc 1.7 requires matplotlib<3.0.0,>=2.1.1, but you have matplotlib 3.0.3 which is incompatible.
mxnet 1.5.0 requires requests<3,>=2.20.0, but you have requests 2.19.1 which is incompatible.

cat >> ~/.bashrc
export PATH="/wsu/home/fh/fh85/fh8591/.local/bin:$PATH"
source ~/.bashrc

pip3 install scipy -U
pip3 install matplotlib -U
pip3 install multiqc -U
pip3 install requests -U

pip3 install tensorqtl
pip3 install importlib-metadata -U

conda create -n tensorqtl_p3.11 python=3.11
conda activate tensorqtl_p3.11
pip3 install tensorqtl
python3 -m tensorqtl --help
#ModuleNotFoundError: No module named 'pandas_plink'
pip3 install pandas_plink

python3 -m tensorqtl --help
Warning: 'rfunc' cannot be imported. R with the 'qvalue' library and the 'rpy2' Python package are needed to compute q-values.
pip3 install 'rpy2<=3.5.12'

conda config --add channels conda-forge
conda config --set channel_priority strict 
conda install r-essentials

conda install -c bioconda bioconductor-qvalue

python3
from rpy2.robjects.packages import importr

OSError: cannot load library '/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib/libR.so': libRblas.so: cannot open shared object file: No such file or directory

export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH

python3
from rpy2.robjects.packages import importr
#works!

python3 -m tensorqtl --help

conda env export > /wsu/home/groups/piquelab/cindy/tensorqtl.environment.yaml

conda env create --name tensorqtl_p3.11_env --file=/wsu/home/groups/piquelab/cindy/tensorqtl.environment.yaml
conda activate tensorqtl_p3.11_env
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
python3 -m tensorqtl --help | head


data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
mkdir -p ${out_path}
mkdir -p ${data_path}/logs

cluster="C4"
treat="CTRL"
chr=1
PC=2
#pheno="${data_path}/counts/phenotypes.${cluster}.${treat}.sort.bed.gz"
geno="${data_path}/residuals/vcf_plink/ref.ac1.${cluster}.${treat}.${chr}" # prefix for plink triplet files
# prefix for output file
out_prefix="${cluster}_${treat}_${chr}_PC${PC}_tensorqtl"
pheno=${data_path}/residuals/phenotypes.$cluster.$treat.PC1-$PC.residuals_voom_PCregress2step.sort.bed.gz

echo $cluster $treat $chr $PC
echo $geno
echo $pheno
echo $out_prefix

#test
python3 -m tensorqtl ${geno} ${pheno} ${out_prefix} \
    --mode cis \
    --fdr 0.1 \
    -o ${out_path}
#    --chunk_size 'chr' \

#test with non-pc regressed residualscluster="C4"
treat="CTRL"
chr=1
PC=2
for chr in {2..22};do
  echo running chr $chr
geno="${data_path}/residuals/vcf_plink/ref.ac1.${cluster}.${treat}.${chr}" # prefix for plink triplet files
# prefix for output file
out_prefix="${cluster}_${treat}_${chr}_PC${PC}_tensorqtlr"
pheno=${data_path}/residuals/phenotypes.$cluster.$treat.residuals_voom.sort.bed.gz
covar="${data_path}/residuals/covariates/${cluster}.${treat}.PC1-$PC.covariates_voom.txt"

echo $cluster $treat $chr $PC
echo $geno
echo $pheno
echo $covar
echo $out_prefix

python3 -m tensorqtl ${geno} ${pheno} ${out_prefix} \
    --mode cis \
    --covariates ${covar} \
    --fdr 0.1 \
    -o ${out_path}
done


#HPC test
data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
cluster="C0"
treat="CTRL"
PC=2
mkdir -p ${out_path}
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
#for cluster in C0 C4; do
for chr in {1..22};do
  for PC in {1..10}; do
  njobs=`squeue -u fh8591 -r| wc -l`
  maxjobs=500
  while [ "$njobs" -gt "$maxjobs" ];do 
  echo waiting for jobspace, sleeping ...
  sleep 300 
  njobs=`squeue -u fh8591 -r| wc -l`
  done #end while
      if [ -f "${out_path}/${cluster}_${treat}_${chr}_PC${PC}_tensorqtlr.cis_qtl.txt.gz" ] ; then
      echo " Output already exists. Skipping..."
    else
    echo submitting cluster $cluster chr $chr PC $PC
    sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",PC="$PC",data_path="$data_path",out_path="$out_path" ${data_path}/tensorQTL/src/run_tensorqtl.sh 
    sleep 1
  fi
  done #PC end
done #chr end
done #cluster end

#check output has all cols (19)
for i in tensorQTL/output/*.txt.gz; do  
  numcol=`less $i | head -n1 | awk '{print NF}'`; 
  echo $i has $numcol; done

#check what has run
treat="CTRL"
rm ${out_path}/not_completed.txt ${out_path}/completed.txt
rm ${out_path}/not_completed.txt ${out_path}/not_completed.txt
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt `;do
for chr in {1..22};do
  for PC in {1..10}; do
    out_prefix="${cluster}_${treat}_${chr}_PC${PC}_tensorqtlr" # prefix for output file
    out="${out_path}/${out_prefix}.cis_qtl.txt.gz"
    #echo $out
    if [ -f "$out" ] ; then
      #echo " Output already exists."
      echo "${cluster}_${treat}_${chr}_${PC}" >> ${out_path}/completed.txt
    else
      echo $out not run
      echo "${cluster}_${treat}_${chr}_${PC}" >> ${out_path}/not_completed.txt
    fi
done
    #sleep 1
done
done
#25812 - 10305
less ${out_path}/completed.txt | wc -l #440
less ${out_path}/not_completed.txt | wc -l #0

#combining 


############################################################
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt | grep -v "C4"`;do
for chr in {1..22};do
  for PC in {1..20}; do

python3 -m tensorqtl ${geno} ${pheno} ${out} \
    --mode cis \
    --chunk_size 'chr' \
    --fdr 0.1
    -o ${out_path}

    --covariates ${covar} \
--chunk_size CHUNK_SIZE
                        For cis-QTL mapping, load genotypes into CPU memory in chunks of chunk_size variants, or by chromosome if chunk_size is
                        'chr'.
  --fdr FDR             FDR for cis-QTLs
-o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory


