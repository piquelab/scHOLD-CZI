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
pip install lineax qtl

#want to create chunks to run (200-300 genes recommended per chunk) in A single-column (no header) file specifying gene identifiers
    split -1 300 -d data.csv chunk_

#Phenotypes are provided in BED format, with a single header line starting with # and the first four columns corresponding to: chr, start, end, phenotype_id, with the remaining columns corresponding to samples (the identifiers must match those in the genotype input). The BED file can specify the center of the cis-window (usually the TSS), with start == end-1, or alternatively, start and end positions, in which case the cis-window is [start-window, end+window]

#Covariates are provided as a tab-delimited tsv file dataframe (samples x covariates) with column headers.

#Genotypes can be provided in PLINK1 bed/bim/fam format. (We will accomodate other formats in the future verssions)
module load plink/2.0
cluster="C0"
treat="CTRL"
plink2 --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz --make-bed --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat
