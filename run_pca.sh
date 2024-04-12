#!/bin/bash

module load misc

VCF=$VCF
OUT_FOLDER=$OUT_FOLDER

cd $OUT_FOLDER

##
plink2 --vcf $VCF \
        --make-bed --out ref_bi \
        --max-alleles 2

##prunning
plink2 --bfile ref_bi --double-id --allow-extra-chr \
   --set-missing-var-ids @:# \
   --indep-pairwise 50 10 0.1 --out ref

###
plink2 --bfile ref_bi --double-id --allow-extra-chr \
   --set-missing-var-ids @:# \
   --extract ref.prune.in \
   --make-bed --pca  --out ref