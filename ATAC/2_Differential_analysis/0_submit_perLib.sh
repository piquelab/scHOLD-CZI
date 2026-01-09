#!/bin/bash/

cd $PWD
opt=option_nFeature15K
cat HOLD.lib.txt | grep -v HOLD7 |  \
while read lib;
do 
###
   sbatch -q highmem --mem=100G --time=10-24:00:00 -N 1-1 -n 1 --job-name=pseudobulk_${lib} --output=slurm_pseudobulk_${opt}_${lib}.outs --wrap "
   module unload gnu7;
   module load gnu9 R;
   Rscript --vanilla 0_generate_pseudobulk_perLib.R -cond ${lib} -opt ${opt} > zzz_pseudobulk_${opt}_${lib}.Rout"
   echo ${opt} ${lib}
   sleep 1;
###
done 
