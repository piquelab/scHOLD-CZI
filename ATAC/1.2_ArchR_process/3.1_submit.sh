#!/bin/bash/

cd $PWD

opt="option_nFeature15K"
ncpu=5
cat cluster.txt | \
while read cl;
do
   sbatch -q highmem --mem=500G --time=3-24:00:00 -N 1-1 -n ${ncpu} --job-name=getPeak_${cl} --output=slurm_getPeak_${opt}_${cl}.outs --wrap "
   module unload gnu7;
   module load gnu9 R;
   Rscript --vanilla 3.1_getPeak.R ${cl} ${opt} ${ncpu} > zzz_getPeak_${opt}_${cl}.Rout"
   echo cluster ${cl} ${opt}
   sleep 1;
done 

   
