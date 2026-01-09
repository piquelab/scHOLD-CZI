#!/bin/bash/

cd $PWD

opt="option_nFeature15K"
ncpu=5
cat HOLD.lib.txt | grep -v HOLD7 | \
while read lib;
do
   sbatch -q highmem --mem=500G --time=3-24:00:00 -N 1-1 -n ${ncpu} --job-name=getPeak_${lib} --output=slurm_getPeak_${opt}_${lib}.outs --wrap "
   module unload gnu7;
   module load gnu9 R;
   Rscript --vanilla 3.0_getPeak.R ${lib} ${opt} ${ncpu} > zzz_getPeak_${opt}_${lib}.Rout"
   echo library ${lib} ${opt}
   sleep 1;
done 

   
