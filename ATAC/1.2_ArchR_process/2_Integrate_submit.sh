#!/bin/bash/

cd $PWD

for i in {1..5};
do
   sbatch -q highmem --mem=500G --time=3-24:00:00 -N 1-1 -n 3 --job-name=integrate_RNA_res${i} --output=slurm_integrate_res${i}.outs --wrap "
   module unload gnu7;
   module load gnu9 R;
   Rscript --vanilla 2_Integrate_scRNA.R $i > zzz_res${i}.Rout"
   echo rna resolution${i}
   sleep 1
done 

   
