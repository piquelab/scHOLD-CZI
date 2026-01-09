#!/bin/bash/

cd $PWD

for i in {1..14}
do   
###
   sbatch -q express -p erprp  --mem=10G --time=4-24:00:00 -N 1-1 -n 1 --job-name=enrich_comb${i} --output=slurm_enrich_comb${i}.outs --wrap "
   module unload gnu7;
   module load gnu9 R;
   Rscript --vanilla 3.0_enrich_motif.R ${i} > zzz_enrich_comb${i}.Rout"
   echo cluster_var${i}
   sleep 1;

###
done 
