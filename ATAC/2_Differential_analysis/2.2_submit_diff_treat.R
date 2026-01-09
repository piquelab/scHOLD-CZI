#!/bin/bash/

cd $PWD
option=Cluster_res0.07
th=20

cat ./Cluster/${option}.cl.txt | sed -n '1,$p' | \
while read cl;
do 

###
   sbatch -q express -p erprp  --mem=20G --time=4-24:00:00 -N 1-1 -n 1 --job-name=diff_motif_${cl} --output=slurm_diff_motif_cluster${cl}.outs --wrap "
   module unload gnu7;
   module load gnu9 R;
   Rscript --vanilla 2.2_diffMotifs_treat.R  -cl ${cl}  > zzz_motif_treat_th${th}_cluster_${cl}.Rout"
   echo th ${th} cluster ${cl}
   sleep 1;
###
## done

done 
