#!/bin/bash/

cd $PWD
option=option_nFeature15K_cluster_res0.12
treat=CTRL ## LPS, CTRL 
## adj_model=combat_seq ## no_adjust
## num_core=1 ### default number 1, suggest the range 1~4
split_var=psychosocial_top10.txt
## split_var=cytokines_var11.txt 
# split_var=PFAS_vars.txt 

cat ./Cluster/${option}.txt | sed -n '9,$p' | \
while read cl;
do 

## cat split_psycho_files.txt | \
## while read split_var;
## do  
###
   sbatch -q express -p erprp  --mem=10G --time=4-24:00:00 -N 1-1 -n 1 --job-name=motif_psycho_${cl} --output=slurm_motif_psycho_cluster${cl}.outs --wrap "
   module unload gnu7;
   module load gnu9 R;
   Rscript --vanilla 2.1_differentialMotifs.R -opt ${option} -cl ${cl} -psycho ${split_var} -cond ${treat} > zzz_motif_${split_var}_cluster${cl}.Rout"
   # Rscript --vanilla 1_differentialPeaks.R ${option} ${th} ${cl} ${split_var} ${treat} combat_seq > zzz_DESeq_psycho_th${th}_cluster_${cl}_${split_var}_${treat}.Rout
   echo ${option} cluster ${cl} condition ${treat} ${split_var}
   sleep 1;
###
## done

done 
