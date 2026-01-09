#!/bin/bash/

cd $PWD
option=option_nFeature15K_cluster_res0.12
treat=LPS ### CTRL, LPS
num_core=4 ### default number 1, suggest the range 1~4
# split_var=psychosocial_variables

cat ./Cluster/${option}.txt | sed -n '1,8p' | \
while read cl;
do 

split_var=psychosocial_top10.txt
## split_var=cytokines_var11.txt
##split_var=PFAS_vars.txt

# cat split_psycho_files.txt | \
# while read split_var;
# do  
###
sbatch -q highmem  --mem=200G --time=4-24:00:00 -N 1-1 -n ${num_core} --job-name=DESeq_psycho_${cl}_${treat}_${split_var} --output=slurm_DESeq_psycho_${option}_cluster${cl}_${treat}_${split_var}.outs --wrap "
   module unload gnu7;
   module load gnu9 R;
   Rscript --vanilla 1_differentialPeaks.R -cl ${cl} -psycho ${split_var} -opt ${option} -cond ${treat} -core ${num_core} > zzz_DESeq_psycho_${option}_cluster_${cl}_${treat}_${split_var}.Rout"
   # Rscript --vanilla 1_differentialPeaks.R ${option} ${th} ${cl} ${split_var} ${treat} combat_seq > zzz_DESeq_psycho_th${th}_cluster_${cl}_${split_var}_${treat}.Rout
   echo  ${option} cluster ${cl}  ${split_var} ${treat}
   sleep 1;
###
# done
done 
