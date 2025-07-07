#Now I have this script which actually loops through the chunks + PCs using a job array, 
#increases the sleep time and attempts a few times per chunk in case the job fails due to a memory access issue or some other error:
#fastqtl_chunks.txt is the file that gives the number of PCs/chunks and then the two lines of code below make the 
#log directory and submit the job
mkdir -p /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs
mkdir -p /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/src/

for i in $(seq 1 2); do
    for j in $(seq 1 2); do
        echo "$i $j"
  done
done > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/fastqtl_test_chunks.txt
sed -n '2,5p' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/endoffileerror${runnum}.txt > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/endoffileerror${runnum}_test.txt
runnum=2
rm /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/test/*
sbatch --export=cluster="${cluster}",treat="${treat}",method="${method}",runnum="${runnum}" /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/src/run_fastqtl_test_array.sh

for i in $(seq 1 20); do
    for j in $(seq 1 30); do
        echo "$i $j"
  done
done > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/fastqtl_chunks.txt

for i in $(seq 1 20); do
    for j in $(seq 1 50); do
        echo "$i $j"
  done
done > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/fastqtl_chunks50.txt

#treat="CTRL"
treat="RNA-CTRL"
method="voom"
runnum=1
#fastqtlfolder=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL
fastqtlfolder=/rs/rs_grp_schold/CZI/RNA/analysis/fastQTL
#CLUSTERS=`awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`
CLUSTERS=$(echo C0 C1 C2 C3 C4 C5) 
for cluster in $CLUSTERS;do
#for cluster in C0 C1 C10 C2 C3 C4 C5 C7 C8 C9;do
njobs=`squeue -u fh8591 -r| wc -l`
maxjobs=500
#if [ "$njobs" -gt "$maxjobs" ]; then
while [ "$njobs" -gt "$maxjobs" ];do 
echo waiting for jobspace, sleeping ...
sleep 500 #most jobs run in 800 sec, 1000 sleep was too little -- start should be 2000, then lower for reruns
njobs=`squeue -u fh8591 -r| wc -l`
done
echo submitting $cluster
sbatch --export=FASTQTLFOLDER="$fastqtlfolder",cluster="${cluster}",treat="${treat}",method="${method}",runnum="${runnum}" ${fastqtlfolder}/src/run_fastqtl_array50.sh 
#fi
done
sleep 2000
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
#for cluster in C0 C1 C10 C2 C3 C4 C5 C7 C8 C9;do
njobs=`squeue -u fh8591 -r| wc -l`
maxjobs=500
#if [ "$njobs" -gt "$maxjobs" ]; then
while [ "$njobs" -gt "$maxjobs" ];do 
echo waiting for jobspace, sleeping ...
sleep 500 #most jobs run in 800 sec, 1000 sleep was too little -- start should be 2000, then lower for reruns
njobs=`squeue -u fh8591 -r| wc -l`
done
echo submitting $cluster
sbatch --export=cluster="${cluster}",treat="${treat}",method="${method}",runnum="${runnum}" /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/src/run_fastqtl_array50.sh 
#fi
done


count=0
for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/*.out | awk '$6 == "Jun" && $7 >= 25 {print $9}'`; do 
if grep -q "Running time:" $i; then
    count=$((count+1))
fi
done

for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/*.out | awk '$6 == "Jun" && $7 >= 19 {print $9}'`; do 
j=`grep "Chunk processed" $i | sed 's/^[^0-9]*//' | sed 's/ .*//'1 | uniq`
cluster=`grep "Scanning phenotype data" $i | cut -d . -f2 | uniq`
treat=`grep "Scanning phenotype data" $i | cut -d . -f3 | uniq`
cluster=`grep "Scanning phenotype data" $i | cut -d . -f2 | uniq`
pcnum=`grep "Scanning phenotype data" $i | cut -d . -f4 | cut -d - -f2 | uniq`
#   jobid=`echo $i | sed 's/.*\///' | sed 's/.out//'`
    #echo $cluster $treat $pcnum $j
    if grep -q "Output already exists" $i; then
        rm $i ${i%.*}.err
        continue
    fi
    if grep -q "Running time:\|Output already exists" $i; then
#        #echo "file ran to completion"
        echo $cluster $treat $pcnum $j >> /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/filesran.txt
        mv $i /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/completed/
        mv ${i%.*}.err /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/completed/
    else
    file=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk${j}.txt.gz
    if gzip -t "$file" >/dev/null 2>&1; then #this says the file is complete 
            echo $cluster $treat $pcnum $j >> /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/filesran.txt
            mv $i /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/completed/
            mv ${i%.*}.err /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/completed/
        else
    #echo $cluster $treat PC$pcnum chunk$j within $i didnt finish
    echo $cluster $treat $pcnum $j >> /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/needtorun_conditions.txt
    mv $i /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/
    mv ${i%.*}.err /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/
    fi
    fi 
done 
for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/fastqtl_*.out | awk '$6 == "Jun" && $7 >= 19 {print $9}'`; do 
if grep -q "arrayoutsidePCxchunkcombos\|array outside PC x chunk combos" $i; then
    rm $i ${i%.*}.err
    continue
elif grep -q "arrayoutsidePCxchunkcombos\|array outside PC x chunk combos" ${i%.*}.err; then
    rm $i ${i%.*}.err
    continue
fi
done
less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/filesran.txt | sort | uniq > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/filesran_unique.txt
less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/needtorun_conditions.txt | sort | uniq > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/needtorun_unique.txt

R
library(data.table)
runnum=2
f <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/needtorun_unique.txt",fill=4)
c <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/filesran_unique.txt",fill=4)
f <- transform(f, treat="CTRL", pcnum=as.integer(ifelse(V2=="CTRL",V3, V2)), chunk=ifelse(V2=="CTRL",V4, V3))
c <- transform(c, treat="CTRL", pcnum=as.integer(ifelse(V2=="CTRL",V3, V2)), chunk=ifelse(V2=="CTRL",V4, V3))
f <- f[,c(1,5,6,7)]
c <- c[,c(1,5,6,7)]
notXY1 <- unique(merge(f,c,all.x = TRUE)[!merge(f,c)])
#f <- transform(f, combo=paste(V1,V2,V3,V4,sep="_"))
#c <- transform(c, combo=paste(V1,V2,V3,V4,sep="_"))
#fnotc <- subset(f, !combo %in% c$combo)
fwrite(notXY1, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/endoffileerror",runnum,".txt"), quote=F,sep='\t',row.names=F,col.names=F)
dim(notXY1)
table(notXY1$V1)
q()
n
#cat /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/endoffileerror${runnum}.txt | wc -l


i=fastQTL/logs/failedruns/fastqtl_30742092_968.out
mkdir /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/bednotfound
for i in /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/*.out;do
if grep -q "bed file not readable" $i; then
    #echo bed didnt read in
    mv $i /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/bednotfound/
    mv ${i%.*}.err /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/bednotfound/
fi
done

for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/*.out | awk '$6 == "Jul" && $7 > 1 {print $9}'`; do 
    if grep -q "Output already exists" $i; then
        rm $i ${i%.*}.err
        continue
    fi
if grep -q "arrayoutsidePCxchunkcombos\|array outside PC x chunk combos" $i; then
    rm $i ${i%.*}.err
    continue
elif grep -q "arrayoutsidePCxchunkcombos\|array outside PC x chunk combos" ${i%.*}.err; then
    rm $i ${i%.*}.err
    continue
fi
done

#remove old logs
for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/fastqtl_*.out | awk '$6 == "Jun" && $7 <= 15 {print $9}'`; do
    rm $i ${i%.*}.err
done

#this reads through failed runs output, compiles the errors and then removes those due to gsl (not useful to QC)
rm /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/compiled_error.txt
for i in `ls -al /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/fastqtl_*.out | awk '$6 == "Jun" && $7 >= 16 {print $9}'`; do 
node=`head -n1 $i | cut -d "]" -f2 | cut -d ":" -f2 | sed 's/ //g'`
error=`head -n1 ${i%.*}.err`
jobname=`echo ${i##*/} | cut -d "." -f1`
echo $node $jobname error=${error}>> /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/compiled_error.txt
done
less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/compiled_error.txt | sort | uniq > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/compiled_erroru.txt
#echo $runnum $treat $method >> /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/tracknodestest.txt
#non gsl errors
grep -q "gsl" /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/compiled_erroru.txt > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/failedruns/compiled_erroru_nongsl.txt

#torun just missing jobs 
runnum=2
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
    grep "$cluster" /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/endoffileerror${runnum}.txt | sort | uniq > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/endoffileerror${runnum}_${cluster}.txt
done

#rm /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/*
runnum=3
treat="CTRL"
method="voom"
fastqtlfolder=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL
#for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
for cluster in `less ${fastqtlfolder}/filenotcomplete.txt | cut -d" " -f1 | sort | uniq`; do
njobs=`squeue -u fh8591 -r| wc -l`
maxjobs=500
#if [ "$njobs" -gt "$maxjobs" ]; then
while [ "$njobs" -gt "$maxjobs" ];do 
echo waiting for jobspace, sleeping ...
sleep 300 #most jobs run in 800 sec, 1000 sleep was too little -- start should be 2000, then lower for reruns
njobs=`squeue -u fh8591 -r| wc -l`
done
echo submitting $cluster
#sbatch --export=cluster="${cluster}",treat="${treat}",method="${method}",runnum="${runnum}" ${fastqtlfolder}/src/run_fastqtl_array50.sh 
sbatch --export=cluster="${cluster}",treat="${treat}",method="${method}",runnum="${runnum}" ${fastqtlfolder}/src/run_fastqtl_array50_extramem.sh 
done

sleep 3500
less ${fastqtlfolder}/filecomplete.txt | wc -l
less ${fastqtlfolder}/filenotcomplete.txt | wc -l
prevclusters=`less ${fastqtlfolder}/filenotcomplete.txt | cut -d" " -f1 | sort | uniq`
rm ${fastqtlfolder}/filecomplete.txt
rm ${fastqtlfolder}/filenotcomplete.txt
treat="CTRL"
#for cluster in `awk 'NR>1 {print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
for cluster in $prevclusters; do
#for cluster in $CLUSTERS;do
echo running $cluster 
for pcnum in $(seq 1 20); do
for j in $(seq 1 50); do
    file=${fastqtlfolder}/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk${j}.txt.gz
    if gzip -t "$file" >/dev/null 2>&1; then #this says the file is complete 
            echo $cluster $treat $pcnum $j >> ${fastqtlfolder}/filecomplete.txt
    else
            echo $cluster $treat $pcnum $j >> ${fastqtlfolder}/filenotcomplete.txt
    fi
done
done
done
less ${fastqtlfolder}/filecomplete.txt | wc -l
less ${fastqtlfolder}/filenotcomplete.txt | wc -l
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
    grep "$cluster" ${fastqtlfolder}/filenotcomplete.txt > ${fastqtlfolder}/logs/endoffileerror${runnum}_${cluster}.txt
done

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
less ${fastqtlfolder}/logs/endoffileerror${runnum}_${cluster}.txt | wc -l
done

for i in {0..9};do
    for j in {0..9}; do
    rm fastQTL/logs/*33${j}${i}*.out
    rm fastQTL/logs/*33${j}${i}*.err
done
done

#run single chunks left
cluster="C10"
treat="CTRL"
pcnum=7
j=29
echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
    
ensg=ENSG00000175550
BED_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz
less $BED_FILE | grep -v $ensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg2=ENSG00000198851
allensg=$(echo ${ensg}\|$ensg2)
ensgname=$(echo ${ensg}\_$ensg2)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50

cluster="C8"
less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/endoffileerror${runnum}_${cluster}.txt
treat="CTRL"
pcnum=3
j=35
echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg=ENSG00000113312
BED_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz
less $BED_FILE | grep -v $ensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg2=ENSG00000134516
allensg=$(echo ${ensg}\|$ensg2)
ensgname=$(echo ${ensg}\_$ensg2)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg3=ENSG00000072803
allensg=$(echo ${ensg}\|$ensg2\|$ensg3)
ensgname=$(echo ${ensg}\_$ensg2\_$ensg3)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg4=ENSG00000198624
allensg=$(echo ${ensg}\|$ensg2\|$ensg3\|$ensg4)
ensgname=$(echo ${ensg}\_$ensg2\_$ensg3\_$ensg4)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50

cluster="C4"
less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/endoffileerror${runnum}_${cluster}.txt
treat="CTRL"
pcnum=1
j=23
echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg=ENSG00000065135
BED_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz
less $BED_FILE | grep -v $ensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg2=ENSG00000143401
allensg=$(echo ${ensg}\|$ensg2)
ensgname=$(echo ${ensg}\_$ensg2)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg3=ENSG00000158710
allensg=$(echo ${ensg}\|$ensg2\|$ensg3)
ensgname=$(echo ${ensg}\_$ensg2\_$ensg3)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg4=ENSG00000143198
allensg=$(echo ${ensg}\|$ensg2\|$ensg3\|$ensg4)
ensgname=$(echo ${ensg}\_$ensg2\_$ensg3\_$ensg4)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50

pcnum=2
j=35
echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg=ENSG00000164405
BED_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz
less $BED_FILE | grep -v $ensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg2=ENSG00000131507
allensg=$(echo ${ensg}\|$ensg2)
ensgname=$(echo ${ensg}\_$ensg2)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg3=ENSG00000175416
allensg=$(echo ${ensg}\|$ensg2\|$ensg3)
ensgname=$(echo ${ensg}\_$ensg2\_$ensg3)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50

pcnum=3
j=12
echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg=ENSG00000134313
BED_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz
less $BED_FILE | grep -v $ensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg2=ENSG00000119185
allensg=$(echo ${ensg}\|$ensg2)
ensgname=$(echo ${ensg}\_$ensg2)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg3=ENSG00000218739
allensg=$(echo ${ensg}\|$ensg2\|$ensg3)
ensgname=$(echo ${ensg}\_$ensg2\_$ensg3)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50

pcnum=10
j=35
echo cluster $cluster treat $treat PCnum $pcnum chunk $j running
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz \
    --permute 1000 10000 --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg=ENSG00000164405
BED_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_${method}_PCregress2step.sort.bed.gz
less $BED_FILE | grep -v $ensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensg}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg2=ENSG00000113658
allensg=$(echo ${ensg}\|$ensg2)
ensgname=$(echo ${ensg}\_$ensg2)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50
ensg3=ENSG00000169045
allensg=$(echo ${ensg}\|$ensg2\|$ensg3)
ensgname=$(echo ${ensg}\_$ensg2\_$ensg3)
less $BED_FILE | egrep -v $allensg | bgzip > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz
tabix -p bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
BED_GENERM_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step_rem${ensgname}.sort.bed.gz 
/wsu/home/groups/piquelab/apps/el7/misc/bin/fastQTL \
    --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
    --bed "$BED_GENERM_FILE" \
    --permute 1000 10000 \
    --window 1e6 \
    --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz \
    --chunk $j 50








#did everything finish?
runnumplus=4
rm fastQTL/endoffileerror${runnumplus}.txt
cat fastQTL/endoffileerror${runnum}.txt | while read i || [[ -n $i ]];do
    cluster=$(echo "$i" | awk '{print $1}')
    treat=$(echo "$i" | awk '{print $2}')
    pcnum=$(echo "$i" | awk '{print $3}' | sed 's/PC//')
    j=$(echo "$i" | awk '{print $4}')
if ! zcat  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk$j.txt.gz
then
    echo $cluster $treat $pcnum $j >> fastQTL/endoffileerror${runnumplus}.txt
fi
done 

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
    for i in $(seq 1 20); do
    echo combining cluster $cluster pcnum $i
        for j in $(seq 1 50); do
        zcat  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations_PCregress2step.chunk*.txt.gz | gzip -c >  /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$i.permutations_PCregress2step.eQTL.txt.gz
        done
    done
done

#############################################################################################3
###########################################################################################

#!/bin/bash
#SBATCH --job-name=FastQTL_array         # Job name
#SBATCH --output=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/fastqtl_%A_%a.out  # Standard output log per array task
#SBATCH --error=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/logs/fastqtl_%A_%a.err   # Standard error log per array task
#SBATCH --array=0-599                    # Adjust based on number of (PC, chunk) combinations
#SBATCH --ntasks=2                       # Number of tasks (FastQTL uses multithreading optionally)
#SBATCH --mem=12G                        # Memory requirement
#SBATCH --time=2-18:00:00                # Max runtime (2 days, 18 hours)
#SBATCH -q primary                       # Queue/partition

### -------- USER CONFIGURABLE VARIABLES -------- ###
# Input configuration file with (PC CHUNK) pairs — 600 lines for 20 PCs × 30 chunks
CHUNK_LIST="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/fastqtl_chunks.txt"

# File path templates
#COV_DIR="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GEPCA/covariates"
VCF_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz
N_CHUNKS=30        # Total number of chunks

# Retry logic
MAX_ATTEMPTS=5
SLEEP_TIME=10

### -------- LOAD MODULES -------- ###
module load misc2

### -------- PARSE INPUTS -------- ###
# Get the i-th line for this task from the list (SLURM_ARRAY_TASK_ID is 0-based)
IJ=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CHUNK_LIST")
PC=$(echo $IJ | awk '{print $1}')
CHUNK=$(echo $IJ | awk '{print $2}')

# Construct file paths dynamically
BED_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/residuals/phenotypes_fastqtl.$cluster.$treat.PC1-$PC.residuals_${method}_PCregress2step.sort.bed.gz
OUT_FILE=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$PC.permutations_PCregress2step.chunk${CHUNK}.txt.gz

### -------- LOGGING -------- ###
echo "[$(date)] Starting FastQTL on node: $(hostname)"
echo "[$(date)] PC=$PC, CHUNK=$CHUNK"
echo "[$(date)] bed file: $BED_FILE"
echo "[$(date)] Output file: $OUT_FILE"

### -------- SKIP IF OUTPUT EXISTS -------- ###
if [ -f "$OUT_FILE" ]; then
    echo "[$(date)] Output already exists. Skipping..."
    exit 0
fi

### -------- MAIN EXECUTION LOOP WITH RETRIES -------- ###
ATTEMPT=1
while [ $ATTEMPT -le $MAX_ATTEMPTS ]; do
    if [ ! -r "$BED_FILE" ]; then
        echo "[$(date)] Attempt $ATTEMPT: bed file not readable. Retrying in $SLEEP_TIME seconds..."
        sleep $SLEEP_TIME
    else
        echo "[$(date)] bed file found. Running FastQTL..."
        fastQTL \
            --vcf "$VCF_FILE" \
            --bed "$BED_FILE" \
            --permute 1000 10000 \
            --window 1e6 \
            --out "$OUT_FILE" \
            --chunk "$CHUNK" "$N_CHUNKS" && break

        echo "[$(date)] FastQTL exited with error. Retrying in $SLEEP_TIME seconds..."
        sleep $SLEEP_TIME
    fi
    ATTEMPT=$((ATTEMPT + 1))
done

### -------- FINAL CHECK -------- ###
if [ $ATTEMPT -gt $MAX_ATTEMPTS ]; then
    echo "[$(date)] Job failed after $MAX_ATTEMPTS attempts for PC=$PC, chunk=$CHUNK"
    exit 1
else
    echo "[$(date)] Job completed successfully for PC=$PC, chunk=$CHUNK"
    exit 0
fi