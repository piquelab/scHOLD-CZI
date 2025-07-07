@rpique I've been running fastQTL using the manually regressed PCs. I have noticed some runs (~500) end early with the error: 
gsl: beta.c:44: ERROR: domain error
Default GSL error handler invoked.
I set to rerun with higher memory and time just in case it was a node issue. The rerun has come across the same error. Online I saw that it could be due to the low variability of some of the transcript expression levels between the samples or they are always the same genotype for a variant. 
This example I put here the histogram of the residuals going into fastQTL. The gene expression variability for the failed gene is a bit less than the gene that passed but it isn't super obvious. The genotypes that passed are: 
 115532 0|0
   6501 0|1
   6596 1|0
   6803 1|1
failed: 
 264716 0|0
   9553 0|1
   9045 1|0
   7196 1|1
This example gene passed in 2 other clusters so I think the genotype is not the issue. The variability when it passed is greater, although not by a ton. 
Do you have any suggestions for how I should handle these situations? ie. should I introduce an additional filter? remove genes that failed? 

failed: 
  77171 0|0
   2349 0|1
   2418 1|0
   1662 1|1
passed:
  77171 0|0
   2349 0|1
   2418 1|0
   1662 1|1





cat slurm/rerun/check.txt | while read i || [[ -n $i ]]; do 
ensg=$(echo "$i" | awk '{print $1}')
bed=$(echo "$i" | awk '{print $2}')
less $bed | grep $ensg | cut -f1-4 | awk "$3>$4"

#current fail
# C0 PC1-7 chunk 22 ENSG00000148019

#check to see if the failed gene ever runs 
library(data.table)
library(qqman)
library(qvalue)
library(ggplot2)
library(tidyverse)
library(plyr)

cluster="C0"
treat="CTRL"
pcnum="7"
ensg="ENSG00000148019"
passfail="fail"

FDR <- 0.1

myDir <- "/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/" #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep("permutations_10PCregress2step.eQTL.txt.gz", filenames)] #pick specific files from list
#filenames <- filenames[grepl(paste0(cluster,".",treat), filenames)] #pick specific files from list
data_names <- gsub(".permutations.eQTL.txt.gz", "", filenames) #remove file ending
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = F)[,analysis:=data_names[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
#combining data
all_PCs <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(all_PCs) <- data_names
# add colnames:
for(i in 1:length(all_PCs)){
    colnames(all_PCs[[i]]) <- c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval","analysis")
}
all_PCsdf <- as.data.frame(Reduce(function(x, y) rbind(x, y),all_PCs,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
s <- subset(all_PCsdf, pid==ensg)

#C0.CTRL.PC1-10

#check variability of gene expression residuals 
#checked C3.CTRL.PC1-13 ENSG00000146414 (passed)
#checked ENSG00000164506 failed in C3.CTRL.PC1-13 (then passed in C4.CTRL.PC1-7)
library(data.table)
library(ggplot2)
cluster="C0"
treat="CTRL"
pcnum="7"
ensg="ENSG00000148019"
passfail="fail"

cluster="C0"
treat="CTRL"
pcnum="10"
passfail="passed"

bed <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.",cluster,".",treat,".PC1-",pcnum,".residuals_qnorm_10PCregress2step.sort.bed.gz"))
g <- subset(bed, ensgene==ensg)
png(width = 8, height = 8, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/",passfail,"gene_",ensg,"_hist.png"), pointsize=12, 
  bg = "transparent", units = "in", res = 1200)
ggplot(data=melt(g[,-c(1:4)]), aes(value)) + 
ggtitle(paste0(ensg," ", passfail," in ",cluster," ", treat, ".PC1-",pcnum))+
theme_bw()+
geom_histogram(bins=50)
dev.off()


#check variability of genotypes
module load bcftools/1.19
cluster="C0"
treat="CTRL"
pcnum="7"
ensg="ENSG00000148019"
passfail="fail"

cluster="C0"
treat="CTRL"
pcnum="10"
passfail="passed"

region=`less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step.sort.bed.gz | grep $ensg |awk '{print $1":"$2"-"$3}'`
tabix /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz $region -h > $ensg.$cluster.$treat.PC1-$pcnum.$passfail.vcf
bcftools view $ensg.$cluster.$treat.PC1-$pcnum.$passfail.vcf | bcftools query -f'[%CHROM:%POS:%REF:%ALT:%SAMPLE:%GT\n]' | cut -f6 -d":" | sort | uniq -c





for i in /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/rerun/*.out; do 
    if grep -q $ensg $i; then
        if grep -q "Running time:" $i; then
        echo "file ran to completion"
    run=`grep "Scanning phenotype data" $i | sed 's/Scanning phenotype data in \[\/rs\/rs_grp_scaloft\/scALOFT_2024\/cindy_analysis\/fastQTL\/phenotypes_fastqtl.//' | sed 's/.residuals_qnorm_10PCregress2step.sort.bed.gz\]//'`
    j=`grep "Chunk processed" $i | sed 's/* Chunk processed //' | sed 's/ \/.*//' | sed 's/ //'`
    cluster=`echo $run | sed 's/[.].*//'`
    treat=`echo $run | sed 's/^[^\.]*\.//' | sed 's/[.].*//'`
    pcnum=`echo $run | sed 's/.*[-]//'`
    echo $cluster $treat PC$pcnum chunk$j within $i finished
    echo $cluster $treat $pcnum $j $i >> slurm/ENSG00000164506_finished.txt
    else 
    echo $cluster $treat PC$pcnum chunk$j within $i didnt finish
    echo $cluster $treat PC$pcnum chunk$j within $i >> slurm/ENSG00000164506_failed.txt
    fi 
    fi
done 
#it does run so what if I try running interactively the one that failed
grep -A5 "ENSG00000164506" /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/slurm/rerun/29756340.out
/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.C4.CTRL.PC1-7.residuals_qnorm_10PCregress2step.sort.bed.gz

cluster="C3"
treat="CTRL"
pcnum=13
j=30
fastQTL --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz \
--bed /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.$cluster.$treat.PC1-$pcnum.residuals_qnorm_10PCregress2step.sort.bed.gz \
--permute 1000 10000 --window 1e6 \
--out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/$cluster.$treat.PC1-$pcnum.permutations_10PCregress2step.chunk$j.txt.gz \
--chunk $j 30

#check maf from vcf
module swap gnu9 gnu7/7.3.0
module load vcftools/0.1.16
vcftools --gzvcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz --freq --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.freq.txt 
head -n 1000 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.freq.txt.frq > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.freq.test

library(data.table)
cluster="C9"
treat="CTRL"
gtinfo <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.",cluster,".",treat,".freq.test"),skip=1)
gtinfo <- transform(gtinfo, A1=sapply(strsplit(V5,":"),function(y)y[1]),A1freq=as.numeric(sapply(strsplit(V5,":"),function(y)y[2])),A2=sapply(strsplit(V6,":"),function(y)y[1]),A2freq=as.numeric(sapply(strsplit(V6,":"),function(y)y[2])),chrnum=gsub("chr","",V1))
gtaf0 <- subset(gtinfo, A1freq==0 | A2freq==0)
gtafnot0 <- subset(gtinfo, A1freq>0 & A2freq>0)
snpremid <- transform(gtafnot0, ID=paste(chrnum,V2,A1,A2,sep=":"))
fwrite(snpremid[,"ID"],quote=F, row.names=F,col.names=F,file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/",cluster,".",treat,".snpremlist.txt"))
gtinfo <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.",cluster,".",treat,".freq.txt.frq"),sep="\t",sep2="\:",header=T)

less /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz | head -n 1000 > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.test.vcf
#bcftools view -i 'ID !~ @/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/$cluster.$treat.snpremlist.txt' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.test.vcf -O z /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.test.vcf.gz
bcftools view -q 0.001:minor /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.test.vcf > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.testno0maf.vcf

#bcftools view -i 'ID !~ @/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/$cluster.$treat.snpremlist.txt' -O z /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.test.vcf -o /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz
#bcftools query -f'[%CHROM:%POS:%REF:%ALT:%SAMPLE:%GT:%RAF:%AF\n]' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz | head -n 20000 > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.$cluster.$treat.gtinfo20k.txt
#gtinfo <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.",cluster,".",treat,".gtinfo.txt"),sep=":",header=F)





#trying C0 with the filter to see if more runs - go from X snps to X snps
bcftools view -q 0.001:minor /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.C0.$treat.filtered.vcf.gz > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.C0.$treat.no0maf.vcf
bgzip /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.C0.$treat.no0maf.vcf && tabix -p vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.C0.$treat.no0maf.vcf.gz
bcftools view /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.C0.$treat.filtered.vcf.gz | wc -l
bcftools view /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/vcf/ref.ac1.C0.$treat.no0maf.vcf.gz | wc -l
26486194-25083138= 1,403,056

#run alt (nomaf0 vcf) vs alt2 (original vcf) and see how much completion
runnum=1
treat="CTRL"
method="voom"
fastqtlfolder=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL
cluster=C0
njobs=`squeue -u fh8591 -r| wc -l`
maxjobs=500
#if [ "$njobs" -gt "$maxjobs" ]; then
while [ "$njobs" -gt "$maxjobs" ];do 
echo waiting for jobspace, sleeping ...
sleep 300 #most jobs run in 800 sec, 1000 sleep was too little -- start should be 2000, then lower for reruns
njobs=`squeue -u fh8591 -r| wc -l`
done
echo submitting $cluster
sbatch --export=cluster="${cluster}",treat="${treat}",method="${method}",runnum="${runnum}" ${fastqtlfolder}/src/run_fastqtl_array50_altout.sh

njobs=`squeue -u fh8591 -r| wc -l`
while [ "$njobs" -gt "$maxjobs" ];do 
echo waiting for jobspace, sleeping ...
sleep 300 #most jobs run in 800 sec, 1000 sleep was too little -- start should be 2000, then lower for reruns
njobs=`squeue -u fh8591 -r| wc -l`
done
echo submitting $cluster
sbatch --export=cluster="${cluster}",treat="${treat}",method="${method}",runnum="${runnum}" ${fastqtlfolder}/src/run_fastqtl_array50_altout2.sh


for pcnum in $(seq 1 20); do
for j in $(seq 1 50); do
    file=${fastqtlfolder}/permutations/alt/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk${j}.txt.gz
    if gzip -t "$file" >/dev/null 2>&1; then #this says the file is complete 
            echo $cluster $treat $pcnum $j >> ${fastqtlfolder}/altfilecomplete.txt
    else
            echo $cluster $treat $pcnum $j >> ${fastqtlfolder}/altfilenotcomplete.txt
    fi
done
done
less ${fastqtlfolder}/altfilecomplete.txt | wc -l #358
less ${fastqtlfolder}/altfilenotcomplete.txt | wc -l #642

for pcnum in $(seq 1 20); do
for j in $(seq 1 50); do
    file=${fastqtlfolder}/permutations/alt2/$cluster.$treat.PC1-$pcnum.permutations_PCregress2step.chunk${j}.txt.gz
    if gzip -t "$file" >/dev/null 2>&1; then #this says the file is complete 
            echo $cluster $treat $pcnum $j >> ${fastqtlfolder}/alt2filecomplete.txt
    else
            echo $cluster $treat $pcnum $j >> ${fastqtlfolder}/alt2filenotcomplete.txt
    fi
done
done
less ${fastqtlfolder}/alt2filecomplete.txt | wc -l #335
less ${fastqtlfolder}/alt2filenotcomplete.txt | wc -l #665
