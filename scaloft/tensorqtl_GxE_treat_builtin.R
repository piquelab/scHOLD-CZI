#save residuals for built in tensorqtl (treatment)

library(tidyr)
library(tidyverse) #gave me issues 031725 -- hax node
library(glmnet)
library(data.table)
library(methods)
library(plyr)
library(ggplot2)
library(ggpubr)
library(irlba)

library("AnnotationHub")
ah <- AnnotationHub()
if(length(ah["AH98047"]) == 0) {
  edb <- ah[["AH75011"]]
} else {
  edb <- ah[["AH98047"]]
}
geneIDs <- genes(edb) %>%
  as.data.frame() %>% 
  setDT(keep.rownames = "ensembl_gene_id") %>%
  .[, c("ensembl_gene_id","entrezid","symbol","seqnames","start","end","strand","gene_biotype", "description")]
names(geneIDs)[c(1,2,4,8)] <- c("ensgene","entrez","chr","biotype")
geneIDs <- subset(geneIDs, biotype=="protein_coding")
geneIDs <- transform(geneIDs, chr=as.character(chr),strand=as.character(strand))
#this 1bp positioning for bed file is based on examples from qtltools and tensorqtl
geneIDs <- transform(geneIDs, strand_start=ifelse(strand=="+",start,end),strand_end=ifelse(strand=="+",start+1,end+1))
#geneIDs <- transform(geneIDs, strand_start=start,strand_end=start+1)
geneIDs <- subset(geneIDs, chr %in% c(1:22))

base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
resmethod="voom"
resFolder <- paste0(base,"residuals/")
if (!file.exists(paste0(resFolder,"counts/"))) dir.create(paste0(resFolder,"counts/"), showWarnings=F)
contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA-DEX")) #no lps vs lps-dex as too few ind
x=2
con=contrastdf[x,]
contrast <- paste0(con$treatment,"_vs_",con$control)
opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)
data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
filenames_cov <- list.files(paste0(data_path,"residuals/covariates/")) #file list from directory
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
best_df <- fread(file=paste0(outFolder,"bestPCs_table.txt"))

allsamples <- ldply(lapply(names(counts_ls),function(clus){
  data_control <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/",clus,".",con$control,".residuals_",resmethod,".txt"))
  data_treat <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/",clus,".",con$treatment,".residuals_",resmethod,".txt"))

  colnames(data_control) <- gsub("[.]","-",colnames(data_control))
  colnames(data_treat) <- gsub("[.]","-",colnames(data_treat))

  names(data_control)[c(2:length(colnames(data_control)))] <- paste0(colnames(data_control)[-c(1)],"-",con$control)
  names(data_treat)[c(2:length(colnames(data_treat)))] <- paste0(colnames(data_treat)[-c(1)],"-",con$treatment)
  inds_ctrl <- colnames(data_control)[-c(1)]
  inds_treat <- colnames(data_treat)[-c(1)]

  datam <- merge(data_control,data_treat, by=c("V1"))
  cv_data <- data.frame(Sample_ID=c(inds_ctrl,inds_treat),treats=c(rep(con$control,length(inds_ctrl)),rep(con$treatment,length(inds_treat))))

  #chr | start | end | gene id | ind1 | ind2 ...
  all_counts_bed <- merge(geneIDs[,c("chr","strand_start","strand_end","symbol","ensgene")],datam,by.x="ensgene",by.y="V1")
  all_counts_bed <- all_counts_bed %>% relocate(ensgene,.after =strand_end) %>% dplyr::select(-symbol) # have to use ensgene as there was multi gene symbols
  #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
  names(all_counts_bed)[c(2:4)] <- c("start","end","gene_id")
  #fwrite(all_counts_bed, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"counts/phenotypes.",clus,".",contrast,".bed"))

  best.PCs <- subset(best_df, cluster==clus)$PCs
  covar <- fread(paste0(data_path,"residuals/covariates/",filenames_cov[grep(paste0(clus,".",con$treatment,".PC1-",best.PCs,".covariates_voom.txt"),filenames_cov)]))
  names(covar)[c(2:length(colnames(covar)))] <- paste0(colnames(covar)[-c(1)],"-",con$control)
  pc_signif_pairs <- subset(fread(paste0(outFolder,"results/",clus,".",treat,".best_", best.PCs, ".GEPCs.txt")),qval<0.1)
  phenosub <- subset(all_counts_bed, gene_id %in% pc_signif_pairs$phenotype_id)
  
  vcfind <- colnames(fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/","vcf/ref.ac1.",clus,".",con$treatment,".vcfheader.txt")))
  # run PCA on residuals
  PCs <- prcomp_irlba(t(data_treat[,-1]), n=best.PCs)
  summary(PCs)
  mypcs <- as.data.frame(PCs$x)
  rownames(mypcs) <-inds_treat

  covs <- as_tibble(t(mypcs))
  #order to vcf
  #covsord <- covs[,match(vcfind, colnames(covs)) ]
  cvs <- cbind(id=colnames(mypcs),covs)
#    fwrite(cvs, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"covariates/PCcovariates_",resmethod,".",clus,".",con$treatment,".txt"))

  bothcovar <- cbind(covar,cvs[,-1])

  phenoind <- colnames(phenosub[,-c(1:4)])
  covarind <- colnames(bothcovar)[-1]
  common_ind <- intersect(phenoind, covarind)

  bothcovarsubcov <- subset(bothcovar, select = c(colnames(covar)[1],common_ind))
  phenoindsubcov <- subset(phenosub, select = c(colnames(phenosub)[c(1:4)],common_ind))
  cv_datasubcov <- subset(cv_data, Sample_ID %in% common_ind)
  cv_datasubcov <- transform(cv_datasubcov, treats=ifelse(treats==con$control,0,1))

  covvarv <- cv_datasubcov[,-1, drop=F]
  rownames(covvarv) <- cv_datasubcov$Sample_ID #NEED ROWNAMES
  names(covvarv)[1] <- contrast
  fwrite(covvarv,file=paste0(data_path,"residuals/covariates/",clus,".",contrast,".for_tensorqtl_int.txt"),sep="\t",col.names=T,row.names=T, quote=F)
  fwrite(bothcovarsubcov,paste0(data_path,"residuals/covariates/",clus,".",contrast,".PC1-",best.PCs,".covariates_voom.txt"),sep="\t",col.names=T,row.names=F, quote=F)
  fwrite(phenoindsubcov,paste0(data_path,"residuals/","phenotypes.",clus,".",contrast,".residuals_voom.eQTLonly.bed"),sep="\t",col.names=T,row.names=F, quote=F)
  samples <- data.frame(cluster=clus,samples=common_ind)
  return(samples)
}), data.frame)

allsamplesu <- unique(allsamples[,-1,drop=F])
fwrite(allsamplesu, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(resFolder,"vcf/sample_list.",contrast,".txt"))

less /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.gsubRI.txt | grep -v "^##" > /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.sampleheaderonly.gsubRI.txt

vcfheader <- fread("/rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.sampleheaderonly.gsubRI.txt")
dat <- select(vcfheader,starts_with("AL-"))
ctrlvcfind <- data.frame(samples=paste0(colnames(dat),"-",con$control))
treatvcfind <- data.frame(samples=paste0(colnames(dat),"-",con$treatment))
fwrite(ctrlvcfind, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(resFolder,"vcf/sample_list.",con$control,".txt"))
fwrite(treatvcfind, sep='\t', quote=F, row.names=F, col.names=F, file=paste0(resFolder,"vcf/sample_list.",con$treatment,".txt"))

module load bcftools/1.19
control="CTRL"
treatment="PHA"
bcftools reheader --samples /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/sample_list.$control.txt -o /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$control.vcf.gz /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf
bcftools reheader --samples /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/sample_list.$treatment.txt -o /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.vcf.gz /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf
tabix -p vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$control.vcf.gz
tabix -p vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.vcf.gz
bcftools merge /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$control.vcf.gz /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.vcf.gz -o /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.$control.vcf.gz
tabix -p vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.$control.vcf.gz

plink2 --vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.$treatment.$control.vcf.gz --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink/ref.ac1.$treatment.$control

#add # to header and sort and gzip
module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0
method="voom"
#ALOFT
var="PHA_vs_CTRL"
for cluster in C0 C1 C10 C11 C2 C3 C4 C5 C6 C7 C8 C9;do
  i=`ls -1 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/phenotypes.$cluster.$var.residuals_${method}.eQTLonly.bed | grep -v 'sort'`
  echo "running " $i
  #less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sortBed -i"}' > ${i%.*}.sort.bed
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sort -k1,1 -k2,2n "}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
done
module swap gnu7/7.3.0 gnu9


conda init bash
source ~/.bashrc
module swap gnu7/7.3.0 gnu9

conda activate tensorqtl_p3.11_env
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
module load R/4.3.2 #old: module load r/4.2.0

data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output"
mkdir -p ${outFolder}/interaction
out_path="${outFolder}/interaction"
var="PHA_vs_CTRL"
cluster="C1"
control="CTRL"
treatment="PHA"

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
PC=`grep -w $cluster ${outFolder}/bestPCs_table.txt | cut -f2`
geno="${data_path}/residuals/vcf_plink/ref.ac1.$treatment.$control" # prefix for plink triplet files
out_prefix="${cluster}_${var}_tensorqtlint" # prefix for output file
pheno="${data_path}/residuals/phenotypes.$cluster.${var}.residuals_voom.eQTLonly.sort.bed.gz"
covar="${data_path}/residuals/covariates/${cluster}.${var}.PC1-$PC.covariates_voom.txt"
int="${data_path}/residuals/covariates/${cluster}.${var}.for_tensorqtl_int.txt"

if [ ! -s "${out_path}/${out_prefix}.cis_qtl_top_assoc.txt.gz" ] ; then

### -------- LOGGING -------- ###
echo "[$(date)] Starting tensorQTL interaction on node: $(hostname)"
echo "Contrast=$var, CLUSTER=$cluster, PC=$PC"
echo "pheno file: $pheno"
echo "geno file: $geno"
echo "covar file: $covar"
echo "interactions file: $int"
echo "Output file: ${out_path}/${out_prefix}"
echo "settings: mode=cis_nominal;best_only;fdr=0.1"

python3 -m tensorqtl ${geno} ${pheno} ${out_prefix} \
    --covariates ${covar} \
    --interaction ${int} \
    --best_only \
    --fdr 0.1 \
    --mode cis_nominal \
    -o ${out_path}

echo end on "[$(date)]"
fi
done


outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/interaction/"
clusters <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt")

myDir <- outFolder #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep("tensorqtlint.cis_qtl_top_assoc.txt.gz", filenames)] #pick specific files from list
int <- ldply(lapply(clusters$cluster,function(cluster){
  #vl <- ldply(lapply(1:length(contrastdf$control),function(x){
  vl <- ldply(lapply(2,function(x){
    con=contrastdf[x,]
    v=paste0(con$treatment,"_vs_",con$control)
    cat("running",cluster,v,"\n")
    filenamesc <- filenames[grepl(paste0(cluster,".",v), filenames)] #pick specific files from list
    if(length(filenamesc)>0){
    data_names <- gsub("_tensorqtlint.cis_qtl_top_assoc.txt.gz", "", filenamesc) #remove file ending
    df <- fread(file.path(myDir, filenamesc[1]),header = T)[,analysis:=data_names[1]] #read in specific files and set the df object names. can dro0p unwanted columns
    #table(df$analysis,df$pval_adj_bh<0.1)
    #pval_emt = pval_gi * tests_emt and pval_adj_bh=p.adjust(pval_emt,method="BH")
    df <- transform(df, genotype_pval_emt=pval_g*tests_emt, variable_pval_emt=pval_i*tests_emt)
    df <- transform(df, genotype_padj=p.adjust(genotype_pval_emt,method="BH"), variable_padj=p.adjust(variable_pval_emt,method="BH"),
      cluster=cluster,treatment=con$treatment,control=con$control,contrast=v)
    df <- df %>% relocate(c(cluster,treatment,control,contrast))  # have to use ensgene as there was multi gene symbols
    #table(singelchrrun$analysis,singelchrrun$genotype_padj<0.1)
    #table(singelchrrun$analysis,singelchrrun$variable_padj<0.1)
    return(df)
    } 
  }),data.frame)
}),data.frame)
#counts_ls[sapply(counts_ls, is.null)] <- NULL

fwrite(int, paste0(outFolder,"tensorint_treatments.GxE_abundance_perclus.txt"), sep='\t', quote=F, row.names=F)

tested <- as.data.frame(table(int$cluster,int$contrast))
intsig <- subset(int,pval_adj_bh<0.1 )
tablesig <- as.data.frame(table(intsig$cluster,intsig$contrast))
outtable <- merge(tested,tablesig,by=c("Var1","Var2"))
colnames(outtable) <- c("cluster","contrast","tested","int_FDR10")
fwrite(outtable, paste0(outFolder,"tensorint_treatments.GxE_summary_perclus.txt"), sep='\t', quote=F, row.names=F)

outtablem <- melt(outtabledf)
names(outtablem)[4] <- "count"
p <- ggplot(outtablem, aes(fill=count, y=value, x=cluster)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(.~contrast,ncol=4)
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/tensorint_treatments.GxE_summary_bar.png")
png(width = 12, height = 10, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

intsig <- transform(intsig, direction=if_else(b_gi<0,"down","up"))
intsigc <- plyr::count(intsig, c("cluster","contrast","direction"))
intsigcdf <- merge(intsigc, intsig)
intsigcdf <- transform(intsigcdf, DEG_direction= if_else(direction=="down",-(freq),freq))
p <- ggplot(unique(intsigcdf[,c("contrast","cluster","DEG_direction","direction")]), aes(fill=direction, y=DEG_direction, x=cluster)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(.~contrast,ncol=4,scales="free_y")+
    labs(y="# Interaction eGenes")
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/tensorint_treatments.GxE_summary_bar_degonly.png")
png(width = 13, height = 10, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()


library(ggrastr)

lapply(split(int,int$contrast),function(v){
    
    v <- transform(v, cluster=as.factor(cluster))
    v <- v %>%
   group_by(cluster)%>%
   arrange(pval_gi) %>%
   mutate(observed=-log10(pval_gi), expected=-log10(ppoints(length(pval_gi))))

p0 <- ggplot(v, aes(x=expected, y=observed, color=cluster))+
    geom_point()+
    geom_abline(color="grey")+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    ggtitle(paste0(unique(v$contrast)))+
    theme_bw()
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/tensorint_treatments._GxE_abundance_pvalues_",unique(v$contrast),".pcl_qqplot.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()
p0 <- ggplot(v, aes(x=expected, y=observed, color=cluster))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    facet_grid(.~cluster, scales="free_y")+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    ggtitle(paste0(unique(v$contrast)))+
    theme_bw()

    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/tensorint_treatments.GxE_abundance_pvalues_",unique(v$contrast),".facetcluster_pcl_qqplot.png")
png(width = 12, height = 5, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()

})
