R
library(data.table)
library(qvalue)
library(ggplot2)
library(plyr)
library(tidyverse)

cluster="C4"
treat="CTRL"
PC=2
ci=0.95

baseoutFolder=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,"ALL.0.2.50.DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)

outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/"
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

lapply(names(counts_ls),function(cluster){
#  lapply(treatments,function(treat){
    cat("running ",cluster,treat,"\n")
#phenotype_id    chrom   num_var variant_id      pos     tss_distance    ma_count        af      beta_shape1     beta_shape2     beta_converged  opt_status      true_nc pval_nominal    slope   slope_se        pval_beta       alpha_cov       model_converged
jaxqtlout <- fread(paste0(outFolder,cluster,"_",treat,".jaxqtl_NB.cis_score.tsv.gz"))
jaxqtlout <- transform(jaxqtlout, model_converged=as.numeric(model_converged),beta_converged=as.numeric(beta_converged),
	pval_beta=as.numeric(pval_beta), slope=as.numeric(slope),slope_se=as.numeric(slope_se), pval_nominal=as.numeric(pval_nominal))
converged <- subset(jaxqtlout,model_converged > 0 & beta_converged > 0)
converged <- transform(converged, padj=p.adjust(pval_beta, method="BH"))
table(converged$padj<0.1)
#165
pc_results_bp <- converged %>% select(phenotype_id, variant_id, pval_beta) %>% filter(!is.na(pval_beta)) %>%
            arrange(pval_beta) %>%
            mutate(r=rank(pval_beta, ties.method = "random"),
                   pexp=r/length(pval_beta),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)))

pc_results_bp$group <- "Corrected p-value"
pc_results_bp <- pc_results_bp %>% dplyr::rename(pval = pval_beta) 

png(width = 12, height = 12, file=paste0(outFolder,"figures/",cluster,".",treat,"_eGene_qqplotpermuted_pvalue.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
    p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(cluster," ",treat," eGene QQ Plot")) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()

converged$pval_beta[converged$pval_beta<1e-20] <- 1e-20
pc_results_bp <- converged %>% select(phenotype_id, variant_id, pval_beta) %>% filter(!is.na(pval_beta)) %>%
            arrange(pval_beta) %>%
            mutate(r=rank(pval_beta, ties.method = "random"),
                   pexp=r/length(pval_beta),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)))

pc_results_bp$group <- "Corrected p-value"
pc_results_bp <- pc_results_bp %>% dplyr::rename(pval = pval_beta) 

png(width = 12, height = 12, file=paste0(outFolder,"figures/",cluster,".",treat,"_eGene_qqplotpermuted_pvalue_clipped.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
    p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(cluster," ",treat," eGene QQ Plot")) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()

#fastqtl comparison
fastqtl <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/permutations/",cluster,".",treat,".PC1-",PC,".permutations_PCregress2step.eQTL.txt.gz"))
colnames(fastqtl) <- c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
fastqtl$bqval <- qvalue(fastqtl$bpval)$qvalues
res <- sum(fastqtl$bqval<0.1,na.rm =TRUE)


#############
#loop over clusters
if (!file.exists(paste0(outFolder,"results/"))) dir.create(paste0(outFolder,"results/"), showWarnings=F)

bestPCtable <- lapply(names(counts_ls),function(cluster){
#  lapply(treatments,function(treat){
    cat("running ",cluster,treat,"\n")

myDir <- outFolder #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".gene_list.*PC.*.cis_score.tsv.gz", filenames)] #pick specific files from list
filenames <- filenames[grepl(paste0(cluster,".",treat), filenames)] #pick specific files from list
data_names <- gsub(".cis_score.tsv.gz", "", filenames) #remove file ending
data_names <- paste0("chr",gsub(paste0(cluster,".",treat,"_"), "", data_names)) #remove file start
PCnames <- gsub(".*\\_", "", data_names) #remove file ending

for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = T)[,PCs:=gsub("PC","",PCnames[i])]) #read in specific files and set the df object names. can dro0p unwanted columns
#combining data
all_PCs <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(all_PCs) <- PCnames
all_PCs_r <- tapply(all_PCs, names(all_PCs), dplyr::bind_rows)

res <- ldply(lapply(all_PCs_r, function(i){
    jaxqtlout <- transform(i, model_converged=as.numeric(model_converged),beta_converged=as.numeric(beta_converged),
        pval_beta=as.numeric(pval_beta), slope=as.numeric(slope),slope_se=as.numeric(slope_se), pval_nominal=as.numeric(pval_nominal))
    converged <- subset(jaxqtlout,model_converged > 0 & beta_converged > 0)
    converged <- transform(converged, padj=p.adjust(pval_beta, method="BH"))
    df <- data.frame(PCs=as.numeric(unique(i$PCs)),converged_genes=dim(converged)[1],eGenes=sum(converged$padj<0.1,na.rm=T),perc_converged=dim(converged)[1]/dim(jaxqtlout)[1])
    return(df)
}), data.frame)
names(res)[1] <- "PCname"

best.PCs <- res[res$eGenes==max(res$eGenes),"PCs"]
best.index <- res[res$eGenes==max(res$eGenes),"PCname"]
res <- res[order(res$PCs),]
if(length(best.PCs)>1){
    best.index <- paste0(cluster,".",treat,".PC1-",1)
    res <- res[order(res$PCs),]
    best.PCs <- best.PCs[1]
}
fwrite(res, file=paste0(outFolder,"results/",cluster,".",treat,".eGenes-per-GEPCs.txt"), sep='\t', quote=F, row.names=F)

# save the best results:
all_PCs_r_best <- all_PCs_r[[best.index]]
all_PCs_r_best <- transform(all_PCs_r_best, model_converged=as.numeric(model_converged),beta_converged=as.numeric(beta_converged),
    pval_beta=as.numeric(pval_beta), slope=as.numeric(slope),slope_se=as.numeric(slope_se), pval_nominal=as.numeric(pval_nominal))
all_PCs_r_best <- subset(all_PCs_r_best,model_converged > 0 & beta_converged > 0)
all_PCs_r_best <- transform(all_PCs_r_best, padj=p.adjust(pval_beta, method="BH"))

fwrite(all_PCs_r_best, paste0(outFolder,"results/",cluster,".",treat,".best_", best.PCs, ".GEPCs.txt"), sep='\t', quote=F, row.names=F)

# subset to significant only:
pc_signif_pairs <- all_PCs_r_best[all_PCs_r_best$padj<0.1,]
pairs <- pc_signif_pairs[,c("phenotype_id","variant_id"),] #geneid and snpid
fwrite(pairs, file=paste0(outFolder,"results/",cluster,".",treat,".",best.index,"_significant_topeeQTL_pairs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
snp_region <- transform(pc_signif_pairs,chr=paste0("chr",sapply(strsplit(variant_id,":"),function(y)y[1])),pos=sapply(strsplit(variant_id,":"),function(y)y[2]))
snp_region_o <- snp_region[,c("chr","pos")] 
fwrite(snp_region_o, file=paste0(outFolder,"results/",cluster,".",treat,".",best.index,"_significant_topeeQTL_snps_region.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
##save SNP IDs
snps <- pc_signif_pairs[,"variant_id"]
fwrite(snps, file=paste0(outFolder,"results/",cluster,".",treat,".",best.index,"_significant_topeeQTL_snps.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

all_PCs_r_best$pval_beta[all_PCs_r_best$pval_beta<1e-20] <- 1e-20
pc_results_bp <- all_PCs_r_best %>% select(phenotype_id, variant_id, pval_beta) %>% filter(!is.na(pval_beta)) %>%
            arrange(pval_beta) %>%
            mutate(r=rank(pval_beta, ties.method = "random"),
                   pexp=r/length(pval_beta),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(pval_beta)-r)))

pc_results_bp$group <- "Corrected p-value"
pc_results_bp <- pc_results_bp %>% dplyr::rename(pval = pval_beta) 

png(width = 10, height = 10, file=paste0(outFolder,"figures/",cluster,".",treat,".",best.index,"_eGene_qqplotpermuted_pvalue_clipped.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
    p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(cluster," ",treat," ",best.index," eGene QQ Plot with SNP:eGenes = ",max(res$eGenes))) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()
    return(res[res$eGenes==max(res$eGenes),])
})
names(bestPCtable) <- names(counts_ls)

#######################################################3
##########################################################
#making QC tables
bestPCtableu <- ldply(bestPCtable, data.frame)
best_df <- ldply(lapply(names(counts_ls),function(c){
    cat("running", c, "\n")
    best.PCs <- subset(bestPCtableu, .id==c)$PCs
    pheno <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/counts/phenotypes.",c,".",treat,".sort.bed.gz"))
    pc_signif_pairs <- fread(paste0(outFolder,"results/",c,".",treat,".best_", best.PCs, ".GEPCs.txt"))
    df <- data.frame(cluster=c,PCs=best.PCs,numInd=length(colnames(pheno))-4,testedgenes=dim(pc_signif_pairs)[1],eGenes_10=dim(pc_signif_pairs[pc_signif_pairs$padj<0.1,])[1],eGenes_5=dim(pc_signif_pairs[pc_signif_pairs$padj<0.05,])[1])
    return(df)
}),data.frame)
fwrite(best_df, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"bestPCs_table.txt"))


########################### 
data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/"
treat="CTRL"
mkdir ${out_path}/dosages

module swap gnu9 gnu7/7.3.0
module load bcftools/1.9

#per cluster dosages
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
    echo running $cluster $treat
    filein=${out_path}/results/${cluster}.${treat}.PC*_significant_topeeQTL_snps_region.txt
    bcftools view -R $filein /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' > ${out_path}/results/${cluster}.${treat}_eQTL_dosages.txt
    bcftools view -h -R $filein /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz |sed '/^##/d' > ${out_path}/results/${cluster}.${treat}_eQTL_h.txt
done

module swap gnu7/7.3.0 gnu9


R 
library(data.table)
#library(tidyverse)
library(qvalue)
library(readr)
library(plyr);library(dplyr)
library(ggplot2)
#library(boom) #not installing for some reason! 
#https://www.tidymodels.org/learn/develop/broom/
#https://rdrr.io/cran/broom/src/R/utilities.R#sym-warn_on_subclass
#https://github.com/tidymodels/broom/blob/main/R/stats-lm.R
tidy.lm <- function(
  x,
  conf.int = FALSE,
  conf.level = 0.95,
  exponentiate = FALSE,
  ...
) {
  warn_on_subclass(x, "tidy")

  ret <- as_tibble(summary(x)$coefficients, rownames = "term")
  colnames(ret) <- c("term", "estimate", "std.error", "statistic", "p.value")

  coefs <- stats::coef(x)

  if (length(coefs) != nrow(ret)) {
    # summary(x)$coefficients misses rank deficient rows (i.e. coefs that
    # summary.lm() sets to NA), catch them here and add them back. This join is
    # costly, so only do it when necessary.
    coefs <- tibble::enframe(coefs, name = "term", value = "estimate")
    ret <- left_join(coefs, ret, by = c("term", "estimate"))
  }

  if (conf.int) {
    ci <- broom_confint_terms(x, level = conf.level)
    ret <- dplyr::left_join(ret, ci, by = "term")
  }

  if (exponentiate) {
    ret <- exponentiate(ret)
  }

  ret
}
warn_on_subclass <- function(x, tidier) {
  if (length(class(x)) > 1 && class(x)[1] != "glm") {
    subclass <- class(x)[1]
    dispatched_method <- class(x)[class(x) %in% c("glm", "lm")][1]

    rlang::warn(
      paste0(
        "The `", tidier, "()` method for objects of class `",
        subclass,
        "` is not maintained by the broom team, and is only supported through ",
        "the `",
        dispatched_method,
        "` tidier method. Please be cautious in interpreting and reporting ",
        "broom output.\n"
      ),
      .frequency = "once",
      .frequency_id = subclass
    )
  }
}

base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
normmethod="voom" 
FDR <- 0.1
cluster="C0"
treat="CTRL"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
glmnetfolder=paste0(baseoutFolder,"glmnet/")
outFolder=paste0(glmnetfolder,normmethod,"/")
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/results/"
countsoutFolder=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(countsoutFolder,"ALL.0.2.50.DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)

variables <- c("pedu", "pincme", 
                "psesl", 
                "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasg", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq", 
                "genPC1", "genPC2", "genPC3",  "Sex", "cage1", "ceth1", "cwght1", "chght1", "csex1",
                "cgpd5", "cgpd", "cbpd"
                )
 
variable_names <- c("Parental Education", "Parental Income",
                "Subjective SES", 
                "Neighborhood Stress", "Self-disclosure", 
                "Perceived responsiveness", "Parental Warmth", "YR Parent Child Conflict",  "Risky Family", "Youth Depression", 
                "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Stimulated IL-5 cortisol trt", "Stimulated IL-13 cortisol trt", "Stimulated IFNy cortisol trt",    
                "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", 
                "DD Asthma symptoms", "Nightly Asthma", "Asthma Frequency", "Asthma Severity", 
                "Nightly Awakenings",  "Sleep Duration", "Sleep quality",
                "New Genotype PC1", "New Genotype PC2", "New Genotype PC3","Sex_alph", "Age", "Ethnicity", "Weight", "Height", "sex",
                "Female Menarche Status", "Female Puberty Score", "Male puberty score"
                )
variables_df <- data.frame(variable=variables, description=variable_names)
variables_dftorun <- subset(variables_df, !variable %in% c("csex1","Sex","cage1","Wave","genPC1","genPC2","genPC3"))

# load the needed files:
# 1. "bed" file with gene expression residuals:
GEfull <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/phenotypes.",cluster,".",treat,".residuals_",normmethod,".bed"))
GE_a <- GEfull[,-1] #maybe this is different in txt not bed?
GE_a <- GEfull[,-c(1:3)] #maybe this is different in txt not bed?
rownames(GE_a) <- GEfull$gene_id

# 3. txt file with list of testable pairs:
clusters <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt")
filenames <- list.files(resultsFolder) #file list from directory
filenames1 <- filenames[grep("_significant_topeeQTL_pairs", filenames)] #pick specific files from list

# 2. txt file with signatures:
preds <- fread(paste0(outFolder,"all-predictions.txt"))

int <- ldply(lapply(names(counts_ls),function(c) {
    #2 dosages from run of QTL

    dosbed <- fread(file=paste0(base,"jaxQTL/output/results/",c,".",treat,"_eQTL_dosages.txt"))
    samples <- fread(paste0(base,"jaxQTL/output/results/",c,".",treat,"_eQTL_h.txt"),header=F)
    ind <-  unlist(unname(as.vector(samples)))[-c(1:9)]
    names(dosbed)[c(4:length(colnames(dosbed)))] <- ind
    dosbed <- dosbed[!duplicated(dosbed[,c("V3")]),] #no duplicated genes
    rownames(dosbed) <- dosbed$V3
    dos <- dosbed[,-c(1:2)]

    pairso <- fread(paste0(resultsFolder,filenames1[grep(paste0(c,".",treat),filenames1)]))
    #only samples from that cluster
    cv <- preds[preds$Sample_ID %in% ind,]

    #match up all samples
    doskeep <- colnames(dos) %in% c("V3",cv$Sample_ID)
    dos <- dos[,..doskeep]
    gekeep <- colnames(GE_a) %in% c("gene_id",cv$Sample_ID)
    GE <- GE_a[,..gekeep]
    cv <-cv[match(colnames(dos)[-1], cv$Sample_ID),]
    rownames(cv) <- cv$Sample_ID

    vl <- ldply(lapply(variables_dftorun$variable,function(v) {
        #v="pedu"
        cat("running ",c,v,"\n")
        cv <- cv[,..v]
        colnames(cv) <- c("var")
        #this skips variables that have all the same value 
        if(length(unique(cv$var)) == 1){
            cat(v,"is all one value, skipping ... \n")
        } else {
        #remove homozygous dosages
        dossub = dos[rowSums(dos[,-1])>1,]

        pairs <- subset(pairso, variant_id %in% dos$V3)
        dossub2 <- subset(dossub, V3 %in% pairs$variant_id)
#        pairs <- subset(pairs, variant_id %in% dossub2$V3)

        # subset GE df to only tested genes and sort the same way as in pairs:
        GEsub <- subset(GE, gene_id %in% pairs$phenotype_id)
        expressionf <- GEsub[match(pairs$phenotype_id, GEsub$gene_id),]
        all(expressionf$gene_id == pairs$phenotype_id) # i think this doesnt match well becuase of mismatching sources of runs
        if(all(expressionf$gene_id == pairs$phenotype_id)==FALSE | is.na(all(expressionf$gene_id == pairs$phenotype_id)) ){
            pairs <- subset(pairs, phenotype_id %in% expressionf$gene_id)
            expressionf <- GEsub[match(pairs$phenotype_id, GEsub$gene_id),]
            all(expressionf$gene_id == pairs$phenotype_id) # shoud be true
        }
        expression <- expressionf[,-1]

#        dossub3 <- subset(dossub2, V3 %in% pairs$variant_id)
#        dosagesf <- dossub2[match(pairs$variant_id, dossub3$V3),]
#        all(dosagesf$V3 == pairs$variant_id) 

        # transpose expression and dosages (because lm takes columns):
        expression <- t(expression)
        dosages <- t(dossub2[,-1])
        colnames(dosages) <- dossub2$V3
        colnames(expression) <- expressionf$gene_id 

        # make a data frame that will take the results but first save original:
        if(dim(pairs)[1]<0){
            cat("pairs is empty \n")
        }else{
        cat("running model \n")
        # loop with the best number of PCs - 0: -- I'm guessing she meant to only run this using the best fastqtl PC run
        exl <- ldply(lapply(1:ncol(expression),function(i){
            gene=colnames(expression)[i]
            psub <- pairs[pairs$phenotype_id==gene,]$variant_id
            if(psub %in% colnames(dosages) ){
             model <- try(lm(expression[,i]~dosages[,psub]*cv$var))
            } else {
                model <- NA 
                class(model) <- c("try-error", class(model))
                cat("variant missing :")
            }
          #model <- try(lm(expression[,i]~dosages[,i]*cv$var))
          if(inherits(model, "try-error")){
            cat("An error occurred in lm in",i, "row \n")
            #pairs[i,c(6:17)] <- NA
          } else {
          tmodel <- as.data.frame(tidy.lm(model))
          pvalues <- t(as.data.frame(tmodel[,5]))
          colnames(pvalues) <- c("Intercept_pval","dosage_pval","metagene_pval","interaction_pval")
          effects <- t(as.data.frame(tmodel[,2]))
          colnames(effects) <- c("Intercept_effect","dosage_effect","metagene_effect","interaction_effect")
          ses  <- t(as.data.frame(tmodel[,3]))
          colnames(ses) <- c("Intercept_SE","dosage_SE","metagene_SE","interaction_SE")
          pairs_p <- cbind(pairs[i,],pvalues,effects,ses)
        #cat("multiple test correct row", i,"\n")
        # multiple test correct:
        pairs_p$Intercept_padj <- p.adjust(pairs_p$Intercept_pval,method="BH")
        pairs_p$dosage_padj <- p.adjust(pairs_p$dosage_pval,method="BH")
        pairs_p$metagene_padj <- p.adjust(pairs_p$metagene_pval,method="BH")
        pairs_p$interaction_padj <- p.adjust(pairs_p$interaction_pval,method="BH")
        pairs_p <- transform(pairs_p, cluster=c,treat=treat,variable=v)
        pairs_p <- pairs_p %>% relocate(c(cluster,treat,variable),.after =variant_id)  # have to use ensgene as there was multi gene symbols
        return(pairs_p)
    }
    }),data.frame) #expression loop
    }     
    }
    }),data.frame) #var loop
}),data.frame) #cluster loop
fwrite(int, paste0(resultsFolder,treat,".GxE_abundance_perclus.txt"), sep='\t', quote=F, row.names=F)

tested <- as.data.frame(table(int$cluster,int$variable))
intsig <- subset(int,interaction_padj<0.1 )
tablesig <- as.data.frame(table(intsig$cluster,intsig$variable))
outtable <- merge(tested,tablesig,by=c("Var1","Var2"))
colnames(outtable) <- c("cluster","variable","tested","int_FDR10")
outtabledf <- merge(outtable,variables_dftorun)
fwrite(outtabledf[,c("cluster","variable","description","tested","int_FDR10")], paste0(resultsFolder,treat,".GxE_summary_perclus.txt"), sep='\t', quote=F, row.names=F)

outtablem <- melt(outtabledf)
names(outtablem)[4] <- "count"
p <- ggplot(outtablem, aes(fill=count, y=value, x=cluster)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(.~description,ncol=4)
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/figures/GxE_summary_bar.png")
png(width = 12, height = 10, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

intsig <- transform(intsig, direction=if_else(interaction_effect<0,"down","up"))
intsigc <- plyr::count(intsig, c("cluster","variable","direction"))
intsigcdf <- merge(intsigc, intsig)
intsigcdf <- merge(intsigcdf,variables_dftorun)
intsigcdf <- transform(intsigcdf, DEG_direction= if_else(direction=="down",-(freq),freq))
p <- ggplot(unique(intsigcdf[,c("description","cluster","DEG_direction","direction")]), aes(fill=direction, y=DEG_direction, x=cluster)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(.~description,ncol=4)+
    labs(y="# Interaction eGenes")
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/figures/GxE_summary_bar_degonly.png")
png(width = 12, height = 11, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

## make a qqplot and p-value histogram:

library(ggrastr)
i="CTRL"

lapply(split(int,int$variable),function(v){
    
    v <- transform(v, cluster=as.factor(cluster))
    v <- v %>%
   group_by(cluster)%>%
   arrange(interaction_pval) %>%
   mutate(observed=-log10(interaction_pval), expected=-log10(ppoints(length(interaction_pval))))

p0 <- ggplot(v, aes(x=expected, y=observed, color=cluster))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    #facet_grid(variable~treat, scales="free_y")+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    ggtitle(paste0(unique(v$variable)))+
    theme_bw()
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/figures/GxE_abundance_pvalues_",unique(v$variable),".pcl_qqplot.png")
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
    ggtitle(paste0(unique(v$variable)))+
    theme_bw()

    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/figures/GxE_abundance_pvalues_",unique(v$variable),".facetcluster_pcl_qqplot.png")
png(width = 12, height = 5, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()

})
