data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"

mkdir ${out_path}/dosages

R 
library(data.table)
library(plyr);library(dplyr)

outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/"
clusters <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt")
best_df <- fread(paste0(outFolder,"bestPCs_table.txt"))
filenames <- list.files(resultsFolder) #file list from directory
filenames1 <- filenames[grep("_significant_topeeQTL_snps", filenames)] #pick specific files from list
data_names <- gsub(".PC.*_significant_topeeQTL_snps.txt", "", filenames1) #remove file ending
for(i in 1:length(filenames1)) assign(data_names[i], fread(file.path(resultsFolder, filenames1[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
eqtlpos <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
u_eqtlpos <- unique(eqtlpos)

fwrite(u_eqtlpos,paste0(resultsFolder,"all_eQTL_coordinates_uniq.txt"),sep="\t",col.names=F,row.names=F, quote=F)
q()
n

module swap gnu9 gnu7/7.3.0
module load bcftools/1.9

bcftools view -i'ID=@/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/all_eQTL_coordinates_uniq.txt' /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' > ${out_path}/results/all_uniq_eQTL_dosages.txt
bcftools view -h -i'ID=@/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/all_eQTL_coordinates_uniq.txt' /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf.gz |sed '/^##/d' > ${out_path}/results/all_uniq_eQTL_h.txt

#per cluster dosages
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
    echo running $cluster $treat
    filein=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/${cluster}.${treat}.PC*_significant_topeeQTL_snps_region.txt
    bcftools view -R $filein /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' > ${out_path}/results/${cluster}.${treat}_eQTL_dosages.txt
    bcftools view -h -R $filein /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz |sed '/^##/d' > ${out_path}/results/${cluster}.${treat}_eQTL_h.txt
done
####

R 
library(data.table)
#library(tidyverse)
library(qvalue)
library(readr)
library(plyr);library(dplyr)
library(ggplot2)

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
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/"
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
shreyaGE <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/IBD-eQTL_rectum_residuals_qnorm.txt", header=T, sep="\t", stringsAsFactors=FALSE)

GEfull <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/phenotypes.",cluster,".",treat,".residuals_",normmethod,".bed"))
GE_a <- GEfull[,-1] #maybe this is different in txt not bed?
GE_a <- GEfull[,-c(1:3)] #maybe this is different in txt not bed?
rownames(GE_a) <- GEfull$gene_id


# 3. txt file with list of testable pairs:
shreyapairs <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/analysis/PC19_significant_topeeQTL_pairs.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/"
clusters <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt")
filenames <- list.files(resultsFolder) #file list from directory
filenames1 <- filenames[grep("_significant_topeeQTL_pairs", filenames)] #pick specific files from list
#filenames1 <- filenames1[grepl(paste0(cluster,".",i),filenames)]
#data_names <- gsub(".PC.*_significant_topeeQTL_pairs.txt", "", filenames1) #remove file ending
#for(i in 1:length(filenames1)) assign(data_names[i], fread(file.path(resultsFolder, filenames1[i]),header = TRUE, sep='\t')[,analysis:=data_names[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
#ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
#names(ddf) <- c(data_names)
#eqtlpos <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#eqtlpos <- transform(eqtlpos, cluster=sapply(strsplit(analysis,"[.]"),function(y) y[1] ), treat=sapply(strsplit(analysis,"[.]"),function(y) y[2] ))

# 2. txt file with signatures:
preds <- fread(paste0(outFolder,"all-predictions.txt"))

int <- ldply(lapply(names(counts_ls),function(c) {
    #2 dosages from run of QTL

    dosbed <- fread(file=paste0(base,"tensorQTL/output/results/",c,".",treat,"_eQTL_dosages.txt"))
    samples <- fread(paste0(base,"tensorQTL/output/results/",c,".",treat,"_eQTL_h.txt"),header=F)
    ind <-  unlist(unname(as.vector(samples)))[-c(1:9)]
    names(dosbed)[c(4:length(colnames(dosbed)))] <- ind
    dosbed <- dosbed[!duplicated(dosbed[,c("V3")]),] #no duplicated genes
    rownames(dosbed) <- dosbed$V3
    dos <- dosbed[,-c(1:2)]

    pairso <- fread(paste0(resultsFolder,filenames1[grep(paste0(c,".",treat),filenames1)]))
    #pairso <- subset(eqtlpos, cluster==c)
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
        # duplicate dosages that are eQTLs for several genes (and order as pairs):
        #dossub <- subset(dos, V3 %in% pairs$variant_id)
        #if(dim(dossub)[1]>0){
        #remove homozygous dosages
        dossub = dos[rowSums(dos[,-1])>1,]

#dos_pos <- paste0(gsub("chr","",dosbed$V1),":",dosbed$V2)
#pairso_pos <- paste0(sapply(strsplit(pairso$variant_id,":"),function(y) y[1]),":",sapply(strsplit(pairso$variant_id,":"),function(y) y[2]))
#eqtlrposc1 <- eqtlrpos[eqtlrpos$cluster=="C1",]
#eqtlrposc1$variant_id <- paste0(gsub("chr","",eqtlrposc1$V1),":",dosbed$V2)
#dos_pos %in% eqtlrposc1$variant_id

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
          tmodel <- as.data.frame(tidy(model))
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
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/GxE_summary_bar.png")
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
    facet_wrap(.~description,ncol=4,scales="free_y")+
    labs(y="# Interaction eGenes")
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/GxE_summary_bar_degonly.png")
png(width = 13, height = 10, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()


## make a qqplot and p-value histogram:

library(ggrastr)
i="CTRL"
#plotDF <- ldply(lapply(split(int,int$variable),function(v){
#    ntest <- nrow(v)
#    v <- v %>%
#   arrange(interaction_pval) %>%
#   mutate(observed=-log10(interaction_pval), expected=-log10(ppoints(ntest)))
#   return(v)
#}),data.frame)

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
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/GxE_abundance_pvalues_",unique(v$variable),".pcl_qqplot.png")
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

    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/GxE_abundance_pvalues_",unique(v$variable),".facetcluster_pcl_qqplot.png")
png(width = 12, height = 5, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()

})




##########stopped here###############


        # save the results:
        write.table(pairs, file=paste0("./GxE_results/GxE_abundance_19PCs_", trait, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

        # report number of significant interactions:
        pairs <- pairs %>% drop_na(interaction_qval)

        signif <- nrow(pairs[pairs$interaction_qval<FDR,])
        tab <- t(c(trait,k, signif))

        write.table(tab, file="./signif_interactions_list_uncorrected.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

        ## make a qqplot and p-value histogram:
        pdf(paste0("./plots/QQplots/GxE_abundance_19PCs_pvalues_QQplot_", trait, ".pdf"))
        library(qqman)
        qq(pairs$interaction_pval)
        dev.off()
        pdf(paste0("./plots/pval_histograms/GxE_abundance_19PCs_pvalues_histogram_", trait, ".pdf"))
        hist(pairs$interaction_pval)
        dev.off()
    })
})
######### Shreya's script ############
################################################################################

# 1/25/2019 JR
##modified by SN 4-29-24
# this file uses lm to test for genotype by microbial abundance interaction for IBD-eQTL project - GE data corrected for 17Cs for ALL taxa that we have microbial data for
library(tidyverse)
library(qvalue)
library(readr)
args = commandArgs(trailingOnly=TRUE)

trait <- "Bacteria.Actinobacteria"
k=0

if (length(args)>0){
    trait <- args[1]
    }

FDR <- 0.1

# load the needed files:
# 1. "bed" file with gene expression:
GE <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/IBD-eQTL_rectum_residuals_qnorm.txt", header=T, sep="\t", stringsAsFactors=FALSE)

cvf <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_eQTL_mapping_updated_genotypePCs_2_9_23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cvf <- as.data.frame(cvf)

cvf <- cvf %>% filter(SampleID %in% colnames(GE))
GE <- GE[,cvf$SampleID]

identical(cvf$SampleID, colnames(GE))

colnames(GE) <- cvf$SUBJECT_ID

wxs_cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD_eQTL_WXS_metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
wxs_cov <- as.data.frame(wxs_cov)

wxs_rectum <- wxs_cov %>% filter(submitted_subject_id %in% cvf$SUBJECT_ID) %>% arrange(Run)

GE <- GE[,wxs_rectum$submitted_subject_id]

# 2. txt file with dosages:
#chr pos ID sampleID1 sampleID2 ...
dosbed <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/IBD-eQTL_rectum_PC19_signif_eGene_dosages.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
dosbed <- dosbed[!duplicated(dosbed[,c("ID")]),]
rownames(dosbed) <- dosbed[,3]
dos <- dosbed[,-c(1:3)]

# 3. txt file with list of testable pairs:
pairs <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/analysis/PC19_significant_topeeQTL_pairs.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# 4. microbial abundance values:
shreyacv <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD-eQTL_Rectum_microbiome_all_taxa_covariates_6-7-23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
shreyacv <- as.data.frame(shreyacv)

#subset dosages and GE to only samples we have microbial data for

dos <- dos[,colnames(dos) %in% shreyacv$SUBJECT_ID]
GE <- GE[,colnames(GE) %in% shreyacv$SUBJECT_ID]

shreyacv <-shreyacv[match(colnames(dos), shreyacv$SUBJECT_ID),]
rownames(shreyacv) <- shreyacv$SUBJECT_ID

shreyacv <- shreyacv %>% select(trait)
colnames(shreyacv) <- c("var") #ind x var


# subset GE df to only tested genes and sort the same way as in pairs:
expression <- GE[pairs$pid,]

# duplicate dosages that are eQTLs for several genes (and order as pairs):
dosages <- dos[pairs$sid,]

# transpose expression and dosages (because lm takes columns):
expression <- t(expression) #ind x gene
dosages <- t(dosages) #ind x snp

# make a data frame that will take the results:
pairs$Intercept_pval <- NA
pairs$dosage_pval <- NA
pairs$metagene_pval <- NA
pairs$interaction_pval <- NA
pairs$Intercept_effect <- NA
pairs$dosage_effect <- NA
pairs$metagene_effect <- NA
pairs$interaction_effect <- NA
pairs$Intercept_SE <- NA
pairs$dosage_SE <- NA
pairs$metagene_SE <- NA
pairs$interaction_SE <- NA

# loop with the best number of PCs - 0:
 for (i in 1:ncol(expression)){ #for each gene
      model <- lm(expression[,i]~dosages[,i]*shreyacv[,1])
      pvalues <- summary(model)$coefficients[,4]
      pairs[i,3] <- pvalues[1]
      pairs[i,4] <- pvalues[2]
      pairs[i,5] <- pvalues[3]
      pairs[i,6] <- pvalues[4]
      effects <- summary(model)$coefficients[,1]
      pairs[i,7] <- effects[1]
      pairs[i,8] <- effects[2]
      pairs[i,9] <- effects[3]
      pairs[i,10] <- effects[4]
      ses  <- summary(model)$coefficients[,2]
      pairs[i,11] <- ses[1]
      pairs[i,12] <- ses[2]
      pairs[i,13] <- ses[3]
      pairs[i,14] <- ses[4]
      }

# multiple test correct:
pairs$Intercept_qval <- qvalue(pairs$Intercept_pval)$qvalues
pairs$dosage_qval <- qvalue(pairs$dosage_pval)$qvalues
pairs$metagene_qval <- qvalue(pairs$metagene_pval)$qvalues
pairs$interaction_qval <- qvalue(pairs$interaction_pval)$qvalues

# save the results:
write.table(pairs, file=paste0("./GxE_results/GxE_abundance_19PCs_", trait, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# report number of significant interactions:
pairs <- pairs %>% drop_na(interaction_qval)

signif <- nrow(pairs[pairs$interaction_qval<FDR,])
tab <- t(c(trait,k, signif))

write.table(tab, file="./signif_interactions_list_uncorrected.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

## make a qqplot and p-value histogram:
pdf(paste0("./plots/QQplots/GxE_abundance_19PCs_pvalues_QQplot_", trait, ".pdf"))
library(qqman)
qq(pairs$interaction_pval)
dev.off()
pdf(paste0("./plots/pval_histograms/GxE_abundance_19PCs_pvalues_histogram_", trait, ".pdf"))
hist(pairs$interaction_pval)
dev.off()