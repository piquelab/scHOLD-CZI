module swap gnu9 gnu7/7.3.0
module load bcftools/1.9


bcftools query -f '%CHROM\t%POS0\t%END\t%ID[\t%DS]\n' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.maf10.vcf.gz > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/maf10_dosages.txt
bcftools view -h /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.maf10.vcf.gz |sed '/^##/d' > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/maf10_h.txt

module swap gnu7/7.3.0 gnu9

R 
library(data.table)
#library(tidyverse)
library(qvalue)
library(readr)
library(plyr);library(dplyr)
library(ggplot2)
#in r/4.5.2 
#install.packages("valr", repos = c('https://rnabioco.r-universe.dev', 'https://cloud.r-project.org'))
library(valr)

base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.1
dimset=50
normmethod="voom" 
FDR <- 0.1
cluster="C0"
treat="CTRL"
filter <- "CTRLonly" #ALOFT used
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
glmnetfolder=paste0(baseoutFolder,"glmnet/")
outFolder=paste0(glmnetfolder,normmethod,"/")
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/results_allsnps/"
if (!file.exists(resultsFolder)) dir.create(resultsFolder, showWarnings=F)
vcfFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/"
clusters <- fread(paste0(baseoutFolder,"ALL.0.1.50.cluster_celltype.txt"))
opfn <- paste0(baseoutFolder,"ALL.",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
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

cats <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/aloft_variables_categories.txt")
catsrm <- subset(cats, category %in% c("other","puberty"))
allvarscorr <- fread(file=paste0(outFolder,project,".",resset,".",dimset,".",treat,".GLMnet-correlations.txt"))
allvarscorr_s <- subset(allvarscorr, !variable %in% c(catsrm$variable) & !variable=="")
names(allvarscorr_s)[5] <- "var_explained"
allvarscorr_1 <- subset(allvarscorr_s, var_explained>=0.05)

# load the needed files:
PCs <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/","bestPCs_table.txt"))

# 2. txt file with signatures:
preds <- fread(paste0(outFolder,"all-predictions.txt"))

int <- ldply(lapply(names(counts_ls),function(c) {
    #2 dosages from run of QTL
    dosbed <- fread(file=paste0(vcfFolder,"maf10_dosages.txt"))
    samples <- fread(paste0(vcfFolder,"maf10_h.txt"),header=F)
    dosind <-  unlist(unname(as.vector(samples)))[-c(1:9)]
    names(dosbed)[c(5:length(colnames(dosbed)))] <- dosind
    names(dosbed)[c(1:4)] <- c("chr","start","end","id")
    dosbed <- dosbed[!duplicated(dosbed[,c("id")]),] #no duplicated genes
    #dosbed$variant_id <- paste(dosbed$chr,dosbed$start,dosbed$end,sep=":")
    dos <- dosbed %>% select(-c(chr,start,end)) #%>% relocate(c(variant_id)) 

    # 1. "bed" file with gene expression residuals:
    GEfull <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals_ctrlonly/phenotypes.",cluster,".",treat,".residuals_",normmethod,".bed"))

    PCC <- subset(PCs, cluster==c)$PCs
    pairso <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/",c,"_",treat,"_PC",PCC,"_tensorqtlr.cis_qtl.txt.gz"))

    GEind <- colnames(GEfull)[-c(1:3)]
    cvind <- preds$Sample_ID
    common_samples <- Reduce(intersect, list(dosind,GEind,cvind))

    #subset all samples
    cv <- preds[preds$Sample_ID %in% common_samples,]
    doskeep <- colnames(dos) %in% c("id",common_samples)
    dos <- dos[,..doskeep]
    gekeep <- colnames(GEfull) %in% c("gene_id",common_samples)
    GE <- GEfull[,..gekeep]
    #GE <- GE %>% mutate(variant_id=paste0("chr",chr,":",start,":",end)) %>% select(-c(chr,start,end)) %>% relocate(c(variant_id)) 

    vl <- ldply(lapply(unique(allvarscorr_1$variable),function(v) {
        #v="pedu"
        cat("running ",c,v,"\n")
        cols <- c("Sample_ID",v)
        cv_var <- cv[,..cols]
        names(cv_var)[2] <- c("var")
        #this skips variables that have all the same value 
        if(length(unique(cv_var$var)) == 1){
            cat(v,"is all one value, skipping ... \n")
        } else {

        #remove homozygous dosages
        dossub = dos[rowSums(dos[,-1])>1,]

        common_variants <- Reduce(intersect, list(dossub$id,pairso$variant_id))

        pairs <- subset(pairso, variant_id %in% common_variants)
        dossub2 <- subset(dossub, id %in% common_variants)

        # subset GE df to only tested genes and sort the same way as in pairs:
        common_genes <- Reduce(intersect, list(pairs$phenotype_id,GE$gene_id))
        
        pairs <- subset(pairs, phenotype_id %in% common_genes)
        GEsub <- subset(GE, gene_id %in% common_genes)

        #match order
        colorder <- match(c("gene_id",colnames(dossub2)[-1]), colnames(GEsub))
        GEsub <- GEsub[,..colorder]

        roworder <- match(colnames(dossub2)[-1], cv_var$Sample_ID) 
        cv_var <-cv_var[roworder,]
        rownames(cv_var) <- cv_var$Sample_ID

        expressionf <- GEsub[match(pairs$phenotype_id, GEsub$gene_id),]
        cat("genes match=",all(expressionf$gene_id == pairs$phenotype_id),"\n") # 
        expression1 <- expressionf[,-1]

        # transpose expression and dosages (because lm takes columns):
        expression <- t(expression1)
        dosages <- t(dossub2[,-1])
        colnames(dosages) <- dossub2$id
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
            #if(psub %in% colnames(dosages) ){
             model <- try(lm(expression[,i]~dosages[,psub]*cv_var$var))
            #} else {
            #    model <- NA 
            #    class(model) <- c("try-error", class(model))
            #    cat("variant missing :")
            #}
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

#int <- fread(paste0(resultsFolder,treat,".GxE_abundance_perclus.txt"))
int <- subset(int, variable %in% allvarscorr_1$variable)
tested <- as.data.frame(table(int$cluster,int$variable))
intsig <- subset(int,interaction_padj<0.1 )
tablesig <- as.data.frame(table(intsig$cluster,intsig$variable))
outtable <- merge(tested,tablesig,by=c("Var1","Var2"))
colnames(outtable) <- c("cluster","variable","tested","int_FDR10")
outtabledf <- merge(outtable,variables_dftorun)
fwrite(outtabledf[,c("cluster","variable","description","tested","int_FDR10")], paste0(resultsFolder,treat,".GxE_summary_perclus.txt"), sep='\t', quote=F, row.names=F)
