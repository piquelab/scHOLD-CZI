data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/"
treat="CTRL"
mkdir ${out_path}/dosages

R 
library(data.table)
library(plyr);library(dplyr)

base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.1
dimset=50
filter <- "CTRLonly" #ALOFT used

outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/"
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/results/"
best_df <- fread(paste0(outFolder,"bestPCs_table.txt"))
clusters <- fread(paste0(base,method,"_pseudobulk_ctrl/",filter,"/",project,".",resset,".",dimset,".cluster_celltype.txt"))

filenames <- list.files(resultsFolder) #file list from directory
filenames1 <- filenames[grep("_significant_topeeQTL_snps.txt", filenames)] #pick specific files from list
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

bcftools view -i'ID=@/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/results/all_eQTL_coordinates_uniq.txt' /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' > ${out_path}/results/all_uniq_eQTL_dosages.txt
bcftools view -h -i'ID=@/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/results/all_eQTL_coordinates_uniq.txt' /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.filtered.vcf.gz |sed '/^##/d' > ${out_path}/results/all_uniq_eQTL_h.txt

#per cluster dosages
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/CTRLonly/ALL.0.1.50.cluster_celltype.txt`;do
    echo running $cluster $treat
    filein=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/results/${cluster}.${treat}.PC*_significant_topeeQTL_snps_region.txt
    bcftools view -R $filein /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' > ${out_path}/results/${cluster}.${treat}_eQTL_dosages.txt
    bcftools view -h -R $filein /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz |sed '/^##/d' > ${out_path}/results/${cluster}.${treat}_eQTL_h.txt
done

#save genotype (for example plots)
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf > ${data_path}/ref.ac1.gsubRI.genotype.txt

#100125 decisiion to use MAF 0.1 as Justyna used
#per cluster dosages
bgzip /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.maf10.vcf && tabix -p vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.maf10.vcf.gz

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/CTRLonly/ALL.0.1.50.cluster_celltype.txt | egrep -v "C9|C10"`;do
    echo running $cluster $treat
    filein=/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/results/${cluster}.${treat}.PC*_significant_topeeQTL_snps_region.txt
    bcftools view -R $filein /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.maf10.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID[\t%DS]\n' > ${out_path}/results/${cluster}.${treat}.maf10_eQTL_dosages.txt
    bcftools view -h -R $filein /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.maf10.vcf.gz |sed '/^##/d' > ${out_path}/results/${cluster}.${treat}.maf10_eQTL_h.txt
done

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.maf10.vcf > ${data_path}/ref.maf10.genotype.txt
bcftools view -h /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.maf10.vcf.gz |sed '/^##/d' > ${data_path}/ref.maf10_h.txt

module swap gnu7/7.3.0 gnu9

####

R 
library(data.table)
#library(tidyverse)
library(qvalue)
library(readr)
library(plyr);library(dplyr)
library(ggplot2)
library(broom)

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
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
glmnetfolder=paste0(baseoutFolder,"glmnet/")
outFolder=paste0(glmnetfolder,normmethod,"/")
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/results/"
#countsoutFolder=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/")
#opfn <- paste0(countsoutFolder,"ALL.0.2.50.DESeq_countlists_wavefilt.icfilt.RData")
#load(opfn)

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

# 3. txt file with list of testable pairs:
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
    dosbed <- fread(file=paste0(resultsFolder,c,".",treat,".maf10_eQTL_dosages.txt"))
    samples <- fread(paste0(resultsFolder,c,".",treat,".maf10_eQTL_h.txt"),header=F)
    dosind <-  unlist(unname(as.vector(samples)))[-c(1:9)]
    names(dosbed)[c(4:length(colnames(dosbed)))] <- dosind
    dosbed <- dosbed[!duplicated(dosbed[,c("V3")]),] #no duplicated genes
    rownames(dosbed) <- dosbed$V3
    dos <- dosbed[,-c(1:2)]

    pairso <- fread(paste0(resultsFolder,filenames1[grep(paste0(c,".",treat),filenames1)]))
    #pairso <- subset(eqtlpos, cluster==c)
    #only samples from that cluster

    # 1. "bed" file with gene expression residuals:
    GEfull <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals_ctrlonly/phenotypes.",c,".",treat,".residuals_",normmethod,".bed"))
    GE_a <- GEfull[,-1] #maybe this is different in txt not bed?
    GE_a <- GEfull[,-c(1:3)] #maybe this is different in txt not bed?
    rownames(GE_a) <- GEfull$gene_id

    GEind <- colnames(GE_a)[-c(1)]
    cvind <- preds$Sample_ID
    common_samples <- Reduce(intersect, list(dosind,GEind,cvind))

    #subset all samples
    cv <- preds[preds$Sample_ID %in% common_samples,]
    doskeep <- colnames(dos) %in% c("V3",common_samples)
    dos <- dos[,..doskeep]
    gekeep <- colnames(GE_a) %in% c("gene_id",common_samples)
    GE <- GE_a[,..gekeep]

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
        # duplicate dosages that are eQTLs for several genes (and order as pairs):
        #dossub <- subset(dos, V3 %in% pairs$variant_id)
        #if(dim(dossub)[1]>0){
        #remove homozygous dosages
        dossub = dos[rowSums(dos[,-1])>1,]

        common_variants <- Reduce(intersect, list(dossub$V3,pairso$variant_id))

        pairs <- subset(pairso, variant_id %in% common_variants)
        dossub2 <- subset(dossub, V3 %in% common_variants)

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
        #if(all(expressionf$gene_id == pairs$phenotype_id)==FALSE | is.na(all(expressionf$gene_id == pairs$phenotype_id)) ){
        #    pairs <- subset(pairs, phenotype_id %in% expressionf$gene_id)
        #    expressionf <- GEsub[match(pairs$phenotype_id, GEsub$gene_id),]
        #    all(expressionf$gene_id == pairs$phenotype_id) # shoud be true
        #}
        expression1 <- expressionf[,-1]

#        dossub3 <- subset(dossub2, V3 %in% pairs$variant_id)

#        dosagesf <- dossub2[match(pairs$variant_id, dossub3$V3),]
#        all(dosagesf$V3 == pairs$variant_id) 

        # transpose expression and dosages (because lm takes columns):
        expression <- t(expression1)
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

outtablem <- melt(outtabledf)
names(outtablem)[4] <- "count"
p <- ggplot(outtablem, aes(fill=count, y=value, x=cluster)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(.~description,ncol=4)
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_summary_bar.png")
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
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_summary_bar_degonly.png")
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
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_abundance_pvalues_",unique(v$variable),".pcl_qqplot.png")
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

    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_abundance_pvalues_",unique(v$variable),".facetcluster_pcl_qqplot.png")
png(width = 12, height = 5, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()

})


cluster_celltype <- fread(paste0(base,method,"_pseudobulk_ctrl/",filter,"/",project,".",resset,".",dimset,".cluster_celltype.txt"))
cats <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/aloft_variables_categories.txt")

int <- merge(int,cluster_celltype,by="cluster")
int <- transform(int, set_combo=paste0(variable,"_",cell_type))
int_sig <- subset(int, interaction_padj<0.1)
int_sig_cat <- merge(int_sig,cats,by=c("variable"))
int_sig_cat <- subset(int_sig_cat, !category %in% c("other","puberty"))
int_sig_cat <- subset(int_sig_cat, variable %in% allvarscorr_1$variable)

#upset plot
library(UpSetR)
library(ComplexUpset)
library(tidyr)

set_list <- split(int_sig_cat[,c("phenotype_id")],int_sig_cat$set_combo)
dftest <- lapply(set_list,unlist)

binary_data <- unique(int_sig_cat[,c("phenotype_id","cell_type")]) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = cell_type, values_from = value, values_fill = 0)

png(width = 10, height = 6, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_upset_celltype.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
upset(
  binary_data,
  intersect=unique(int_sig_cat$cell_type),
  #sets = c("Env1", "Env2", "Env3"), # Specify which sets to include and their order
  #order.by = "freq",               # Order intersections by frequency (most common first)
  themes = upset_default_themes(text = element_text(size = 20)),                # Increase the font size for better readability
  min_size = 10,                    #include only intersections with a size (number of elements) greater than or equal to a certain threshold
  #mainbar.y.label = "Number of shared QTLs", # Custom label for the intersection bar plot
  #sets.x.label = "Total QTLs per environment" # Custom label for the set size bar plot
)
dev.off()

binary_data <- unique(int_sig_cat[,c("phenotype_id","description")]) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = description, values_from = value, values_fill = 0)

png(width = 14, height = 17, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_upset_variable.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
upset(
  binary_data,
  intersect=unique(int_sig_cat$description),
  #sets = c("Env1", "Env2", "Env3"), # Specify which sets to include and their order
  #order.by = "freq",               # Order intersections by frequency (most common first)
  themes = upset_default_themes(text = element_text(size = 25)),                # Increase the font size for better readability
  min_size = 8,                    #include only intersections with a size (number of elements) greater than or equal to a certain threshold
  #mainbar.y.label = "Number of shared QTLs", # Custom label for the intersection bar plot
  #sets.x.label = "Total QTLs per environment" # Custom label for the set size bar plot
)
dev.off()

binary_data <- unique(int_sig_cat[,c("phenotype_id","category")]) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = category, values_from = value, values_fill = 0)

png(width = 13, height = 10, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_upset_category.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
upset(
  binary_data,
  intersect=unique(int_sig_cat$category),
  #sets = c("Env1", "Env2", "Env3"), # Specify which sets to include and their order
  #order.by = "freq",               # Order intersections by frequency (most common first)
  themes = upset_default_themes(text = element_text(size = 20)),                # Increase the font size for better readability
  min_size = 10,                    #include only intersections with a size (number of elements) greater than or equal to a certain threshold
  #mainbar.y.label = "Number of shared QTLs", # Custom label for the intersection bar plot
  #sets.x.label = "Total QTLs per environment" # Custom label for the set size bar plot
)
dev.off()

#per variable
binary_data <- int_sig[int_sig$variable=="pedu",c("phenotype_id","cell_type")] %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = cell_type, values_from = value, values_fill = 0)

#table with categories (variable | colorbar | GxE eGenes)
int_sigc <- merge(int_sig,cats,by="variable")
int_sigc <- subset(int_sigc, variable %in% allvarscorr_1$variable)
intsigdf <- ddply(int_sigc, c("description","category"), plyr::summarize,
    'GxE eGenes'=length(unique(phenotype_id)))

intsigdf <- ddply(int_sigc, c("description","category","cell_type"), plyr::summarize,
    'GxE eGenes'=length(unique(phenotype_id)))
intsigdfc <- reshape2::dcast(intsigdf,description + category ~ cell_type,value.var="GxE eGenes")
fwrite(intsigdfc[order(intsigdfc$category),],sep='\t', quote=F, row.names=F, col.names=T,paste0(resultsFolder,"GxE_table.",resset,".",dimset,".txt"))

#overlap with asthma genes
#from wei et al single cell atac
atac_genes <- fread("/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/analysis2/4_TWAS_smr/5_pub.outs/2_supp_tables/TableS5_1_asthma-risk-genes_ALOFT.txt.gz")
#from gtex eqtl
gtex_genes <- fread("/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/analysis2/4_TWAS_smr/5_pub.outs/2_supp_tables/TableS6_1_asthma-risk-genes_gtex.txt.gz")


int_atac <- merge(int,atac_genes,by.x="phenotype_id",by.y="Gene") #FDR_intact

t <- table(GxE=int_atac$interaction_padj<0.1,ATAC=int_atac$FDR_intact<0.1, useNA = "always")[-3,-3]
df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
  CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))

           truetrue oddsratio      pval    CI_low  CI_high
odds ratio       52 0.8011669 0.1508398 0.5872311 1.073407

#overall
int_atac$interaction_pval[int_atac$interaction_pval<1e-20] <- 1e-20
int_atac <- transform(int_atac, cluster=as.factor(cluster),INTACT_sig=ifelse(FDR_intact<0.1,"INTACT_sig","INTACT_not_sig"))
v <- int_atac %>%
    group_by(INTACT_sig)%>%
    arrange(interaction_pval) %>%
    mutate(observed=-log10(interaction_pval), expected=-log10(ppoints(length(interaction_pval))))

p0 <- ggplot(v, aes(x=expected, y=observed, color=INTACT_sig))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    #facet_wrap(.~cluster, scales="free_y",ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    geom_text(
    aes(x = 0, y = Inf, label = paste0("OR=",round(t1$oddsratio,2),", p=",round(t1$pval,2))),
    hjust = -0.5, # Right-align the text
    vjust = 1.5, # Top-align the text
    size = 5,
    color = "black"
    ) +
    #ggtitle()+
    theme_bw()
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_INTACTenrich/GxE_pvalues_atac_enrich_qqplot.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()

#per variable
int_atac <- merge(int,atac_genes,by.x="phenotype_id",by.y="Gene") #FDR_intact
int_atac <- merge(int_atac,variables_df,by="variable")

lapply(split(int_atac,int_atac$variable),function(v){
    t <- table(GxE=v$interaction_padj<0.1,ATAC=v$FDR_intact<0.1, useNA = "always")[-3,-3]
    df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
    t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
      CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))
    if(t1$pval<0.1){
    cat(unique(v$description),": OR=",round(t1$oddsratio,2),", p=",round(t1$pval,2),"\n")
        v$interaction_pval[v$interaction_pval<1e-20] <- 1e-20
    v <- transform(v, cluster=as.factor(cluster),INTACT_sig=ifelse(FDR_intact<0.1,"INTACT_sig","INTACT_not_sig"))
    v <- v %>%
        group_by(INTACT_sig)%>%
        arrange(interaction_pval) %>%
        mutate(observed=-log10(interaction_pval), expected=-log10(ppoints(length(interaction_pval))))

    p0 <- ggplot(v, aes(x=expected, y=observed, color=INTACT_sig))+
        geom_point()+
        geom_abline(color="grey")+
        xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
        ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
        geom_text(
        aes(x = 0, y = Inf, label = paste0("OR=",round(t1$oddsratio,2),", p=",round(t1$pval,2))),
        hjust = -0.5, # Right-align the text
        vjust = 1.5, # Top-align the text
        size = 5,
        color = "black"
        ) +
        ggtitle(paste0(unique(v$description)))+
        theme_bw()
        figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_INTACTenrich/GxE_pvalues_",unique(v$variable),"_atac_enrich_qqplot.png")
    png(width = 8, height = 8, file=figfn, pointsize=12, 
          bg = "transparent", canvas = "white", units = "in", res = 600)
    print(p0)
    dev.off()

    lapply(split(v,v$cluster),function(c) {

    t <- table(GxE=c$interaction_padj<0.1,ATAC=c$FDR_intact<0.1, useNA = "always")[-3,-3]
    df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
    t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
      CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))
    if(t1$pval<0.1){
        c <- c %>%
            group_by(cluster,INTACT_sig) %>%
            arrange(interaction_pval) %>%
            mutate(observed=-log10(interaction_pval), expected=-log10(ppoints(length(interaction_pval))))

        p0 <- ggplot(c, aes(x=expected, y=observed, color=INTACT_sig))+
            geom_point()+
            geom_abline(color="grey")+
            xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
            ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
            geom_text(
            aes(x = 0, y = Inf, label = paste0("OR=",round(t1$oddsratio,2),", p=",round(t1$pval,2))),
            hjust = -0.5, # Right-align the text
            vjust = 1.5, # Top-align the text
            size = 5,
            color = "black"
            ) +
            ggtitle(paste0(unique(v$description)," ",unique(c$cell_type)))+
            theme_bw()
            figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_INTACTenrich/GxE_pvalues_",unique(c$variable),"_",unique(c$cluster),"_atac_enrich_qqplot_bycluster.png")
        png(width = 8, height = 8, file=figfn, pointsize=12, 
              bg = "transparent", canvas = "white", units = "in", res = 600)
        print(p0)
        dev.off()
    }
    })
    }
})

#with gtex
int_gtex <- merge(int,gtex_genes,by.x="phenotype_id",by.y="Gene") #FDR_intact

t <- table(GxE=int_gtex$interaction_padj<0.1,gtex=int_gtex$FDR_intact<0.1, useNA = "always")[-3,-3]
df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
  CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))

#overall
int_gtex$interaction_pval[int_gtex$interaction_pval<1e-20] <- 1e-20
int_gtex <- transform(int_gtex, cluster=as.factor(cluster),INTACT_sig=ifelse(FDR_intact<0.1,"INTACT_sig","INTACT_not_sig"))
v <- int_gtex %>%
    group_by(INTACT_sig)%>%
    arrange(interaction_pval) %>%
    mutate(observed=-log10(interaction_pval), expected=-log10(ppoints(length(interaction_pval))))

p0 <- ggplot(v, aes(x=expected, y=observed, color=INTACT_sig))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    #facet_wrap(.~cluster, scales="free_y",ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    geom_text(
    aes(x = 0, y = Inf, label = paste0("OR=",round(t1$oddsratio,2),", p=",round(t1$pval,2))),
    hjust = -0.5, # Right-align the text
    vjust = 1.5, # Top-align the text
    size = 5,
    color = "black"
    ) +
    #ggtitle()+
    theme_bw()
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_INTACTenrich/GxE_pvalues_gtex_enrich_qqplot.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()

#per variable
int_gtex <- merge(int,gtex_genes,by.x="phenotype_id",by.y="Gene") #FDR_intact
int_gtex <- merge(int_gtex,variables_df,by="variable")

lapply(split(int_gtex,int_gtex$variable),function(v){
    t <- table(GxE=v$interaction_padj<0.1,gtex=v$FDR_intact<0.1, useNA = "always")[-3,-3]
    df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
    t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
      CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))
    if(t1$pval<0.1){
    cat(unique(v$description),": OR=",round(t1$oddsratio,2),", p=",round(t1$pval,2),"\n")
        v$interaction_pval[v$interaction_pval<1e-20] <- 1e-20
    v <- transform(v, cluster=as.factor(cluster),INTACT_sig=ifelse(FDR_intact<0.1,"INTACT_sig","INTACT_not_sig"))
    v <- v %>%
        group_by(INTACT_sig)%>%
        arrange(interaction_pval) %>%
        mutate(observed=-log10(interaction_pval), expected=-log10(ppoints(length(interaction_pval))))

    p0 <- ggplot(v, aes(x=expected, y=observed, color=INTACT_sig))+
        geom_point()+
        geom_abline(color="grey")+
        xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
        ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
        geom_text(
        aes(x = 0, y = Inf, label = paste0("OR=",round(t1$oddsratio,2),", p=",round(t1$pval,2))),
        hjust = -0.5, # Right-align the text
        vjust = 1.5, # Top-align the text
        size = 5,
        color = "black"
        ) +
        ggtitle(paste0(unique(v$description)))+
        theme_bw()
        figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_INTACTenrich/GxE_pvalues_",unique(v$variable),"_gtex_enrich_qqplot.png")
    png(width = 8, height = 8, file=figfn, pointsize=12, 
          bg = "transparent", canvas = "white", units = "in", res = 600)
    print(p0)
    dev.off()

    lapply(split(v,v$cluster),function(c) {

    t <- table(GxE=c$interaction_padj<0.1,gtex=c$FDR_intact<0.1, useNA = "always")[-3,-3]
    df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
    t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
      CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))
    if(t1$pval<0.1){
        c <- c %>%
            group_by(cluster,INTACT_sig) %>%
            arrange(interaction_pval) %>%
            mutate(observed=-log10(interaction_pval), expected=-log10(ppoints(length(interaction_pval))))

        p0 <- ggplot(c, aes(x=expected, y=observed, color=INTACT_sig))+
            geom_point()+
            geom_abline(color="grey")+
            xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
            ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
            geom_text(
            aes(x = 0, y = Inf, label = paste0("OR=",round(t1$oddsratio,2),", p=",round(t1$pval,2))),
            hjust = -0.5, # Right-align the text
            vjust = 1.5, # Top-align the text
            size = 5,
            color = "black"
            ) +
            ggtitle(paste0(unique(v$description)," ",unique(c$cell_type)))+
            theme_bw()
            figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_INTACTenrich/GxE_pvalues_",unique(c$variable),"_",unique(c$cluster),"_gtex_enrich_qqplot_bycluster.png")
        png(width = 8, height = 8, file=figfn, pointsize=12, 
              bg = "transparent", canvas = "white", units = "in", res = 600)
        print(p0)
        dev.off()
    }
    })
    }
})


#example plots
#geno <- fread(paste0(base,"ref.ac1.gsubRI.genotype.txt"))
#vcfind <- fread(cmd = 'grep -v "^##" /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.headeronly.gsubRI.txt')

geno <- fread(paste0(base,"ref.maf10.genotype.txt"))
genoind <- fread(paste0(base,"ref.maf10_h.txt"))

names(geno) <- c("chr","pos","variant_id","ref","alt",colnames(genoind)[-c(1:9)])

int_gtexsig <- subset(int_gtex, FDR_intact<0.1 & interaction_padj<0.1)
int_gtexsig <- merge(int_gtexsig,cats,by=c("variable","description"))
int_gtexsig <- int_gtexsig[order(int_gtexsig$interaction_padj,int_gtexsig$FDR_intact),]

#decided 93025 8pm to focus on atac intersection
int_atacsig <- subset(int_atac, FDR_intact<0.1 & interaction_padj<0.1)
int_atacsig <- merge(int_atacsig,cats,by=c("variable","description"))
int_atacsig <- int_atacsig[order(int_atacsig$interaction_padj,int_atacsig$FDR_intact),]


library(edgeR)
library(stringr)

for (cat in c("asthma","psychosocial")){
    cat("running", cat, "\n")
int_atacsig_asthma <- subset(int_atacsig, category==cat)
#int_atacsig_asthma_top <- int_atacsig_asthma[, head(.SD, 1), by=symbol]
if(dim(int_atacsig_asthma)[1]>75){
    int_atacsig_asthma_top <- int_atacsig_asthma[c(1:75),]
} else{
    int_atacsig_asthma_top <- int_atacsig_asthma
}
#int_atacsig_asthma_top10 <- int_atacsig_asthma_top[c(1:20),]

#atac_genes_sig <- subset(atac_genes, FDR_intact<0.1)
#int_gtexsig_asthma_atac <- subset(int_gtexsig_asthma, symbol %in% atac_genes_sig$symbol)

#int_gtexsig1 <- plyr::count(int_gtexsig_asthma,c("symbol","cluster")) #only one variable per gene_cluster combo
#int1 <- int_gtexsig1[int_gtexsig1$freq<3,] #more unique results

gcplist <- lapply(split(int_atacsig_asthma_top,int_atacsig_asthma_top$symbol),function(g){
    #g=split(int_atacsig_asthma_top,int_atacsig_asthma_top$symbol)[[1]]
    #for (c in g$cluster){
    cplist <- lapply(g$cluster,function(c) {
        #c=g$cluster[1]
        exc <- subset(g, cluster==c)
        GEfull <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals_ctrlonly/phenotypes.",c,".",treat,".residuals_",normmethod,".bed"))
        #for(var in exc$variable){
        plist <- lapply(exc$variable, function(var){
        #var=exc$variable[1]
        ex1 <- subset(exc,variable==var)

        gene=unique(ex1$symbol)
        gene_id <- unique(ex1$phenotype_id)
        snp_id <- unique(ex1$variant_id)

        #if(!isTRUE(file.size(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/example_GxE_ATAConly/",cat,"/GxE_example_",var,"_",c,"_",gene,".png")) > 0)){
        cat("running",gene,var,c,"\n") 
#var="ctasfq"
#cluster="C0"

#combatrun="income_PCs_sex_age_adjusted" #CTRL only doesnt need treatment
#opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",c,".",combatrun,".RData")
#load(opfn)
#log2cpm <- log2(cpm(adjusted))
#sample_cols <- sapply(strsplit(colnames(log2cpm),"_"),function(y)y[3])
#log2cpm_df <- as.data.frame(log2cpm)
#colnames(log2cpm_df) <- sample_cols
#log2cpm_df$symbol <- rownames(log2cpm_df)
#log2cpm_dfm <- melt(log2cpm_df,id.vars="symbol")
#names(log2cpm_dfm)[2:3] <- c("Sample_ID","log2cpm_gene_expression")

#cv
#cluster_metadata_sce <- metadata_ls[[c]]
#cluster_metadata <- data.frame(cluster_metadata_sce)
cols <- c("Sample_ID",var)
cv_var <- preds[,..cols]

#gene expression
GEfull_ex <- subset(GEfull, gene_id==ex1$phenotype_id)

#genotype
gtbed_snp_info <- geno[,1:5]
colkeep <- colnames(geno) %in% c(colnames(gtbed_snp_info),cv_var$Sample_ID)
#genoind <- geno[,..colkeep]
gtbed <- geno[,..colkeep]
#cbind(gtbed_snp_info,genoind)

#dosages
dosbed <- fread(file=paste0(resultsFolder,c,".",treat,".maf10_eQTL_dosages.txt"))
samples <- fread(paste0(resultsFolder,c,".",treat,".maf10_eQTL_h.txt"),header=F)
ind <-  unlist(unname(as.vector(samples)))[-c(1:9)]
names(dosbed)[c(4:length(colnames(dosbed)))] <- ind
dosbed <- dosbed[!duplicated(dosbed[,c("V3")]),] #no duplicated genes --did not find any
dos <- dosbed[,-c(1:3)]

gtind <- colnames(gtbed)[-c(1:5)]
dosind <- colnames(dos)
GEind <- colnames(GEfull_ex)[-c(1:4)]
#cvind <- cluster_metadata$Sample_ID
cvind <- cv_var$Sample_ID

common_samples <- Reduce(intersect, list(gtind,dosind,GEind,cvind))

colkeep <- colnames(GEfull) %in% common_samples
GEfull_df <- GEfull_ex[,..colkeep] 
colkeep <- colnames(dos) %in% common_samples
dos <- dos[,..colkeep] 
#cluster_metadata <- subset(cluster_metadata, Sample_ID %in% common_samples)
cv <- subset(cv_var, Sample_ID %in% common_samples)
colkeep <- colnames(gtbed) %in% c(colnames(gtbed_snp_info),common_samples)
gtbed <- gtbed[,..colkeep] 

roworder <- match(colnames(gtbed)[-c(1:5)], cv$Sample_ID) 
cv <-cv[roworder,]
rownames(cv) <- cv$Sample_ID

colorder <- match(colnames(gtbed)[-c(1:5)], colnames(GEfull_df))
GEfull_df <- GEfull_df[,..colorder]
colorder <- match(colnames(gtbed)[-c(1:5)], colnames(dos))
dos <- dos[,..colorder]
rownames(dos) <- dosbed$V3

#model for slope
cv_trait <- unlist(cv[,..var])

gt_data <- gtbed[gtbed$variant_id == snp_id,]
dos_data <- dos[rownames(dos)==snp_id, ]

genotype_num <- t(gt_data[,-c(1:5)])
expression_data <- t(GEfull_df)
dosages <- t(dos_data)

model <- lm(expression_data~dosages*cv_trait)
Intercept <- summary(model)$coefficients[1,1]
dosage_beta <- summary(model)$coefficients[2,1]
metagene_beta <- summary(model)$coefficients[3,1]
interaction_beta  <- summary(model)$coefficients[4,1]
in0 <- Intercept
in1 <- Intercept+dosage_beta
in2 <- Intercept+2*dosage_beta
slop0 <- metagene_beta
slop1 <- metagene_beta+interaction_beta
slop2 <- metagene_beta+2*interaction_beta

mycolors <- c("#c25757", "#e9c61d","#3A68AE")
names(mycolors) <- c(paste0(gt_data$ref, "/", gt_data$ref), paste0(gt_data$ref, "/", gt_data$alt),paste0(gt_data$alt, "/", gt_data$alt))

df <- as.data.frame(cbind(genotype_num, expression_data, dosages))
colnames(df) <- c("genotype_num", "expression_data", "dosage")
df$expression_data <- as.numeric(df$expression_data)
df <- transform(df, genotype=ifelse(genotype_num=="0|0",paste0(unique(gt_data$ref),"/",unique(gt_data$ref)),
    if_else(genotype_num=="1|1",paste0(unique(gt_data$alt),"/",unique(gt_data$alt)),paste0(unique(gt_data$ref),"/",unique(gt_data$alt)))))
new_order <- c(paste0(unique(gt_data$ref),"/",unique(gt_data$ref)), paste0(unique(gt_data$ref),"/",unique(gt_data$alt)), paste0(unique(gt_data$alt),"/",unique(gt_data$alt)))
df$genotype <- factor(df$genotype, levels = new_order)
df$var <- cv_trait

gtcount <- plyr::count(df,"genotype")
if(!all(gtcount$freq>9)){
    cat("example doesnt have enough counts \n")
} else {

figuredir="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/example_GxE_ATAConly/"
if (!file.exists(paste0(figuredir))) dir.create(paste0(figuredir), showWarnings=F)
if (!file.exists(paste0(figuredir,cat))) dir.create(paste0(figuredir,cat), showWarnings=F)

png(paste0(figuredir,cat,"/GxE_example_",var,"_",c,"_",gene,".png"))
p <- ggplot(df, aes(x=var, y= expression_data, color = genotype)) +
    geom_point(alpha=0.7)+
    geom_abline(intercept=in0,slope=slop0,color="#c25757", size=1)+
    geom_abline(intercept=in1,slope=slop1,color="#e9c61d", size=1)+
    geom_abline(intercept=in2,slope=slop2,color="#3A68AE", size=1)+
    theme_classic() +
    scale_colour_manual(values=mycolors) +
    xlab(ex1$description) +
    ylab(paste(gene, "Gene Expression Residuals")) +
    ggtitle(paste0(unique(ex1$cell_type),"\t",gene,"\t",snp_id,"\n",
        "interaction effect= ",round(ex1$interaction_effect,2),", padj= ",round(ex1$interaction_padj,2),", TWAS z-score= ",round(ex1$zscore_twas,2))) +
    theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), axis.title.x = element_text(size = rel(0.7)))
print(p)
dev.off()
}
})
})
})
}

#table of the number of sig overlap
t <- table(GxE=int_atac$interaction_padj<0.1,ATAC=int_atac$FDR_intact<0.1, useNA = "always")[-3,-3]
df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
  CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))
t1

#venn 
library(Vennerable)
intsig_0.05 <- subset(intsig, variable %in% allvarscorr_1$variable)
d_list <- list(GxE_eGene=unique(intsig_0.05$phenotype_id),Risk_gene=unique(atac_genes[atac_genes$FDR_intact<0.1,]$Gene)) 
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_INTACTenrich/","INTACT_ATAC_GxE_venn_0.05.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
dev.off()

d_list <- list(GxE_eGene=unique(intsig_0.05$phenotype_id),Risk_gene=unique(gtex_genes[gtex_genes$FDR_intact<0.1,]$Gene)) 
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/GxE_INTACTenrich/","INTACT_GTEx_GxE_venn_0.05.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
dev.off()

##############################################


fwrite(vldfc, quote=F, sep="\t",row.names=F, file=paste0(ctrltensoroutFolder, "INTACT_ATAC_GxE_intersection_unique.txt"))



vl <- lapply(split(int_atac,int_atac$variable),function(v){
    t <- table(GxE=v$interaction_padj<0.1,ATAC=v$FDR_intact<0.1, useNA = "always")[-3,-3]
    df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
    t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
      CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))

    cl <- list.rbind(lapply(split(v,v$cluster),function(c) {
        t <- table(GxE=c$interaction_padj<0.1,ATAC=c$FDR_intact<0.1, useNA = "always")[-3,-3]
        tt <- reshape2::dcast(as.data.frame(t), GxE ~ ATAC)
        colnames(tt) <- c("","Not_risk_gene","Risk_gene")
        tt[,1] <- c("Not_GxE_eGene","GxE_eGene")
        tt <- transform(tt, variable=unique(c$variable), cell_type=unique(c$cell_type))
        return(tt)
    }))
    return(cl)
})
vldf <- list.rbind(vl)
fwrite(vldf, quote=F, sep="\t",row.names=F, file=paste0(ctrltensoroutFolder, "INTACT_ATAC_GxE_intersection_tables.txt"))

vl <- lapply(split(int_atacsig,int_atacsig$variable),function(v){
    cl <- list.rbind(lapply(split(v,v$cluster),function(c) {
    t <- data.frame(variable=unique(c$variable), cell_type=unique(c$cell_type),GxEandINTACTrisk=length(unique(c$symbol)))
        return(t)
    }))
    return(cl)
})
vldf <- list.rbind(vl)
vldfc <- reshape2::dcast(vldf, variable ~ cell_type)
vldfc[is.na(vldfc)] <- 0
fwrite(vldfc, quote=F, sep="\t",row.names=F, file=paste0(ctrltensoroutFolder, "INTACT_ATAC_GxE_intersection.txt"))





    tt1 <- reshape2::dcast(as.data.frame(t), GxE ~ ATAC)
    colnames(tt1) <- c("","Not_risk_gene","Risk_gene")
    tt1[,1] <- c("Not_GxE_eGene","GxE_eGene")
    tt1 <- transform(tt1, variable=paste0(unique(c$variable),"1"), cell_type=unique(c$cell_type))

    rbind(tt,tt1)



    lapply(split(v,v$cluster),function(c) {

    t <- table(GxE=c$interaction_padj<0.1,gtex=c$FDR_intact<0.1, useNA = "always")[-3,-3]
    df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
    t1 <- unique(data.frame(truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
      CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))


int_atacsig <- subset(int_atac, FDR_intact<0.1 & interaction_padj<0.1)
int_atacsig <- merge(int_atacsig,cats,by=c("variable","description"))




####
res_dfm <- melt(GEfull[,-c(1:3)],id.vars="gene_id")
names(res_dfm) <- c("phenotype_id","Sample_ID","gene_expression_residual")


ex1gene_geno <- merge(ex1, geno, by="variant_id") 
ex1gene_genom <- melt(ex1gene_geno[,-c(3,5:44)],id.vars=c("symbol","variant_id","phenotype_id"))
ex1gene_genom <- transform(ex1gene_genom, genotype=ifelse(value=="0|0",paste0(unique(ex1gene_geno$ref),"/",unique(ex1gene_geno$ref)),
    if_else(value=="1|1",paste0(unique(ex1gene_geno$alt),"/",unique(ex1gene_geno$alt)),paste0(unique(ex1gene_geno$ref),"/",unique(ex1gene_geno$alt)))))
new_order <- c(paste0(unique(ex1gene_geno$ref),"/",unique(ex1gene_geno$ref)), paste0(unique(ex1gene_geno$ref),"/",unique(ex1gene_geno$alt)), paste0(unique(ex1gene_geno$alt),"/",unique(ex1gene_geno$alt)))
ex1gene_genom$genotype <- factor(ex1gene_genom$genotype, levels = new_order)
names(ex1gene_genom)[4:5] <- c("Sample_ID","geno")

df <- merge(ex1gene_genom,log2cpm_dfm,by=c("symbol","Sample_ID"))
#df <- merge(ex1gene_genom,res_dfm,by=c("phenotype_id","Sample_ID"))
df2 <- merge(df,cluster_metadata[,c("Sample_ID",var)],by=c("Sample_ID"))

my_cols <- c('darkred','purple','darkblue')
#my_cols <- c('A/A'='darkred','A/G'='purple','G/G'='darkblue')

p <- ggplot(df2, aes(x = !!sym(var), y = log2cpm_gene_expression, color = genotype)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)+
  stat_cor(label.x.npc = "middle", label.y.npc = "top", method="pearson",cor.coef.name = "R", size=5, label.sep="", aes(label = ..r.label..), r.digits=2,na.rm=T)+
  scale_color_manual(values = my_cols)+
  labs(color = "Genotype",x=unique(ex1$description))+
  ggtitle(paste0(unique(ex1$cell_type),"\n",unique(df2$symbol)," ",unique(df2$variant_id)))+
  theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
  axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
  axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
  legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) 

figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/example_GxE_ATAConly/",cat,"/GxE_example_",var,"_",c,"_",gene,".png")
#png(width = 8, height = 8, file=figfn, pointsize=12, 
#      bg = "transparent", canvas = "white", units = "in", res = 600)
#print(p)
#dev.off()

#} #if plot exists
return(p)
})
return(plist)
})
#figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/example_GxE_ATAConly/",cat,"/GxE_example_",gene,".png")
#png(width = 8, height = 8, file=figfn, pointsize=12, 
#      bg = "transparent", canvas = "white", units = "in", res = 600)
#grid.arrange(grobs = cplist, ncol = 3) # Arrange in 3 columns
#dev.off()
return(cplist)
})
red.list <- lapply(rapply(gcplist, enquote, how="unlist"), eval)

figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/example_GxE_ATAConly/",cat,"/GxE_example_genes.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
grid.arrange(grobs = gcplist, ncol = 3) # Arrange in 3 columns
dev.off()
}


#check traf1
atac_genes[atac_genes$symbol=="TRAF1",]
atac_genes[atac_genes$symbol=="RUVBL1",]
atac_genes[atac_genes$symbol=="TTC19",]
atac_genes[atac_genes$symbol=="TUBB2A",]
atac_genes[atac_genes$symbol=="CD226",]
atac_genes[atac_genes$symbol=="NT5C2",]
atac_genes[atac_genes$symbol=="MAN1A2",]
atac_genes[atac_genes$symbol=="PSMD5",]
atac_genes[atac_genes$symbol=="SLFN5",]
atac_genes[atac_genes$symbol=="SNED1",]
atac_genes[atac_genes$symbol=="GLIS3",]


int_atac[int_atac$symbol=="TUBB2A" & variable=="cddstf",]






#look good
GxE_example_cdatot_C0_RNASET2.png
#RNASET2 activity is required for rAf RNASET2-induced M2 polarization of macrophages and suggests an important immune regulatory role for Af RNASET2 in ABPA pathogenesis.
#not sig in atac
GxE_example_ctasfq_C5_KDELR2.png
#Studies have identified a link between the KDELR2 gene and asthma, specifically associating it with mucus hypersecretion in the airways. The connection relates to the gene's function within the endoplasmic reticulum (ER) and its role in cellular stress responses. 
GxE_example_ctasfq_C0_KDELR2.png
https://www.ebi.ac.uk/gwas/variants/rs6796
atac_genes[atac_genes$symbol=="KDELR2"]
gtex_genes[gtex_genes$symbol=="KDELR2"]
#a positive z-score means higher gene expression is associated with increased trait values, while a negative z-score suggests higher gene expression is linked to decreased trait values
GxE_example_aBPFPM_C1_HPR.png
GxE_example_cssdh_C0_ABO.png
#some link?
