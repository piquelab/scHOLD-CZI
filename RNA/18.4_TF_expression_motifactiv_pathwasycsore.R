####
####
library(tidyverse)
library(Matrix)
## library(DESeq2)
## library(biobroom)
library(data.table)
library(motifmatchr)
library(GenomicRanges)
library(SummarizedExperiment)
##
library(cowplot)
library(RColorBrewer)
library(scales)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(ggrastr)

###
###
rm(list=ls())



setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/ATAC/")
outdir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/ATAC/pathway/"

dataDir <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

pathwaydir <- "/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/GO/pathway_scores/"

############################
### show examples
##########################


outdir2 <- "./Examples/"
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings = F, recursive = T)

####
### need provide cluster, variable, motif names  

#var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Perceived Stress")

#cl_rna <- "C4" #"C0"
#cl_atac <- "C3" #"C0"
#celltype <- "Monocyte"
#var0 <-  "ISEL_Mean" #"PSS_all_mean"  ### "ISEL_Mean"
#var00 <- var_full[var0]



cluster_rna <- c("C0")#, "C4")
cluster_atac <- c("C0")#, "C3")
celltypes <- c("CD4+ T cell")#, "Monocyte")
variable <-  c("PSS_all_mean", "ISEL_Mean")
variable_name <-c("Psychological Stress", "Social Support")


#fname = "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/TF/00_C0_DAMs_opposite_sign.txt"
#motifs <- read.table(fname, sep="\t", header=T, quote='"', comment="")
#motifs <- motifs %>% filter(sign(zval_diff) == -1)

#motifs <- c("IRF7_157", "STAT1..STAT2_36")

#pname1 <- "Antiviral mechanism by IFN-stimulated genes" #"Cytosolic sensors of pathogen-associated DNA"#, "Interferon Signaling" #"Regulation of IFNG signaling" #"Interferon Signaling"
#pname2 <- "Interleukin-6 signaling" #"ISG15 antiviral mechanism" #"Interleukin-12 family signaling" #"Interferon gamma signaling"

#dir <- "/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/cell20filt/SES_PCs_sex_age_and_treats_generem/figures/GO/enrichPath/for_motifs/"

fname = paste0(outdir, "02_TFmotif_reactome_pathway_score_lm_pvalue_ranked_na_removed_n31_uniq9.txt")
pathmotif <- read.table(fname, sep="\t", header=TRUE, quote='"', comment="")

motifs <- pathmotif$motif_name
pname <- pathmotif$PATHNAME

#motifs <- motifs[-9]

#grep("BCL6", motif2, value=T)


    for(i in 1:length(cluster_rna)) {
        #tryCatch({

        cl_rna <- cluster_rna[i]
        cl_atac <- cluster_atac[i]
        celltype <- celltypes[i]

         for(j in 1:length(variable)) {
                var0 <- variable[j]
                var00 <- variable_name[j]


            for(m in 1:length(motifs)) {
                 motif0 <- motifs[m]
                 motif00 <- unlist(str_split(gsub("_.*", "", motif0), "\\.\\."))
                 pname1 <- pname[m]

###
### motif list
#res_motif <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz", header=T)

#res_motif <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz", header=T)
res_motif <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz", header=T)

res_motif$motif_name <- res_motif$gene

motif2 <- sort(unique(res_motif$motif_name))
## motif2 <- motif2[grepl("IRF|STAT|REL", motif2)]
## motif2 <- motif2[grepl("IRF", motif2)]

res_motif$cluster_atac <- res_motif$Cluster

res_motif <- res_motif %>% filter(Cluster %in% c("C0", "C1", "C2", "C4", "C6"))


res_motif$cluster_rna <- "NA"
res_motif$cluster_rna[res_motif$Cluster == "C0"] <- "C0"
res_motif$cluster_rna[res_motif$Cluster == "C1"] <- "C1"
res_motif$cluster_rna[res_motif$Cluster == "C2"] <- "C2"
#res_motif$cluster_rna[res_motif$Cluster == "C3"] <- "C3"
res_motif$cluster_rna[res_motif$Cluster == "C4"] <- "C3"
res_motif$cluster_rna[res_motif$Cluster == "C6"] <- "C4"



res2_motif <- res_motif%>%
    filter(motif_name%in%motif2, cluster_rna==cl_rna, cluster_atac==cl_atac, psycho_variable==var0,
           #padj_odds<0.1, padj_diff<0.1)%>%
           padj_t<0.1)%>%
    arrange(padj_t)


###
### example TF
 #motif0 <- "STAT1..STAT2_36"
 #motif00 <- unlist(str_split(gsub("_.*", "", motif0), "\\.\\."))
#motif0 <- "STAT1_41" #"STAT1..STAT2_36"  ##"IRF2_3"
#motif00 <- unlist(str_split(gsub("_.*", "", motif0), "\\.\\."))

res2_motif%>%filter(motif_name==motif0)


### all DEGs 
fname =  paste(dataDir, "allres_allvars_noCombat_DESeq.txt", sep = "")
resall1 <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
res_rna <- resall1 %>% filter(padj < 0.1)

 
res <- res_rna %>% filter(var %in% c("ISEL_Mean", "PSS_all_mean"))


###
### gene expression data 
#fn_rna <- "../2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/Ali_RNA.diff.DESeq.txt.gz"
#res_rna <- fread(fn_rna, header=T, data.table=F)
 
res2_rna <- res_rna%>%dplyr::rename(gene=identifier)%>%
    filter(cluster==cl_rna, var==var0)%>%
    arrange(pvalue) ### 113 gene

res2_rna%>%filter(gene%in%motif00)


 
###
#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/Cluster_res0.07/1.3_YtX_ave.th20.clean.rds
### motif activity
fn0 <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/2_motif.outs/option_nFeature15K_cluster_res0.12/1.3_YtX_ave.th20.clean.rds"
motif_dat <- read_rds(fn0)
x0 <- as.data.frame(str_split(colnames(motif_dat), "_", simplify = T))
names(x0) <- c("cluster", "treat", "sampleID")
x0$bti <- colnames(motif_dat)
x <- x0%>%dplyr::filter(cluster==cl_atac, treat=="CTRL")
motif_dat2 <- motif_dat[, x$bti]
colnames(motif_dat2) <- x$sampleID

#/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/norm_counts/
###
#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data
### gene expression
#/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/5_normalize_RNA_data/0_normalized_data
fn0 <- paste("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/norm_counts/", cl_rna, "_CTRL_rna.normal.rds", sep="")
rna_norm <- read_rds(fn0)
x0 <- str_split(colnames(rna_norm), "_", simplify = T)
colnames(rna_norm) <- x0[,3]

motif00 <- motif00[motif00 %in% rownames(rna_norm)]
if (length(motif00) == 0) {  
    cat("non of", motif0, "is not present in", celltype, var0, " RNA data\n")
    next
    }



## load pathways score
fn0 <- paste0(pathwaydir, cl_rna, "_reactome_pathway_scores_CTRL_HOLD_normalized_counts.txt", sep="")
pathwyscore <- read.table(fn0, sep="\t", header=T, quote='"', comment="")
rownames(pathwyscore) <- pathwyscore$PATHID

## load number of genes used per pathway
fn0 <- paste0(pathwaydir, "gene_counts_", cl_rna, "_reactome_pathway_scores_CTRL_HOLD.txt", sep="")
pathway_genes <- read.table(fn0, sep="\t", header=T, quote='"', comment="")

fname <- paste0(pathwaydir, "reactomePA_all_geneIDs.txt", sep="")
pathwayname <- read.table(fname, sep="\t", header=T, quote='"', comment="")
#x0_all <- str_split(colnames(TF_score_all), "_", simplify = T)
#colnames(TF_score_all) <- x0_all[, 3]
pathwayname$pathway_name <- gsub("^Homo sapiens: ", "", pathwayname$PATHNAME)

pID1 <- pathwayname$PATHID[pathwayname$pathway_name == pname1]
#pID2 <- pathwayname$PATHID[pathwayname$pathway_name == pname2]

yscore1 <- pathwyscore %>% filter(PATHID == pID1)
        #yscore <- pathwyscore %>% filter(PATHID == "R-HSA-1280215")

yscore_df1 <- pathwyscore %>% filter(PATHID == pID1) %>%
    pivot_longer(cols = -PATHID, names_to = "bti", values_to = "pathway_score") %>%
    dplyr::select(bti, pathway_score) %>%
    mutate(sampleID = str_extract(bti, "HO\\.\\d+") %>% str_replace("\\.", "-"))


#yscore2 <- pathwyscore %>% filter(PATHID == pID2)
        #yscore <- pathwyscore %>% filter(PATHID == "R-HSA-1280215")

#yscore_df2 <- pathwyscore %>% filter(PATHID == pID2) %>%
#    pivot_longer(cols = -PATHID, names_to = "bti", values_to = "pathway_score") %>%
#    dplyr::select(bti, pathway_score) %>%
#    mutate(sampleID = str_extract(bti, "HO\\.\\d+") %>% str_replace("\\.", "-"))


pathscore1 <- as.matrix(pathwyscore[,-1])
colnames(pathscore1) <- colnames(pathscore1) %>% str_replace_all(".*(HO\\.\\d+).*", "\\1") %>% str_replace("\\.", "-")

#pathscore2 <- as.matrix(pathwyscore[,-1])
#colnames(pathscore2) <- colnames(pathscore2) %>% str_replace_all(".*(HO\\.\\d+).*", "\\1") %>% str_replace("\\.", "-")

cat("## restructured pathway scores ", motif00, "\n")

#cv2 <- merge(yscore_df, cv0, by="sampleID")


###
### covariates
#fn0 <-paste0(atacdir, "psycho_var_dir/0.1_covariates_varSel.txt")
#cv0 <- read.table(fn0, header=T, sep="\t")
cv0 <- read.table("/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/2_Differential_analysis/psycho_var_dir/0.1_covariates_varSel.txt", header=TRUE, sep="\t")
colSel <- c("sampleID", "Batch", "sex_alph", "PC1", "PC2", "age", var0)
cv <- cv0%>%dplyr::select(all_of(colSel))%>%drop_na(all_of(colSel))



###
### combine data together

#shared <- Reduce(intersect, list(colnames(motif_dat2), colnames(rna_norm), colnames(TF_score), cv$sampleID))
shared <- Reduce(intersect, list(colnames(motif_dat2), colnames(rna_norm), yscore_df1$sampleID, cv$sampleID)) #yscore_df2$sampleID, 


           #motif_regulated_score=TF_score[motif0, sampleID], 
           #motif_regulated_score_all=TF_score_all[motif0, sampleID])

cv2 <- cv%>%filter(sampleID%in%shared)%>%
    mutate(motif_activity=motif_dat2[motif0, sampleID],
           pathway_score_1=pathscore1[pID1, sampleID])#, 
           #pathway_score_2=pathscore2[pID2, sampleID])
cv2 <- cbind(cv2, t(rna_norm[motif00, cv2$sampleID, drop=F]))


cat("## merged four data completed: ", motif00, "\n")




###
### motif activity

formu <- as.formula(paste("motif_activity", "~factor(Batch) + factor(sex_alph) + PC1 + PC2 + age"))
lm0 <- lm(formu, data=cv2)

plotDF <- cv2%>%dplyr::select(all_of(c("sampleID", var0)))%>%mutate(y=residuals(lm0))
names(plotDF)[2] <- "x"
lm0 <- lm(y~x, data=plotDF)
coef_x <- summary(lm0)$coefficients["x",]

##
### annotate lm 
b0 <- round(as.numeric(coef_x[1]), digits=3)
p0 <- as.numeric(coef_x[4])
sig0 <- case_when(p0<0.0001~"***",
    p0<0.01~"**", p0<0.05~"*", .default="NS")

eq <- deparse(bquote(italic(beta)==.(b0)~.(sig0)))

eq_x <- max(plotDF$x)*0.7
eq_y <- max(plotDF$y)*0.8


p1 <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=1, color="grey")+
   geom_smooth(method=lm, se=FALSE, linewidth=0.8)+
   annotate("text", label=eq, x=eq_x, y=eq_y, parse=T, size=3)+ 
   xlab(var00)+
   ggtitle(str_wrap(paste(gsub("_.*", "", motif0), "motif activity"), width = 20))+
   theme_bw()+
   theme(plot.title = element_text(hjust=0.5, size=11),
         axis.title.y = element_blank(),
         axis.title.x=element_text(size=9),
         axis.text = element_text(size=9)) 


###
### gene expression
p2_ls <- lapply(motif00, function(gene0){
    
## gene0 <- motif00
formu <- as.formula(paste(gene0, "~factor(Batch) + factor(sex_alph) + PC1 + PC2 + age"))
lm0 <- lm(formu, data=cv2)
 
plotDF <- cv2%>%dplyr::select(all_of(c("sampleID", var0)))%>%mutate(y=residuals(lm0))
names(plotDF)[2] <- "x"
lm0 <- lm(y~x, data=plotDF)
coef_x <- summary(lm0)$coefficients["x",]

cat(gene0, coef_x, "\n")    
##
### annotate lm 
b0 <- round(as.numeric(coef_x[1]), digits=3)
p0 <- as.numeric(coef_x[4])
sig0 <- case_when(p0<0.0001~"***",
    p0<0.01~"**", p0<0.05~"*", .default="NS")

eq <- deparse(bquote(italic(beta)==.(b0)~.(sig0)))

eq_x <- max(plotDF$x)*0.7
eq_y <- max(plotDF$y)*0.8


p2 <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=1, color="grey")+
   geom_smooth(method=lm, se=FALSE, linewidth=0.8)+
   annotate("text", label=eq, x=eq_x, y=eq_y, parse=T, size=3)+     
   xlab(var00)+
   ggtitle(bquote(italic(.(gene0))~"expression"))+
   theme_bw()+
   theme(plot.title = element_text(hjust=0.5, size=11),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=9),
         axis.text.x = element_text(size=9))

 p2
})



###
### TF regulated score DEGs
formu <- as.formula(paste("pathway_score_1", "~factor(Batch) + factor(sex_alph) + PC1 + PC2 + age"))
lm0 <- lm(formu, data=cv2)
 
plotDF <- cv2%>%dplyr::select(all_of(c("sampleID", var0)))%>%mutate(y=residuals(lm0))
names(plotDF)[2] <- "x"
lm0 <- lm(y~x, data=plotDF)
coef_x <- summary(lm0)$coefficients["x",]

##
### annotate lm 
b0 <- round(as.numeric(coef_x[1]), digits=3)
p0 <- as.numeric(coef_x[4])
sig0 <- case_when(p0<0.0001~"***",
    p0<0.01~"**", p0<0.05~"*", .default="NS")
eq <- deparse(bquote(italic(beta)==.(b0)~.(sig0)))

eq_x <- max(plotDF$x)*0.7
eq_y <- max(plotDF$y)*0.8

p3 <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=1, color="grey")+
   geom_smooth(method=lm, se=FALSE, linewidth=0.8)+
   annotate("text", label=eq, x=eq_x, y=eq_y, parse=T, size=3)+     
   xlab(var00)+
   ggtitle(str_wrap(paste(gsub("_.*", "", motif0), ": ", pname1), width = 20))+
   theme_bw()+
   theme(plot.title = element_text(hjust=0.5, size=11),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=9),
         axis.text = element_text(size=9))

### skip p4, second pathway
###

pname1 <- gsub("/", "-", pname1)

###
### save figures
if ( length(p2_ls)==1){
###     
pcomb <- plot_grid(p2_ls[[1]], p1,  p3,  ncol=1) #p4,
cl_atac2 <- gsub("C", "R", cl_atac)
cl_rna2 <- gsub("C", "A", cl_rna)
figfn <- paste(outdir2, "Figure_", motif0, "_", var0, "_", cl_atac2, "_", cl_rna2, "_", celltype, pname1, ".scatter.png", sep="")
ggsave(figfn, pcomb, width=450, height=1700, units="px", dpi=240)
### 
}else{
###
pcomb <- plot_grid(p2_ls[[1]], p2_ls[[2]], p1, p3, ncol=1) #p4, 
cl_atac2 <- gsub("C", "R", cl_atac)
cl_rna2 <- gsub("C", "A", cl_rna)
figfn <- paste(outdir2, "Figure_", gsub("_.*", "", motif0), "_", var0, "_", cl_atac2, "_", cl_rna2, "_", celltype, pname1, ".scatter.png", sep="")
ggsave(figfn, pcomb, width=450, height=2100, units="px", dpi=240)    
}






        }

    }

}









### TF regulated score  
formu <- as.formula(paste("pathway_score_2", "~factor(Batch) + factor(sex_alph) + PC1 + PC2 + age"))
lm0 <- lm(formu, data=cv2)
 
plotDF <- cv2%>%dplyr::select(all_of(c("sampleID", var0)))%>%mutate(y=residuals(lm0))
names(plotDF)[2] <- "x"
lm0 <- lm(y~x, data=plotDF)
coef_x <- summary(lm0)$coefficients["x",]

##
### annotate lm 
b0 <- round(as.numeric(coef_x[1]), digits=3)
p0 <- as.numeric(coef_x[4])
sig0 <- case_when(p0<0.0001~"***",
    p0<0.01~"**", p0<0.05~"*", .default="NS")
eq <- deparse(bquote(italic(beta)==.(b0)~.(sig0)))

eq_x <- max(plotDF$x)*0.7
eq_y <- max(plotDF$y)*0.8

p4 <- ggplot(plotDF, aes(x=x, y=y))+
   geom_point(size=1, color="grey")+
   geom_smooth(method=lm, se=FALSE, linewidth=0.8)+
   annotate("text", label=eq, x=eq_x, y=eq_y, parse=T, size=3)+     
   xlab(var00)+
   ggtitle(str_wrap(paste(gsub("_.*", "", motif0), ": ", pname2), width = 20))+
   theme_bw()+
   theme(plot.title = element_text(hjust=0.5, size=11),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size=9),
         axis.text = element_text(size=9))



