####
####
library(tidyverse)
library(Matrix)
## library(DESeq2)
## library(biobroom)
library(data.table)

##
library(cowplot)
library(RColorBrewer)
library(colorspace)
library(scales)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(ggrastr)

library(dendsort, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###

rm(list=ls())

outdir2 <- "./888_pub.outs/plots_main/"
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)




#################
### Barplots 
#################
 
indir <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"

varSel <- c("ISEL_Mean", "PSS_all_mean")

fn <- paste(indir, "2_psycho_plotData.comb.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
## res <- res%>%filter(Cluster!="C9", psycho_variable%in%varSel)
cl_atac <- paste("C", c(0, 1, 2, 4, 6), sep="")
res <- res%>%filter(Cluster%in%cl_atac, psycho_variable%in%varSel)
 

var0 <- varSel[2]

res%>%filter(padj_t<0.1, psycho_variable==var0)%>%pull(gene)%>%unique()%>%length()


### CD4+
res%>%filter(padj_t<0.1, psycho_variable == var0,
     Cluster%in%c("C0", "C3", "C5", "C7"))%>%pull(gene)%>%unique()%>%length()

### CD8+
res%>%filter(padj_t<0.1, psycho_variable == var0,
     Cluster%in%c("C1"))%>%pull(gene)%>%unique()%>%length()

### NK
res%>%filter(padj_t<0.1, psycho_variable == var0,
     Cluster%in%c("C2"))%>%pull(gene)%>%unique()%>%length()

### Mono
res%>%filter(padj_t<0.1, psycho_variable == var0,
     Cluster%in%c("C4"))%>%pull(gene)%>%unique()%>%length()

### B 
res%>%filter(padj_t<0.1, psycho_variable == var0,
     Cluster%in%c("C6"))%>%pull(gene)%>%unique()%>%length()

## res%>%filter(padj_t<0.1, psycho_variable == var0,
##      Cluster%in%c("C2"))%>%pull(gene)%>%unique()%>%length()

## res%>%filter(padj_t<0.1, psycho_variable == var0,
##      Cluster%in%c("C3"))%>%pull(gene)%>%unique()%>%length()

## res%>%filter(padj_t<0.1, psycho_variable == var0,
##      Cluster%in%c("C4"))%>%pull(gene)%>%unique()%>%length()

DAMs <- res%>%filter(padj_t<0.1)%>%pull(gene)%>%unique()   ## 182 DAMs



###############################
#### heatmap map
##############################


## col2 <- hue_pal()(9) ### default ggplot color
 
varSel <- c("ISEL_Mean", "PSS_all_mean")

cl_atac <- paste("C", c(0, 1, 2, 4), sep="")

indir <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"

###
fn <- paste(indir, "2_psycho_plotData.comb.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
res <- res%>%filter(psycho_variable%in%varSel)%>%
    mutate(comb=paste(Cluster, psycho_variable, sep="_"))
##motif2 <- res%>%filter(padj_t<0.1)%>%pull(gene)%>%unique() ## 182 motifs
motif2 <- unique(res$gene)
res <- res%>%filter(Cluster%in%cl_atac)

res2 <- res%>%filter(gene%in%motif2)

MCl_name <- c("C0"="A0_R0_T CD4+", "C1"="A1_R1_T CD8+", "C2"="A2_R2_NK" , "C3"="A3_R0_T CD4+", "C4"="A4_R3_Monocyte",
              "C5"="A5_R0_T CD4+", "C6"="A6_R4_B", "C7"="A7_R0_T CD4+")

## MCl_name <- c("C0"="A0_R0_CD4+", "C1"="A1_R1_CD8+", "C2"="A2_R2_NK" , "C3"="A3_R3_Monocyte", "C4"="A4_R4_B",
##               "C5"="A5_R0_CD4+", "C6"="A6_R0_CD4+", "C7"="A7_R0_CD4+")

var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")

col_cl <- c("C0"="#FF7F00",
   "C3"=darken("#FF7F00", 0.1), "C5"=lighten("#FF7F00", 0.2), "C7"=lighten("#FF7F00",0.5),
   "C1"="#E6E600", "C2"="#4DAF4A", "C4"="#984EA3", "C6"="#D97986")

## col_cl <- c("C0"="#FF7F00", "C1"="#E6E600", "C2"="#4DAF4A",
##       "C3"="#984EA3", "C4"="#AA4B56", 
##        "C5"=darken("#FF7F00", 0.1), "C6"=lighten("#FF7F00", 0.2), "C7"=lighten("#FF7F00",0.5),  
##        "C8"="#D4B9DA", "C9"=lighten("#D97986", 0.5))

##"#D97986"



col_var <- c("ISEL_Mean"="#3A84D8", "PSS_all_mean"="#E34234")


###
### data for heatmap
mat_df <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=statistic, values_fill=NA)
mat <- mat_df%>%column_to_rownames(var="gene") 

corr_mat <- cor(mat, method="pearson")

## opfn <- paste(outdir2, "2.2_corr_psycho.rds", sep="")
## write_rds(corr_mat, file=opfn)

###
### significance

nvar <- ncol(mat)
sig_mat <- matrix(0, nvar, nvar)

for ( i in 1 : nvar){
    ##
    for ( j in i : nvar){
        ##
        zi <- mat[, i]
        zj <- mat[, j]
        pval <- cor.test(zi, zj, method = "pearson")$p.value
        is_sig <- ifelse(pval < 0.05, 1, 0)
        ##
        sig_mat[i, j] <- is_sig
        sig_mat[j, i] <- is_sig
    }    
}


###
corr2 <- corr_mat*sig_mat

df_col <- data.frame(col_name=colnames(corr2))%>%
    mutate(variable=gsub("C._", "", col_name), var2=var_full[variable],
           cluster=gsub("_.*", "", col_name),
           MCl=MCl_name[cluster], MCl_val=as.numeric(gsub(".*_R|_.*", "", MCl)))%>%
    arrange(var2, MCl_val, .locale="en")

df_col <- df_col%>%mutate(comb_new=paste(var2, MCl, sep="_"))

###
col_sort <- df_col$col_name
corr2 <- corr2[col_sort, col_sort]


###
### set arguments for plots
## var_name <- c("PSS_all_mean"="Perceived Stress", "isel"="Social Stress")
## MCl_name <- c("C0"="CD4+ T(1)", "C1"="CD8+ T", "C2"="NK" , "C3"="Monocyte", "C4"="B",
##               "C5"="CD4+ T (2)", "C6"="CD4+ T (2)", "C7"="CD4+ T (2)")  
## df_col <- df_col%>%
##     mutate(var_full = var_name[variable], cl_name=MCl_name[cluster],
##            col_name_full = paste(var_full, cl_name, sep="_"))


## col2 <- hue_pal()(9)

####
#### annotation color for row and column 

anno_df <- df_col%>%dplyr::select(variable, cluster)

col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, cluster=col_cl),
    show_annotation_name=F, show_legend=F)
row_ha <- rowAnnotation(df = anno_df, col=list(variable=col_var, cluster=col_cl),
    show_annotation_name=F, show_legend=F)


###
### heatmap 
 
mycol <- colorRamp2(seq(-1, 1, length.out=20), colorRampPalette(c("blue", "white", "red"))(20))


 
colnames(corr2) <- df_col$MCl
rownames(corr2) <- df_col$MCl
corr2[corr2==0] <- NA

corr3 <-corr_mat[col_sort, col_sort]

##
p0 <- Heatmap(corr2, col = mycol, na_col = "white",
    cluster_rows=F, cluster_columns=F,           
    show_row_names = T, show_column_names = T,
     column_names_gp = gpar(fontsize=12),
    row_names_gp=gpar(fontsize=12),
    top_annotation = col_ha,
    left_annotation = row_ha,
    show_heatmap_legend=F)
    ###, row_names_max_width=unit(6, "cm")
    ## column_names_rot = -45,
    ## heatmap_legend_param = list(title = "PCC",
    ##     at = seq(-1, 1, by=0.5),
    ##     title_gp = gpar(fontsize=12),
    ##     labels_gp = gpar(fontsize=12),
    ##     grid_width = grid::unit(0.38, "cm"),
    ##     legend_height = grid::unit(6, "cm")))
   
    ## cell_fun=function(j, i, x, y, width, height, fill){
    ##     grid.text(round(corr3[i, j], digits=3), x, y, gp=gpar(fontsize=9))
    ##     })

hlg <- Legend(col_fun=mycol, title="PCC", at=seq(-1, 1, by=0.5),
    title_gp = gpar(fontsize=12),
    labels_gp = gpar(fontsize=12),
 ##   grid_height = grid::unit(0.38, "cm"),
    legend_width = grid::unit(5, "cm"), direction="horizontal")
               
###
### save figures 
figfn <- paste(outdir2, "Figure3B_v2_psycho_corr.heatmap.pdf", sep="")
pdf(figfn, width = 6, height = 6)
p0 <- draw(p0) ##padding=unit(c(2, 2, 2, 10), "mm"))
dev.off()

### png
figfn <- paste(outdir2, "Figure3B_v2_psycho_corr.heatmap.png", sep="")
png(figfn, width = 600, height = 600, res=120)
p0 <- draw(p0) ##padding=unit(c(2, 2, 2, 10), "mm"))
dev.off()


### legend 
figfn <- paste(outdir2, "Figure3B_v2_hlegend.pdf", sep="")
pdf(figfn, width=2.3, height=0.8)
draw(hlg, x=unit(0.5, "npc"), y=unit(0.5, "npc"), just="bottomright")
dev.off()

## png
figfn <- paste(outdir2, "Figure3B_v2_hlegend.png", sep="")
png(figfn, width=250, height=80, res=120)
draw(hlg, x=unit(0.5, "npc"), y=unit(0.5, "npc"), just="bottomright")
dev.off()


## png(figfn, width = 850, height = 700, res = 120)
## p0 <- draw(p0)
## dev.off()



######################################
#### Integrate RNAseq results 
#####################################

 
outdir2 <- "./888_pub.outs/plots_main/"


fn <- "./3_regulatory.outs/1_enrich_DAM.fisher.txt"
res <- read.table(fn, header = T, sep="\t")

varSel <- c("PSS_all_mean", "ISEL_Mean")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")


###
### data for plot
plotDF <- res%>%
    mutate(var_name=var_full[psycho_variable], comb_new=paste(var_name, MCls, sep="_"),
           comb_new_sort=fct_reorder(comb_new, var_name, .desc=T))%>%
    filter(psycho_variable%in%varSel, nTF>=5)


## plotDF%>%filter(pval<0.05, odds>1)%>%dplyr::select(MCls, comb, var, log_odds, se, nDEG, nTF)%>%arrange(var)
col_cl <- c("A0_R0_T CD4+"="#FF7F00",
   "A3_R0_T CD4+"=darken("#FF7F00", 0.1), "A5_R0_T CD4+"=lighten("#FF7F00", 0.2), "A7_R0_T CD4+"=lighten("#FF7F00",0.5),
   "A1_R1_T CD8+"="#E6E600", "A2_R2_NK"="#4DAF4A", "A4_R3_Monocyte"="#984EA3", "A6_R4_B"="#D97986")

## col_cl <- c("A0_R0_CD4+"="#FF7F00", "A1_R1_CD8+"="#E6E600", "A2_R2_NK"="#4DAF4A",
##       "A3_R3_Monocyte"="#984EA3", "A4_R4_B"="#AA4B56", 
##        "A5_R0_CD4+"=darken("#FF7F00", 0.1), "A6_R0_CD4+"=lighten("#FF7F00", 0.2),
##        "A7_R0_CD4+"=lighten("#FF7F00",0.5))

col_rna <- c("R0"="#FF7F00", "R1"="#E6E600", "R2"="#4DAF4A",
       "R3"="#984EA3", "R4"="#AA4B56")



  
plotDF2 <- plotDF%>%
    dplyr::select(comb_new, MCls, comb_new_sort, rna_cluster, atac_cluster, pval, odds, lower, upper,
    log_odds, log_lower, log_upper)

plotDF2%>%arrange(comb_new_sort)

ylab2 <- gsub(".*_A", "A", plotDF2$comb_new_sort)
names(ylab2) <- plotDF2$comb_new_sort

##
p0 <- ggplot(plotDF2, aes(x = log_odds, y = comb_new_sort, color = factor(MCls)))+
    geom_errorbarh(aes(xmax =log_upper, xmin = log_lower), linewidth = 0.5, height = 0.2)+
    geom_point(shape = 19, size = 1.5)+
    scale_colour_manual(values = col_cl, guide="none")+
    geom_vline(aes(xintercept = 0), linewidth = 0.25, linetype = "dashed")+
    xlab("log odds ratio")+xlim(-1, 2)+
    scale_y_discrete(labels=ylab2)+
    theme_bw()+
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12))

###
###
figfn <- paste(outdir2, "Figure3C_enrich_DAM.forest.png", sep="")
ggsave(figfn, p0, width = 600, height = 520, units = "px", dpi=120)

figfn <- paste(outdir2, "Figure3C_enrich_DAM.forest.pdf", sep="")
ggsave(figfn, p0, width = 4.5, height = 5)
 


###
### marginal barplots of DAM
### data for bar plots 
plotDF2 <- plotDF%>%
    dplyr::select(comb_new, MCls, comb_new_sort, rna_cluster, atac_cluster, nTF, nDEG)


###
### barplots of DEGs 

p2 <- ggplot(plotDF2, aes(x = comb_new_sort, y = nDEG, fill = factor(rna_cluster)))+
    geom_bar(stat="identity")+
    coord_flip()+
    scale_fill_manual(values=col_rna, guide="none")+
    scale_y_continuous("#DEGs")+
    theme_bw()+
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12),
          axis.text = element_text(size=12))

###
p3 <- ggplot(plotDF2, aes(x=comb_new_sort, y=nTF, fill=factor(MCls)))+
    geom_bar(stat="identity")+
    coord_flip()+
    scale_fill_manual(values=col_cl, guide="none")+
    scale_y_continuous("#DAMs")+
    theme_bw()+
    theme(axis.title.y = element_blank(),
          axis.text.y =  element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12),
          axis.text = element_text(size=12))

##
## figfn <- paste(outdir2, "Figure3C_DAM.bar.png", sep="")
## ggsave(figfn, p2, width = 280, height = 520, units = "px", dpi=120)
 


##
## figfn <- paste(outdir2, "Figure3C_DEG.bar.png", sep="")
## ggsave(figfn, p3, width = 280, height = 520, units = "px", dpi=120)

figfn <- paste(outdir2, "Figure3C_DEG_DAMs.pdf", sep="")
ggsave(figfn, plot_grid(p2, p3), width=3.5, height=5)


### png
figfn <- paste(outdir2, "Figure3C_DEG_DAMs.png", sep="")
ggsave(figfn, plot_grid(p2, p3), width=420, height=600, units="px", dpi=120)


## x <- plotDF%>%dplyr::select(


#######################################################
#### Figure 4, show heatmpa of 59 example motifs
#######################################################

###
### from 3_regulatory2.R script

rm(list=ls())


###
### read data
 
outdir2 <- "./888_pub.outs/plots_main/"
 

varSel <- c("ISEL_Mean", "PSS_all_mean")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress") 
###
### cluster  
fn0 <- "./3_regulatory.outs/TableS_nDEGs_nTFs_enriched.tsv" ##list_cluster_var_nDEGs0_nTF1_new.txt 
df_cl <- read_tsv(fn0)
df_cl <- df_cl%>%filter(psycho_variable%in%varSel, nDEG>0, nTF>=5)%>%arrange(psycho_variable)

## fn0 <- "./3_regulatory.outs/list_cluster_var_nDEGs_nTF_new.txt"
## x <- read.table(fn0, header=T, sep="\t")

###
### 
fn_enrich <- "./3_regulatory.outs/enrich_motif/2_comb_enrich.motif.txt.gz"
res_enrich <- fread(fn_enrich, header=T, sep="\t")
## res_enrich <- res_enrich%>%filter(psycho_variable%in%varSel)

###
fn2 <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
res <- fread(fn2, header=T, data.table=F, sep="\t")
##res <- res%>%filter(Cluster!="C9")
res <- res%>%filter(psycho_variable%in%varSel) ##, Cluster!="C9")
 

###
#### select 173 motifs, significantly enriched odds>1&padj<0.1,
th <- 0.1
motif_sel <- lapply(1:nrow(df_cl), function(i){
   ## 
   cl_rna <- gsub("^R", "C", df_cl$rna_cluster[i])
   cl_atac <- gsub("^A", "C", df_cl$atac_cluster[i])
    var0 <- df_cl$psycho_variable[i]

   ### 
   enrich_motif <- res_enrich%>%
       filter(cluster_rna==cl_rna, psycho_variable==var0, cluster_atac==cl_atac, odds>1, padj<th)%>%
       pull(motif_name)%>%unique()    
   DAMs <- res%>%filter(Cluster==cl_atac, psycho_variable==var0, padj_t<th)%>%pull(gene)%>%unique()

   motif2 <- intersect(enrich_motif, DAMs)
   motif2
})    

motif_sel <- unique(unlist(motif_sel))  ### 173 motifs (update) ## 5% 78 motifs

 
## motif_sel <- motif_sel[grepl("IRF|STAT|REL", motif_sel)] ### 14 motif
 

###
### make 78 motif table 
res_all <- map_dfr(1:nrow(df_cl), function(i){
    ##
    cl_rna <- gsub("^R", "C", df_cl$rna_cluster[i])
    cl_atac <- gsub("^A", "C", df_cl$atac_cluster[i])
    var0 <- df_cl$psycho_variable[i]

   ### 
   res0_enrich <- res_enrich%>%
       filter(cluster_rna==cl_rna, cluster_atac==cl_atac, psycho_variable==var0)%>%
       dplyr::select(cluster_rna, psycho_variable, motif_name, odds,
                    pval_odds=pval, padj_odds=padj, log_odds, log_lower, log_upper)
   ### 
   res0 <- res%>%filter(Cluster==cl_atac, psycho_variable==var0)%>%
        dplyr::select(cluster_atac=Cluster,
                      motif_name=gene, beta_diff=estimate, std_diff=stderror, zval_diff=statistic,
                      pval_diff=pval_t, padj_diff=padj_t)
   res0 <- res0_enrich%>%inner_join(res0, by="motif_name")
   res0$MCls <- df_cl$MCls[i]
   cat(nrow(res0), "\n")
   res0 
})


###
### save 
## opfn <- paste(outdir2, "TableS_78motif_all.infor.txt", sep="")
## res_sel <- res_all%>%filter(motif_name%in%motif_sel)
## write.table(res_sel, file=opfn, quote=F, sep="\t", row.names=F)

###
### save
res_sel <- res_all%>%filter(odds>1, padj_odds<0.1, padj_diff<0.1)
res2 <- res_sel%>%
    dplyr::select(cluster_atac, cluster_rna, MCls, psycho_variable, motif_name,
    beta_diff, std_diff, zval_diff, pval_diff, padj_diff, odds, pval_odds, padj_odds, log_odds)                            

opfn <- "./888_pub.outs/supp/TableS_sigMotif_DAM_enrich.infor.txt"
write.table(res2, file=opfn, quote=F, sep="\t", row.names=F)




###
###


#################
### dot plots
#################
## MCls_unq <- sort(unique(df_cl$MCls))
## nMCls <- length(MCls_unq)
## df_cl_new <- data.frame(MCls=rep(MCls_unq, times=2))%>%
##       mutate(rna_cluster=gsub("_.*", "", gsub(".*_R", "R", MCls)),
##       atac_cluster=gsub("_R.*", "", MCls),
##       psycho_variable=rep(c("ISEL_Mean", "PSS_all_mean"), each=nMCls))

plotDF <- map_dfr(1:nrow(df_cl), function(i){
   ###
   ###
   cl_rna <- gsub("^R", "C", df_cl$rna_cluster[i])
   cl_atac <- gsub("^A", "C", df_cl$atac_cluster[i])
   var0 <- df_cl$psycho_variable[i]

   ### 
   df0 <- res_enrich%>%
       filter(cluster_rna==cl_rna, cluster_atac==cl_atac, psycho_variable==var0, motif_name%in%motif_sel)%>%
       dplyr::select(motif_name, log_odds, log_odds_se=se, pval_enrich=pval, padj_enrich=padj)

   df2 <- res%>%filter(Cluster==cl_atac, psycho_variable==var0, gene%in%motif_sel)%>%
       dplyr::select(motif_name=gene, beta_eff=estimate, beta_se=stderror, zval=statistic, pval=pval_t, padj=padj_t)

   ##
   ## gene2 <- gsub(motif_sel 
   ## df3 <- res%>%filter(cluster==cl_rna, var==var0, identifier%in% 

   ##
   df0 <- df0%>%left_join(df2, by="motif_name")
   df0$rna_cluster <- df_cl$rna_cluster[i]
   df0$atac_cluster <- df_cl$atac_cluster[i]
   df0$MCls <- df_cl$MCls[i] 
   df0$psycho_var <- var0

   ###
   ###
   df0
})


##
opfn <- paste(outdir2, "Figure4_plot.data.rds", sep="")
write_rds(plotDF, file=opfn)


#############################################
### Heatmap for Effects on motif activity
###############################################
  
plotDF <- read_rds("./888_pub.outs/plots_main/Figure4_plot.data.rds") 
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")
plotDF2 <- plotDF%>%
    mutate(is_sig=as.integer(padj<0.1), zval2 = zval * is_sig,
           var2 = var_full[psycho_var], comb = paste(var2, MCls, sep = "_"))



###
### version-1, heatmap show all values 
mat <- plotDF2%>%
    pivot_wider(id_cols=motif_name, names_from=comb, values_from=zval, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()

rnSel <- gsub("_.*", "", rownames(mat))
motifSel_unq <- rownames(mat)[!duplicated(rnSel)]
mat <- mat[motifSel_unq, ] 

rnSel <- gsub("_.*", "", rownames(mat))
rownames(mat) <- rnSel

rn2 <- sort(unique(rownames(mat)))
comb <- sort(colnames(mat))
mat2 <-  mat[rn2, comb]  


###
### get clusters and it needs run in the first time. Only run once.  
hmap <- Heatmap(mat2, cluster_rows = T, cluster_columns = F, row_split = 4)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)

## ## row_dend <- dendsort(hclust(dist(mat2)))

df_cluster <- NULL
for ( i in 1:length(cl)){
    ##
    df0 <- data.frame(cluster = paste("cluster", i, sep=""),
        motif_name = rownames(mat2)[cl[[i]]])
    df_cluster <- rbind(df_cluster, df0)
}

opfn <- paste(outdir2, "Figure4A_heatmap_row_cluster.tsv", sep="")
write_tsv(df_cluster, file = opfn)


###########################################################
###  heatmap show significant values 10% fdr
###########################################################



###
fn <- paste(outdir2, "Figure4A_heatmap_row_cluster.tsv", sep="")
df_cluster <- read_tsv(fn, show_col_types = F)


###
mat3 <- plotDF2%>%
    pivot_wider(id_cols=motif_name, names_from=comb, values_from=zval2, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()

mat3 <- mat3[motifSel_unq, ]
rnSel <- gsub("_.*", "", rownames(mat3))
rownames(mat3) <- rnSel

rn2 <- sort(rownames(mat3))
comb <- sort(colnames(mat3))
mat3 <-  mat3[rn2, comb]

mat3 <- mat3[df_cluster$motif_name, ]



###
### column annotation


## MCl_name <- c("C0"="A0_R0_CD4+", "C1"="A1_R1_CD8+", "C2"="A2_R2_NK" , "C3"="A3_R3_Monocyte", "C4"="A4_R4_B",
##               "C5"="A5_R0_CD4+", "C6"="A6_R0_CD4+", "C7"="A7_R0_CD4+")

var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")

col_MCls <- c("A0_R0_T CD4+"="#FF7F00", "A3_R0_T CD4+"=darken("#FF7F00", 0.1), 
   "A1_R1_T CD8+"="#E6E600", "A2_R2_NK"="#4DAF4A", "A4_R3_Monocyte"="#984EA3")
   

col_var <- c("Social Support"="#3A84D8", "Psychological Stress"="#E34234")

## col_var <- c("Social Support"="#3A84D8", "Psychological Stress"="#E34234")
## col_MCls <- c("R0"="#E60A0A", "R1"="#F9A602", "R2"="#A65192", "R3"="#0570B0", "R4"="#00AF72", "R5"="#87CEFA")


anno_df <- data.frame(variable=gsub("_A.*", "", colnames(mat3)),
    cluster=gsub(".*_A", "A", colnames(mat3)))

col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, cluster=col_MCls),
    show_annotation_name = F, show_legend = F)                  



##
## hight row names  
 row_hlight <- c("BCL6", "IRF2", "IRF4", "IRF7", "STAT1", "STAT1..STAT2")
## geneSel <- read.table("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/DEGs_TF/motifs/TF-genes_expressed_in_RNA.txt", header=T)
## row_hlight <- gsub("_.*", "", geneSel$motif_name)


## iwhich <- which(gsub("_.*", "", rownames(mat3))%in%row_hlight)
## fontcolors <- rep("black", nrow(mat3))
## fontcolors[iwhich] <- "black"

fontfaces <- rep("plain", nrow(mat3))
## fontfaces[iwhich] <- "bold"


####
### dendrograms    
row_dend <- as.dendrogram(hclust(dist(mat2[rownames(mat3),])))

##just need run one time
opfn <- paste(outdir2, "Figure4A_rowname_reorder.txt", sep="")
write.table(labels(row_dend), opfn, quote=F, row.names=F, col.names=F)

## [rownames(mat3),]
## dendrograms from all zvalue
## library(dendsort)
## row_dend <- dendsort(hclust(dist(mat2[rownames(mat3), ])))


###
### row annotation  
anno_row <- data.frame(motif_name = rownames(mat3))%>%
    left_join(df_cluster, by="motif_name")%>%
    dplyr::select(cluster)
col_cluster <- c("cluster1"="#1b9e77", "cluster2"="#d95f02", "cluster3"="#e7298a", "cluster4"="#7570b3")
## row_ha <- rowAnnotation(rows = anno_text(rownames(mat3), gp=gpar(fontsize=9, fontface=fontfaces, col = fontcolors)),
##     df = anno_row, col=list(cluster=col_cluster), simple_anno_size=unit(0.7, "cm"),
##     show_annotation_name=F, show_legend=F)

row_ha <- rowAnnotation(rows = anno_text(rownames(mat3), gp=gpar(fontsize=9)),
    df = anno_row, col=list(cluster=col_cluster), simple_anno_size=unit(0.7, "cm"),
    show_annotation_name=F, show_legend=F)



###
### setting color
max_val <- max(abs(max(mat3)), abs(min(mat3)))
max_val <- round(max_val, 2)
min_val <- -max_val
mycol <- colorRamp2(seq(min_val, max_val, length.out=20), colorRampPalette(c("blue", "white", "red"))(20))


###
### make heatmap
colnames(mat3) <- gsub(".*_A", "A", colnames(mat3))
## rownames(mat3) <- gsub("_.*", "", rownames(mat3))

p1 <- Heatmap(mat3, col=mycol,
   cluster_rows=row_dend, row_order=rownames(mat3), row_dend_reorder = FALSE, 
   cluster_columns=F,
   column_names_rot = -45, column_names_gp = gpar(fontsize=12), 
   show_row_names = FALSE, ##row_names_gp = gpar(fontsize=9),
   top_annotation = col_ha, right_annotation = row_ha,  
   heatmap_legend_param = list(title = bquote(~italic(Z)~"-score"),
      at = round(seq(min_val, max_val, length.out=4), 1),
      title_gp = gpar(fontsize=12),
      labels_gp = gpar(fontsize=12),      
      grid_width = grid::unit(0.5, "cm"), legend_height = grid::unit(6, "cm")))

###
### save png
figfn <- paste(outdir2, "Figure4A_zscore_DAM.heatmap.png", sep="")
png(figfn, width = 580, height=1300, res = 120)
set.seed(0)
p1 <- draw(p1)
dev.off()
 
figfn <- paste(outdir2, "Figure4A_zscore_DAM.heatmap.pdf", sep="")
pdf(figfn, width = 5.8, height=10.5)
set.seed(0)
p1 <- draw(p1)
dev.off()


## x <- read.table("./3_regulatory.outs/correct_excludeX/motif_list_from_Ali.txt", header=T, sep="\t")



###########################
### log odds ratio
###########################

var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")
## i_inf <- is.infinite((plotDF$log_odds))
## max_val <- max(plotDF$log_odds[!i_inf])
## plotDF$log_odds[i_inf] <- max_val

plotDF2 <- plotDF%>%
    mutate(is_sig=as.integer(padj_enrich<0.1), log_odds=ifelse(is_sig, log_odds, 0), 
           var2 = var_full[psycho_var], comb = paste(var2, MCls, sep = "_"))

 
mat <- plotDF2%>%
    pivot_wider(id_cols=motif_name, names_from=comb, values_from=log_odds, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()
mat <- mat[motifSel_unq,]

rnSel <- gsub("_.*", "", rownames(mat))
rownames(mat) <- rnSel

rn2 <- sort(rnSel)
comb <- sort(colnames(mat))
mat2 <- mat[rn2, comb]


### 
### get the cluster and sort 
fn <- "./888_pub.outs/plots_main/Figure4A_heatmap_row_cluster.tsv"
df_cluster <- read_tsv(fn, show_col_types = F)

fn0 <- "./888_pub.outs/plots_main/Figure4A_rowname_reorder.txt" 
motif_sort <- read.table(fn0, header=F)$V1



mat2 <- mat2[motif_sort, ]




##
## hight row names  
## row_hlight <- c("BCL6", "IRF2", "IRF4", "IRF7", "STAT1", "STAT1..STAT2")
## geneSel <- read.table("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/DEGs_TF/motifs/TF-genes_expressed_in_RNA.txt", header=T)
## row_hlight <- gsub("_.*", "", geneSel$motif_name)


## iwhich <- which(rownames(mat2)%in%row_hlight)
## fontcolors <- rep("black", nrow(mat2))
## fontcolors[iwhich] <- "black"

fontfaces <- rep("plain", nrow(mat2))
## fontfaces[iwhich] <- "bold"


###
### row annotation  
anno_row <- data.frame(motif_name = rownames(mat2))%>%
    left_join(df_cluster, by="motif_name")%>%
    dplyr::select(cluster)
col_cluster <- c("cluster1"="#1b9e77", "cluster2"="#d95f02", "cluster3"="#e7298a", "cluster4"="#7570b3")
## row_ha <- rowAnnotation(rows = anno_text(rownames(mat2), gp=gpar(fontsize=9, fontface=fontfaces, col = fontcolors)),
##     df = anno_row, col=list(cluster=col_cluster),  simple_anno_size=unit(0.7, "cm"),
##     show_annotation_name=F, show_legend=F)

row_ha <- rowAnnotation(rows = anno_text(rownames(mat2), gp=gpar(fontsize=9)),
    df = anno_row, col=list(cluster=col_cluster),  simple_anno_size=unit(0.7, "cm"),
    show_annotation_name=F, show_legend=F)


###
### row annotation 
## anno_row <- data.frame(motif_name = rownames(mat2))%>%
##     left_join(df_cluster, by="motif_name")%>%
##     dplyr::select(cluster)
## col_cluster <- c("cluster1"="#1b9e77", "cluster2"="#d95f02", "cluster3"="#e7298a", "cluster4"="#7570b3")
## row_ha <- rowAnnotation(df = anno_row, col=list(cluster=col_cluster),
##     show_annotation_name=F, show_legend=F)


###
### setting color
max_odds <- max(plotDF$log_odds)
mycol <- colorRamp2(seq(0, max_odds, length.out=20),
                    colorRampPalette(c("white", brewer.pal(n=7, name="RdPu")[2:7]))(20))

###
### column annotation


col_MCls <- c("A0_R0_T CD4+"="#FF7F00", "A3_R0_T CD4+"=darken("#FF7F00", 0.1), 
   "A1_R1_T CD8+"="#E6E600", "A2_R2_NK"="#4DAF4A", "A4_R3_Monocyte"="#984EA3")

## col_MCls <- c("A0_R0_CD4+"="#FF7F00", "A1_R1_CD8+"="#E6E600", "A2_R2_NK"="#4DAF4A",
##       "A3_R3_Monocyte"="#984EA3", "A4_R4_B"="#AA4B56", 
##        "A5_R0_CD4+"=darken("#FF7F00", 0.1), "A6_R0_CD4+"=lighten("#FF7F00", 0.2),
##        "A7_R0_CD4+"=lighten("#FF7F00",0.5))

col_var <- c("Social Support"="#3A84D8", "Psychological Stress"="#E34234")

## col_var <- c("Social Support"="#3A84D8", "Psychological Stress"="#E34234")
## col_cl <- c("R0"="#E60A0A", "R1"="#F9A602", "R2"="#A65192", "R3"="#0570B0", "R4"="#00AF72", "R5"="#87CEFA")

anno_df <- data.frame(variable=gsub("_A.*", "", colnames(mat2)),
    cluster=gsub(".*_A", "A", colnames(mat2)))

col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, cluster=col_MCls),
    show_annotation_name = F, show_legend = F)                        


###
### make a heatmap
colnames(mat2) <- gsub(".*_A", "A", colnames(mat2))
p2 <- Heatmap(mat2, col=mycol, cluster_rows=F, cluster_columns=F, na_col="white",
   column_names_rot = -45, column_names_gp = gpar(fontsize=12),
   show_row_names = FALSE,
   ## row_names_gp = gpar(fontsize=9),
   top_annotation = col_ha, right_annotation = row_ha, 
   heatmap_legend_param = list(title = bquote(log~"odds ratio"),
      at = round(seq(0, max_odds, length.out=3)),
      title_gp = gpar(fontsize=12),
      labels_gp = gpar(fontsize=12),      
      grid_width = grid::unit(0.5, "cm"), legend_height = grid::unit(6, "cm")))

###
### save 
figfn <- paste(outdir2, "Figure4B_odds.heatmap.png", sep="")
png(figfn, width = 580, height=1300, res = 120)
p2 <- draw(p2)
dev.off()

  
figfn <- paste(outdir2, "Figure4B_odds.heatmap.pdf", sep="")
pdf(figfn, width = 5.8, height=10.5)
p2 <- draw(p2)
dev.off()



###
###
## df0 <- data.frame(x=c(-27, 26, 39))%>%mutate(y=x*1.8+32)

## p <- ggplot(df0, aes(x, y))+
##     geom_point(color="green", shape=2, size=3)+
##     geom_text(aes(x=x, y=y, label=y, vjust=-0.5), size=3)+
##     geom_abline(slope=1.8, intercept=32, color="grey", linewidth=1, linetype = "dashed")+
##     xlab(bquote("Temperatue: "~degree * C))+ xlim(-30, 50)+
##     ylab(bquote("Temperature: "~degree * F))+ylim(-20, 120)+    
##     theme_bw()
 
## figfn2 <- paste(outdir2, "Figure_gary.pdf", sep="")
## ggsave(figfn2, p, width=4, height=4)






###
### May put supp 

## ###
## ### scatter plots
## fn <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/heatmaps/zscore_corr/02.5.1.zscore_matrix_unionDEGs_pearson_-RNA-CTRL.SES_PCs_sex_age_and_treats_generem_.txt"
## corr_rna <- read.table(fn, header=T)
## psycho_var <- c("DSES_09", "chronic_sum", "SES", "cytocomp",
##   "PSS_all_mean", "Chol_HDL", "BPd_avg", "isel", "pr_comp")

## colnames(corr_rna) <- psycho_var
## rownames(corr_rna) <- psycho_var

## opfn <- paste(outdir2, "2.2_corr_rna.rds", sep="")
## write_rds(corr_rna, file=opfn)



## ###
## fn0 <- paste(outdir2, "2.2_corr_rna.rds", sep="")
## corr_rna <- read_rds(fn0)


## fn2 <- paste(outdir2, "2.2_corr_psycho.rds", sep="") 
## corr_motif <- read_rds(fn2)


## ## identical(sort(colnames(corr_rna)), sort(colnames(corr_motif)))

## var2 <- sort(colnames(corr_rna))
## nvar <- length(var2)
## var_comb <- NULL
## for (i in 1 : (nvar-1)){
##    for ( j in (i+1) : nvar){ 
##       ###
##       var_comb <-c(var_comb,  paste(var2[i], var2[j], sep="#"))
##    }    
## }


## ###
## ### plot data frame
## plot_data <- NULL
## for (ii in var_comb){
##     ##
##     v1 <- gsub("#.*", "", ii)
##     v2 <- gsub(".*#", "", ii)
##     ##
##     df0 <- data.frame(comb=ii, rr_x_motif = corr_motif[v1, v2], rr_y_rna = corr_rna[v1, v2]) 
##     ##
##     plot_data <- rbind(plot_data, df0)
## }
 
## opfn <- paste(outdir2, "Figure2.4_scatter.plotdata.rds", sep="")
## write_rds(plot_data, opfn)

## ###
## p0 <- ggplot(plot_data, aes(x = rr_x_motif, y = rr_y_rna))+
##    geom_point(size=3, shape=23, color="#F87660")+
##    scale_x_continuous("PCC from motif", limits=c(-1, 1))+
##    scale_y_continuous("PCC from RNA", limits=c(-1, 1))+
##    geom_abline(color="grey")+
##    geom_vline(xintercept=0, color="blue", linetype=2, linewidth=1)+
##    geom_hline(yintercept=0, color="blue", linetype=2, linewidth=1)+ 
##    theme_bw()+
##    theme(axis.title=element_text(size=10),
##          axis.text = element_text(size=10))

## ##
## figfn <- paste(outdir2, "Figure2.4_corr_compare.scatter.png", sep="")
## ggsave(figfn, p0, width=500, height=500, units="px", dpi=120)




###
### 













## ##########################
## ### QQ plots 
## ########################

## fn <- "./psycho_var_dir/psychosocial_top10.txt"
## var2 <- sort(read.table(fn, header=F)$V1)
## ## varSel <- sort(var2)[c(3, 4, 6:9)]



## outdir <- "./2_motif.outs/correct_excludeX/Cluster_res0.07_th20_psycho/"
## fn_ls <- list.files(outdir, ".results.txt.gz$")
## fn_ls <- fn_ls[grepl("CTRL", fn_ls)]

   
## res_all <- NULL
## for ( ii in fn_ls){
##     ###

##     fn0 <- paste(outdir, ii, sep="")
##     res0 <- fread(fn0, header=T, data.table=F)

##     for ( var0 in var2){

##        res <- res0%>%filter(psycho_variable==var0) 
##        ntest <- nrow(res)

##        cat(var0, res$nind[1], ntest, "\n")
        
##        res <- res%>%
##            arrange(pval_t)%>%
##            mutate(observed=-log10(pval_t), expected=-log10(ppoints(ntest)))
##        ###
##        ### 
##        res_all <- rbind(res_all, res)
##    }     
## }
 
## opfn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
## fwrite(res_all, file=opfn, sep="\t", quote=F, na=NA)


## ###
## ### plots

## col2 <- hue_pal()(9) ### default ggplot color

## infn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
## plotDF <- fread(infn, header=T, data.table=F)

## var2 <- sort(unique(plotDF$psycho_variable))
## varSel2 <- var2[c(3, 4, 6:9)]
## plotDF2 <- plotDF%>%filter(psycho_variable%in%varSel2)

## ## motifs <- plotDF2%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()

## ###
## p0 <- ggplot(plotDF2, aes(x=expected, y=observed, color=Cluster))+
##     rasterise(geom_point(size=0.06), dpi=300)+
##     geom_abline(color="grey")+
##     scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
##         "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C9"="#FF61C3"),
##         guide=guide_legend(override.aes=list(size=3)))+
##     facet_wrap(~psycho_variable, scales="free_y", ncol=3)+
##     xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
##     ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           legend.text=element_text(size=9),
##           legend.key.size=grid::unit(0.4, "cm"),
##           axis.title=element_text(size=12),
##           axis.text=element_text(size=10),
##           strip.text=element_text(size=12))

## ###
## ###
## figfn <- paste(outdir2, "Figure2.1_th20_psycho_ctrl.qq.png", sep="")
## ggsave(figfn, p0, width=850, height=500, units="px", dpi=120)
