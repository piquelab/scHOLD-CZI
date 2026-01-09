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
library(scales)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(ggrastr)
library(colorspace)

library(dendsort, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
###

rm(list=ls())

outdir <- "./888_pub.outs/supp/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


###
### for differential analysis
### 2.1_differentialMotifs.R, 2.1_submit_diff.sh and 2.1_summary_diff.R
### outdir, ./2_motif.outs/correct_excludeX/Cluster_res0.07_th20_psycho/ results from differentil analysis
### ./2_motif.outs/correct_excludeX/Summary_res0.07_psycho/ summary results 


###
### integration RNA analysis analysis
### ./3_regulatory.outs/correct_excludeX/1_enrich_DAM.fisher.txt, enrichment of DEG in DAMs (overview)
## ./3_regulatory.outs/list_cluster_var_nDEGs20_nTF1_new.txt
## ./3_regulatory.outs/correct_excludeX/enrich_motif/2_comb_enrich.motif.txt.gz
## ./2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz


###
### 

varSel <- c("ISEL_Mean", "PSS_all_mean")
cl_atac <- c(paste("C", c(0, 1, 2, 4, 6), sep="")) ##paste("C", 0:7, sep="")
 
indir <- "1_DiffPeak.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"

###
fn <- paste(indir, "2_psycho_plotData.comb.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
res <- res%>%filter(psycho_variable%in%varSel, Cluster%in%cl_atac)%>%
    mutate(comb=paste(Cluster, psycho_variable, sep="_"))


MCl_name <- c("C0"="A0_R0_T CD4+", "C1"="A1_R1_T CD8+", "C2"="A2_R2_NK" , "C3"="A3_R0_T CD4+", "C4"="A4_R3_Monocyte",
              "C5"="A5_R0_T CD4+", "C6"="A6_R4_B", "C7"="A7_R0_T CD4+")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")

col_cl <- c("A0_R0_T CD4+"="#FF7F00",
   "A3_R0_T CD4+"=darken("#FF7F00", 0.1), "A5_R0_T CD4+"=lighten("#FF7F00", 0.2),
   "A7_R0_T CD4+"=lighten("#FF7F00",0.5),
   "A1_R1_T CD8+"="#E6E600", "A2_R2_NK"="#4DAF4A", "A4_R3_Monocyte"="#984EA3", "A6_R4_B"="#D97986")

plotDF <- res%>%
    mutate(psycho_full=var_full[psycho_variable],
           MCl=MCl_name[Cluster],
           MCl_val=as.numeric(gsub(".*_R|_.*", "", MCl)),
           MCl2=fct_reorder(MCl, MCl_val))




##################
### Bar plots
####################

sigs <- res%>%group_by(Cluster, psycho_variable)%>%
    summarize(ny=sum(p.adjusted<0.1, na.rm=TRUE), .groups="drop")%>%ungroup()%>%filter(ny>0)


sigs <- sigs%>%
    mutate(MCl=MCl_name[Cluster], var2=var_full[psycho_variable],
           MCl_val=as.numeric(gsub(".*_R|_.*", "", MCl)),
           MCl2=fct_reorder(MCl, MCl_val))

anno_df <- sigs%>%filter(ny>0)
ymax <- max(sigs$ny)*1.2

p <- ggplot(sigs, aes(x=MCl2, y=ny))+
   geom_bar(aes(fill=MCl2), stat="identity")+
   scale_fill_manual(values=col_cl, guide="none")+
   geom_text(data=anno_df, aes(x=MCl2, y=ny, label=ny, vjust=-0.5), size=4.5)+
   facet_wrap(~var2, ncol=2)+
   scale_y_continuous(breaks=seq(0, ymax, by=100), limits=c(0, ymax))+ 
   theme_bw()+
   theme(strip.text=element_text(size=14),
         axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=12),
         axis.text.y=element_text(size=12))

###
### save pdf
figfn <- paste(outdir, "FigS5_DAR.barplots.pdf", sep="")
ggsave(filename=figfn, p, width=7, height=4.2)

### save png
figfn <- paste(outdir, "FigS5_DAR.barplots.png", sep="")
ggsave(figfn, p, width=700, height=420, units="px", dpi=120)
 

###################
### QQ plots 
#################



###
p0 <- ggplot(plotDF, aes(x=expected, y=observed, color=MCl2))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=col_cl,
        guide=guide_legend(override.aes=list(size=3)))+
    facet_wrap(~psycho_full, scales="free_y", ncol=2)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          strip.text=element_text(size=12))

###
###
figfn <- paste(outdir, "FigS5_DAR_psycho_ctrl.qq.png", sep="")
ggsave(figfn, p0, width=800, height=380, units="px", dpi=120)

figfn <- paste(outdir, "FigS5_DAR_psycho_ctrl.qq.pdf", sep="")
ggsave(figfn, p0, width=6, height=3)




########################
### Heatmap for DARs
###########################
## col2 <- hue_pal()(9) ### default ggplot color

varSel <- c("ISEL_Mean", "PSS_all_mean")

cl_atac <- paste("C", c(0, 1, 2, 4), sep="")
 
indir <- "1_DiffPeak.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"

###
fn <- paste(indir, "2_psycho_plotData.comb.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
res <- res%>%filter(psycho_variable%in%varSel, Cluster%in%cl_atac)%>%
    mutate(comb=paste(Cluster, psycho_variable, sep="_"))

DARs <- lapply(cl_atac, function(ii){
    ##
    res2 <- res%>%filter(Cluster==ii, !is.na(p.adjusted))
    unique(res2$gene)
})
DARs <- Reduce(intersect, DARs)  ## 70426

## DARs <- res%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique() ## 670 DARs
## res <- res%>%filter(Cluster%in%cl_atac)

res2 <- res%>%filter(gene%in%DARs)

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

## corr_mat <- cor(mat, method="pearson")

## ## opfn <- paste(outdir2, "2.2_corr_psycho.rds", sep="")
## ## write_rds(corr_mat, file=opfn)

###
### significance

nvar <- ncol(mat)
sig_mat <- matrix(0, nvar, nvar)
corr_mat <- matrix(0, nvar, nvar)
for ( i in 1 : nvar){
    ##
    for ( j in i : nvar){
        ##
        zi <- mat[, i]
        zj <- mat[, j]
        notna <- !is.na(zi)&!is.na(zj)
        zi <- zi[notna]
        zj <- zj[notna]
        test_rr <- cor.test(zi, zj, method = "pearson")
        pval <- test_rr$p.value
        is_sig <- ifelse(pval < 0.05, 1, 0)
        rr <- test_rr$estimate
        ##
        sig_mat[i, j] <- is_sig
        sig_mat[j, i] <- is_sig
        corr_mat[i, j] <- rr
        corr_mat[j, i] <- rr        
    }    
}

colnames(corr_mat) <- colnames(mat)
rownames(corr_mat) <- colnames(mat)
colnames(sig_mat) <- colnames(mat)
rownames(sig_mat) <- colnames(mat)


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
    ### show_heatmap_legend=F)
    ###, row_names_max_width=unit(6, "cm")
    ## column_names_rot = -45,
    heatmap_legend_param = list(title = "PCC",
        at = seq(-1, 1, by=0.5),
        title_gp = gpar(fontsize=12),
        labels_gp = gpar(fontsize=12),
        grid_width = grid::unit(0.38, "cm"),
        legend_height = grid::unit(6, "cm")))
   
    ## cell_fun=function(j, i, x, y, width, height, fill){
    ##     grid.text(round(corr3[i, j], digits=3), x, y, gp=gpar(fontsize=9))
    ##     })

## hlg <- Legend(col_fun=mycol, title="PCC", at=seq(-1, 1, by=0.5),
##     title_gp = gpar(fontsize=12),
##     labels_gp = gpar(fontsize=12),
##  ##   grid_height = grid::unit(0.38, "cm"),
##     legend_width = grid::unit(5, "cm"), direction="horizontal")
               
###
### save figures 
figfn <- paste(outdir, "FigS5b_DARs_psycho_corr.heatmap.pdf", sep="")
pdf(figfn, width = 6.8, height = 6)
p0 <- draw(p0) ##padding=unit(c(2, 2, 2, 10), "mm"))
dev.off()

### png
figfn <- paste(outdir, "FigS5b_DARs_psycho_corr.heatmap.png", sep="")
png(figfn, width = 700, height = 600, res=120)
p0 <- draw(p0) ##padding=unit(c(2, 2, 2, 10), "mm"))
dev.off()


### legend 
## figfn <- paste(outdir, "FigS5_DARs_hlegend.pdf", sep="")
## pdf(figfn, width=2.3, height=0.8)
## draw(hlg, x=unit(0.5, "npc"), y=unit(0.5, "npc"), just="bottomright")
## dev.off()

## ## png
## figfn <- paste(outdir, "FigS5_DARs_hlegend.png", sep="")
## png(figfn, width=250, height=80, res=120)
## draw(hlg, x=unit(0.5, "npc"), y=unit(0.5, "npc"), just="bottomright")
## dev.off()


##################################
### Heatmap for DEG results 
###################################

###
### variable
varSel <- c("cytocomp", "ISEL_Mean", "PSS_all_mean")
clSel <- paste("C", 0:4, sep="")

###
### cluster and variable full name and colors  
MCl_name <- c("C0"="R0_T CD4+", "C1"="R1_T CD8+", "C2"="R2_NK",  "C3"="R3_Monocyte", "C4"="R4_B")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress", "cytocomp"="Cytokines")
col_cl <- c("C0"="#FF7F00", "C1"="#E6E600", "C2"="#4DAF4A", "C3"="#984EA3", "C4"="#D97986")
col_var <- c("ISEL_Mean"="#3A84D8", "PSS_all_mean"="#E34234", "cytocomp"="#ff7f00")


###
indir <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"
fn <- paste(indir, "Ali_RNA_psycho.diff.DESeq.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
res <- res%>%mutate(statistic=logFC/SE)%>%
    dplyr::rename(gene=identifier, p.adjusted=padj, Cluster=cluster, psycho_variable=var) 

res <- res%>%filter(psycho_variable%in%varSel, Cluster%in%clSel)%>%
    mutate(comb=paste(Cluster, psycho_variable, sep="_"))
 
## DARs <- lapply(cl_atac, function(ii){
##     ##
##     res2 <- res%>%filter(Cluster==ii, !is.na(p.adjusted))
##     unique(res2$gene)
## })
## DARs <- Reduce(intersect, DARs)  ## 70426

geneSel <- res%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique() ## 2754
## res <- res%>%filter(Cluster%in%cl_atac)

res2 <- res%>%filter(gene%in%geneSel)


###
### data for heatmap
mat_df <- res2%>%pivot_wider(id_cols=gene, names_from=comb, values_from=statistic, values_fill=NA)
mat <- mat_df%>%column_to_rownames(var="gene") 

## corr_mat <- cor(mat, method="pearson")

## ## opfn <- paste(outdir2, "2.2_corr_psycho.rds", sep="")
## ## write_rds(corr_mat, file=opfn)

###
### significance

nvar <- ncol(mat)
sig_mat <- matrix(0, nvar, nvar)
corr_mat <- matrix(0, nvar, nvar)
for ( i in 1 : nvar){
    ##
    for ( j in i : nvar){
        ##
        zi <- mat[, i]
        zj <- mat[, j]
        notna <- !is.na(zi)&!is.na(zj)
        zi <- zi[notna]
        zj <- zj[notna]
        test_rr <- cor.test(zi, zj, method = "pearson")
        pval <- test_rr$p.value
        is_sig <- ifelse(pval < 0.05, 1, 0)
        rr <- test_rr$estimate
        ##
        sig_mat[i, j] <- is_sig
        sig_mat[j, i] <- is_sig
        corr_mat[i, j] <- rr
        corr_mat[j, i] <- rr        
    }    
}

colnames(corr_mat) <- colnames(mat)
rownames(corr_mat) <- colnames(mat)
colnames(sig_mat) <- colnames(mat)
rownames(sig_mat) <- colnames(mat)


###
corr2 <- corr_mat*sig_mat

df_col <- data.frame(col_name=colnames(corr2))%>%
    mutate(variable=gsub("C._", "", col_name), var2=var_full[variable],
           cluster=gsub("_.*", "", col_name),
           MCl=MCl_name[cluster], MCl_val=as.numeric(gsub(".*R|_.*", "", MCl)))%>%
    arrange(var2, MCl_val, .locale="en")

df_col <- df_col%>%mutate(comb_new=paste(var2, MCl, sep="_"))

###
col_sort <- df_col$col_name
corr2 <- corr2[col_sort, col_sort]



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
    ### show_heatmap_legend=F)
    ###, row_names_max_width=unit(6, "cm")
    ## column_names_rot = -45,
    heatmap_legend_param = list(title = "PCC",
        at = seq(-1, 1, by=0.5),
        title_gp = gpar(fontsize=12),
        labels_gp = gpar(fontsize=12),
        grid_width = grid::unit(0.38, "cm"),
        legend_height = grid::unit(6, "cm")))
   
    ## cell_fun=function(j, i, x, y, width, height, fill){
    ##     grid.text(round(corr3[i, j], digits=3), x, y, gp=gpar(fontsize=9))
    ##     })

## hlg <- Legend(col_fun=mycol, title="PCC", at=seq(-1, 1, by=0.5),
##     title_gp = gpar(fontsize=12),
##     labels_gp = gpar(fontsize=12),
##  ##   grid_height = grid::unit(0.38, "cm"),
##     legend_width = grid::unit(5, "cm"), direction="horizontal")
               
###
### save figures 
figfn <- paste(outdir, "FigS5_DEGs_psycho_corr.heatmap.pdf", sep="")
pdf(figfn, width = 6.8, height = 6)
p0 <- draw(p0) ##padding=unit(c(2, 2, 2, 10), "mm"))
dev.off()

### png
figfn <- paste(outdir, "FigS5_DEGs_psycho_corr.heatmap.png", sep="")
png(figfn, width = 680, height = 600, res=120)
p0 <- draw(p0) ##padding=unit(c(2, 2, 2, 10), "mm"))
dev.off()


############################
### stratified QQ plots ####
############################


cl_atac <- paste("C", c(0, 1, 2, 4, 6), sep="")

indir <- "./1_DiffPeak.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"

fn <- paste(indir, "3.2_summary_grDEG_short.tsv", sep="")
summ <- read_tsv(fn)

fn2 <- paste(indir, "3.2_plotData_grDEG.comb.txt.gz", sep="")
res <- fread(fn2, sep="\t", header=T, data.table=F)
##
res <- res%>%filter(comb%in%summ$comb)
res <- res%>%filter(Cluster%in%cl_atac)


###
MCl_name <- c("C0"="A0_R0_T CD4+", "C1"="A1_R1_T CD8+", "C2"="A2_R2_NK" , "C3"="A3_R0_T CD4+", "C4"="A4_R3_Monocyte",
              "C5"="A5_R0_T CD4+", "C6"="A6_R4_B", "C7"="A7_R0_T CD4+")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress", "cytocomp"="Cytokines")

### names
varSel <- sort(unique(res$psycho_variable))
var_val <- 1:length(varSel)
names(var_val) <- varSel


plotDF <- res%>%
    mutate(var_value=as.numeric(var_val[psycho_variable]),
           var_name=var_full[psycho_variable], 
           MCls = MCl_name[Cluster],
           comb2 = paste(MCls, var_name, sep="_"),
           comb_sort=fct_reorder(comb, var _value))
nlen <- length(unique(plotDF$comb))
p2 <- ggplot(plotDF, aes(x=expected_gr, y=observed_gr, color=factor(is_gr)))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=c("gr0"="#f1b6da", "gr1"="#d01c8b"),
        labels=c("gr0"="Not in DEG", "gr1"="In DEG"),
        guide=guide_legend(override.aes=list(size=2)))+
  facet_grid(var_name~MCls, scales="free")+  ## ncol=7       
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=10),
          legend.key.size=grid::unit(0.4, "cm"),
          ##legend.position="inside",
          ##legend.position.inside=c(0.95, 0.12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          strip.text=element_text(size=10))

### 
###
figfn <- paste(outdir, "FigS5_stratified_grDEG.qq.pdf", sep="")
ggsave(figfn, p2, width=12, height=6)

### png 
figfn2 <- paste(outdir, "FigS5_stratified_grDEG.qq.png", sep="")
ggsave(figfn2, p2, width=1200, height=600, units="px", dpi=120)





####################
### supp tables 
######################

fn <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
res <- fread(fn, header=T, sep="\t")
 
cl_sel <- paste("C", c(0, 1, 2, 4, 6), sep="")
varSel <- c("ISEL_Mean", "PSS_all_mean")

MCl_name <- c("C0"="A0_R0_T CD4+", "C1"="A1_R1_T CD8+", "C2"="A2_R2_NK" , "C3"="A3_R0_T CD4+", "C4"="A4_R3_Monocyte",
              "C5"="A5_R0_T CD4+", "C6"="A6_R4_B", "C7"="A7_R0_T CD4+")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")
 
##
res2 <- res%>%filter(psycho_variable%in%varSel, Cluster%in%cl_sel, psycho_variable%in%varSel)%>%
    mutate(MCls=MCl_name[Cluster], psycho_variable=var_full[psycho_variable])%>%
    dplyr::select(MCls, Cluster, psycho_variable, treat, motif_name=gene,
                  estimate, stderror, statistic, p.value=pval_t, p.adjusted=padj_t, nind)

opfn <- paste(outdir, "TableS04_DAM_results.txt.gz", sep="")
fwrite(res2, file=opfn, sep="\t", row.names=F, col.names=T, quote=F, na=NA)


## var0 <- varSel[1]
## res%>%filter(padj_t<0.1, psycho_variable == var0,
##      Cluster%in%c("C0", "C5", "C6", "C7"))%>%pull(gene)%>%unique()%>%length()

## res%>%filter(padj_t<0.1, psycho_variable == var0,
##      Cluster%in%c("C1"))%>%pull(gene)%>%unique()%>%length()

## res%>%filter(padj_t<0.1, psycho_variable == var0,
##      Cluster%in%c("C2"))%>%pull(gene)%>%unique()%>%length()

## res%>%filter(padj_t<0.1, psycho_variable == var0,
##      Cluster%in%c("C3"))%>%pull(gene)%>%unique()%>%length()

## res%>%filter(padj_t<0.1, psycho_variable == var0,
##      Cluster%in%c("C4"))%>%pull(gene)%>%unique()%>%length()




####################
### QQ plots
#####################




MCl_name <- c("C0"="A0_R0_T CD4+", "C1"="A1_R1_T CD8+", "C2"="A2_R2_NK" , "C3"="A3_R0_T CD4+", "C4"="A4_R3_Monocyte",
              "C5"="A5_R0_T CD4+", "C6"="A6_R4_B", "C7"="A7_R0_T CD4+")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")

col_cl <- c("A0_R0_T CD4+"="#FF7F00",
   "A3_R0_T CD4+"=darken("#FF7F00", 0.1), "A5_R0_T CD4+"=lighten("#FF7F00", 0.2),
   "A7_R0_T CD4+"=lighten("#FF7F00",0.5),
   "A1_R1_T CD8+"="#E6E600", "A2_R2_NK"="#4DAF4A", "A4_R3_Monocyte"="#984EA3", "A6_R4_B"="#D97986")

###
cl_sel <- paste("C", 0:7, sep="")
varSel <- c("ISEL_Mean", "PSS_all_mean")


###
### read data 
fn <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
res <- fread(fn, header=T, sep="\t")

res <- res%>%filter(Cluster%in%cl_sel, psycho_variable%in%varSel)


plotDF <- res%>%
    mutate(psycho_full=var_full[psycho_variable],
           MCl=MCl_name[Cluster],
           MCl_val=as.numeric(gsub(".*_R|_.*", "", MCl)),
           MCl2=fct_reorder(MCl, MCl_val))

###
p0 <- ggplot(plotDF, aes(x=expected, y=observed, color=MCl2))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=col_cl,
        guide=guide_legend(override.aes=list(size=3)))+
    facet_wrap(~psycho_full, scales="free_y", ncol=2)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          strip.text=element_text(size=12))

###
###
figfn <- paste(outdir, "FigS6_DAM_psycho_ctrl.qq.png", sep="")
ggsave(figfn, p0, width=800, height=400, units="px", dpi=120)

figfn <- paste(outdir, "FigS6_DAM_psycho_ctrl.qq.pdf", sep="")
ggsave(figfn, p0, width=6, height=3)



##################
#### Barplots 
##################

indir <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"

varSel <- c("ISEL_Mean", "PSS_all_mean")

fn <- paste(indir, "2_psycho_plotData.comb.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t")
## res <- res%>%filter(Cluster!="C9", psycho_variable%in%varSel)
cl_atac <- paste("C", 0:7, sep="")
res <- res%>%filter(Cluster%in%cl_atac, psycho_variable%in%varSel)
 

sigs <- res%>%group_by(Cluster, psycho_variable)%>%
    summarize(ny=sum(padj_t<0.1), .groups="drop")%>%ungroup()%>%filter(ny>0)

MCl_name <- c("C0"="A0_R0_T CD4+", "C1"="A1_R1_T CD8+", "C2"="A2_R2_NK" , "C3"="A3_R0_T CD4+", "C4"="A4_R3_Monocyte",
              "C5"="A5_R0_T CD4+", "C6"="A6_R4_B", "C7"="A7_R0_T CD4+")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Psychological Stress")
## col_cl <- c("C0"="#E60A0A",  "C1"="#A65192", "C2"="#0570B0", "C3"="#00AF72", "C4"="#87CEFA",
##             "C5"="#F9A602", "C6"="#F9A602", "C7"="#F9A602")

## celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")
## col_cl <- c("C0"="#E69F00","C1"="#984EA3", "C2"="#D55E00", "C3"="#CC79A7", "C4"="#56B4E9",
##             "C5"="#F9A602", "C6"="#F9A602", "C7"="#F9A602")

col_cl <- c("C0"="#FF7F00",
   "C3"=darken("#FF7F00", 0.1), "C5"=lighten("#FF7F00", 0.2), "C7"=lighten("#FF7F00",0.5),
   "C1"="#E6E600", "C2"="#4DAF4A", "C4"="#984EA3", "C6"="#D97986")

sigs <- sigs%>%
    mutate(MCl=MCl_name[Cluster], var2=var_full[psycho_variable],
           MCl_val=as.numeric(gsub(".*_R|_.*", "", MCl)),
           MCl2=fct_reorder(MCl, MCl_val))

anno_df <- sigs%>%filter(ny>0)
ymax <- max(sigs$ny)*1.2

p <- ggplot(sigs, aes(x=MCl2, y=ny))+
   geom_bar(aes(fill=factor(Cluster)), stat="identity")+
   scale_fill_manual(values=col_cl, guide="none")+
   geom_text(data=anno_df, aes(x=MCl2, y=ny, label=ny, vjust=-0.5), size=4.5)+
   facet_wrap(~var2, ncol=2)+
   scale_y_continuous(breaks=seq(0, ymax, by=10), limits=c(0, ymax))+ 
   theme_bw()+
   theme(strip.text=element_text(size=12),
         axis.title=element_blank(),
         axis.text.x=element_text(angle=60, hjust=1, vjust=1, size=10),
         axis.text.y=element_text(size=10))

###
### save
figfn <- paste(outdir, "FigS6_DAM.barplots.pdf", sep="")
ggsave(filename=figfn, p, width=7, height=4.2)

### png
figfn2 <- paste(outdir, "FigS6_DAM.barplots.png", sep="")
ggsave(filename=figfn2, p, width=700, height=420, units="px", dpi=120)




###
### END
