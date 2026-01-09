##
###
library(Matrix)
library(tidyverse)
library(data.table)
## library(Seurat) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
## library(SeuratDisk) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
## library(SeuratData) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
## library(Signac) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
## library(SeuratWrappers) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
## library(SeuratObject) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library/")
## library(ArchR, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(scales)
library(openxlsx)


rm(list=ls())

##option <- "option_estLSI_ncell380K"
option <- "harmony_default"
outdir <- paste("./2_Integrate_output/", option, "/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)






##############################################################
### ATAC with 0.12 resolution and RNA with 0.1 resolution
##############################################################

 
## infn_atac <- "./1_results_all/Batch_correct_Harmony/1.0_cluster.rds"
## df_atac <- read_rds(infn_atac)

## ###
## ### color setting 
## col2 <- hue_pal()(11)
## col_rna <- c(col2[c(1:3, 5, 7, 11)],"grey")
## names(col_rna) <- paste("C", 0:6, sep="")

## infn_meta <- paste(outdir, "1_RNA.res0.1.meta.rds", sep="")
## df0 <- as.data.frame(read_rds(infn_meta))

## identical(rownames(df_atac), rownames(df0))
## if ( ! identical(rownames(df_atac), rownames(df0)) ) df_atac <- df_atac[rownames(df0), ]

## df0 <- cbind(df0, df_atac%>%dplyr::select(Clusters=ATAC_cluster_res.0.12))

## ###
## ### new cluster 
## df_cl <- as.data.frame(table(df0$Clusters))
## names(df_cl) <- c("Cluster", "ncell")
## df_cl <- df_cl%>%dplyr::arrange(dplyr::desc(ncell))
## new_cl <- paste("C", 0:(nrow(df_cl)-1), sep="")
## names(new_cl) <- df_cl$Cluster 
 

## ## ### cell type value
## ## df_MCl <- as.data.frame(table(df0$predicted.id))
## ## names(df_MCl) <- c("Cluster", "ncell")
## ## df_MCl <- df_MCl%>%dplyr::arrange(dplyr::desc(ncell))

## ## MCl_val <- 0:(nrow(df_MCl)-1)
## ## names(MCl_val) <- df_MCl$Cluster



## plotDF <- df0%>%
##     mutate(Clusters_new=new_cl[Clusters],
##            cl_val = as.numeric(gsub("^C", "", Clusters_new)),
##            Clusters_new_sort = fct_reorder(Clusters_new, cl_val),
##            predicted.id = paste("C", predicted.id, sep=""),
##            predicted.id_val = as.numeric(gsub("C", "", predicted.id)),
##            predicted.id_sort = fct_reorder(predicted.id, predicted.id_val))
 
## ## ,   ##as.numeric(gsub("^C", "", Cluster)),
## ##            MCl_value=as.numeric(MCl_val[as.character(predicted.id)]),
## ##            predicted_sort=fct_reorder(predicted.id, MCl_value))

## ###  
## ### p1, atac-umap, colored by atac cluster id 
## p1 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2, colour=as.factor(Clusters_new_sort)))+
##    rasterise(geom_point(size=0.1), dpi=300)+
##    ## facet_wrap(~factor(EXP_sort), ncol=8)+
##    guides(col=guide_legend(override.aes=list(size=1.2), ncol=1))+
##    ggtitle("ATAC-umap colored by ATAC Cluster)")+ 
##    theme_bw()+
##    theme(###strip.text=element_text(size=9),
##          legend.title=element_blank(),
##          plot.title=element_text(hjust=0.5, size=8),
##          ###legend.position.inside=c(0.8, 0.3),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=9),
##          axis.title=element_text(size=9),
##          axis.text=element_text(size=9))

## ###
## ### p2, atac-umap, colored by transferred RNA cluster
## ## col_rna <- c(col2[c(1:5, 9)],"grey")
## ## names(col_rna) <- as.character(0:6)
## p2 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2, colour=as.factor(predicted.id_sort)))+
##    rasterise(geom_point(size=0.1), dpi=300)+
##    ## geom_point(size=0.1)  
##    ## facet_wrap(~factor(EXP_sort), ncol=8)+
##    scale_color_manual(values=col_rna,
##        guide=guide_legend(override.aes=list(size=1.2), ncol=1))+
##    ggtitle("ATAC-umap colored by RNA cluster")+ 
##    theme_bw()+
##    theme(###strip.text=element_text(size=9),
##          legend.title=element_blank(),
##          plot.title=element_text(hjust=0.5, size=8),
##          ###legend.position.inside=c(0.8, 0.3),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=9),
##          axis.title=element_text(size=9),
##          axis.text=element_text(size=9))



## ### 3 heatmap

## df2 <- plotDF%>%dplyr::select(predicted.id, Clusters=Clusters_new)

## ## %>%
## ##     mutate(seurat_clusters=factor(seurat_clusters))

## ###
 
## df_summ <- df2%>%
##    group_by(predicted.id, Clusters)%>%summarize(Freq=n(),.groups="drop")%>%ungroup()%>%
##    group_by(Clusters)%>%
##    mutate(nt=sum(Freq), Perc=Freq/nt)%>%ungroup()

## df_summ <- df_summ%>%mutate(Cluster_sort = fct_reorder(Clusters, desc(nt)))

## p3 <- ggplot(df_summ, aes(x=Cluster_sort, y=predicted.id, fill=Perc))+
##    geom_tile()+
##    scale_fill_gradient("Fraction of cells",
##       low="#ffffcc", high="#e31a1c", limits=c(0,1), n.breaks=5, na.value=NA,
##       guide=guide_legend(keywidth=unit(0.35, "cm"), keyheight=unit(0.8, "cm")))+
##    xlab("Cluster(ATAC)")+ylab("Cluster(RNA)")+
##    ggtitle("Fraction of cells")+ 
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          plot.title=element_text(hjust=0.5, size=9),
##          axis.text.x=element_text(hjust=0.5, size=9),
##          axis.text.y=element_text(size=9),
##          axis.title=element_text(size=9))

## ###
## ### p4, RNA umap, 0.1 resolution, colored by RNA cluster at 0.1 resolution  

## ### input
## infn <- "./2_Integrate_output/scHOLD_RNA_all.outs/sc_rna.meta.rds"
## x <- read_rds(infn)

## ii <- 0.1
## ii2 <- paste("RNA_snn_res.", ii, sep="")
## nnSel <- c(ii2, "umap_1", "umap_2") 
## plotDF_rna <- x%>%dplyr::select(all_of(nnSel))
## names(plotDF_rna) <- c("Cluster", "UMAP_1", "UMAP_2") 


## plotDF_rna <- plotDF_rna%>%mutate(cl_val = as.numeric(Cluster),
##        Cluster2 = paste("C", Cluster, sep=""),                           
##        Cluster2_sort = fct_reorder(Cluster2, cl_val))
## ###
## ### 
## p4 <- ggplot(plotDF_rna, aes(x=UMAP_1, y=UMAP_2, colour=Cluster2_sort))+
##    rasterise(geom_point(size=0.1), dpi=300)+
##    ## geom_point(size=0.1) 
##    ## facet_wrap(~factor(EXP_sort), ncol=8)+
##    scale_color_manual(values = col_rna,
##           guide=guide_legend(title=paste("#", ii, sep=""),
##           override.aes=list(size=1.2), ncol=1))+
##    ggtitle("RNA-umap")+ 
##    theme_bw()+
##    theme(###strip.text=element_text(size=9),
##          plot.title = element_text(hjust=0.5, size=9), 
##          legend.title=element_text(size=8),
##          ###legend.position.inside=c(0.8, 0.3),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=9),
##          axis.title=element_text(size=9),
##          axis.text=element_text(size=9))

## ###
## ## figfn <- paste(outdir, "Figure1_res0.15.umap2.png", sep="")
## ## ggsave(figfn, p0, width=500, height=420, units="px", dpi=120)

## pcomb <- plot_grid(p1, p2, p3, p4, nrow=2, ncol=2, align="hv", axis="tb")
## ## , rel_widths=c(0.95, 0.95, 1.05)) ## 1.25
   



## ###
## ### save figures
## ##
## ## figfn <- paste(outdir, "Figure1.0_atac0.07_rna0.15.comb.annot.png", sep="")
## ## ggsave(figfn, pcomb, width=1200, height=350, units="px", dpi=100)

## figfn <- paste(outdir, "Figure1.0_atac0.12_rna0.1.comb.annot.pdf", sep="")
## pdf(figfn, width = 8.5, height=7)
## pcomb
## dev.off()




## ###########################
## ### save data for final 
## ############################

## ###
## ### df_summ from the above heatmap plot
## summ <- df_summ%>%group_by(Clusters)%>%slice_max(order_by=Freq, n=1)%>%ungroup()%>%
##     arrange(Clusters)%>%
##     dplyr::select(ATAC_cluster=Clusters, RNA_cluster=predicted.id, Freq, nt, Perc)

## summ2 <- summ%>%
##     mutate(ATAC_cluster=gsub("C", "A", ATAC_cluster), RNA_cluster=gsub("C", "R",  RNA_cluster))

## MCls_name <- c("R0"="T CD4+", "R1"="T CD8+", "R2"="NK", "R3"="Monocyte", "R4"="B", "R5"="DC")

## ###
## ###
## summ2 <- summ2%>%
##     mutate(MCls=MCls_name[RNA_cluster], MCl_comb=paste(ATAC_cluster, RNA_cluster, MCls, sep="_"),
##            Perc = round(Perc, 3))
## summ2 <- summ2%>%arrange(desc(nt))

## ###
## ###
## opfn <- paste(outdir, "TableS_summ.tsv", sep="")
## write_tsv(summ2, file=opfn)








####################################################################################################
### (2) Visulization- how well  the clusters match cell type annotation with different resolution 
#####################################################################################################


## option <- "option_estLSI_ncell380K"
option <- "harmony_default"
outdir <- paste("./2_Integrate_output/", option, "/", sep="")

opt2 <- "Batch_correct_Harmony" 
infn_atac <- "./1_results_all/", option, "/1.0_cluster.rds", sep="")
df_atac <- read_rds(infn_atac)


###
### loop for atac resolution
res_ATAC <- c("0.05", "0.06", "0.07", "0.08", "0.1", "0.12", "0.15")
for (kk in res_ATAC){
###    
figfn <- paste(outdir, "Figure2_ATAC", kk, ".pred.heatmap.pdf", sep="")    

kk2 <- paste("ATAC_cluster_res.", kk, sep="")
atac_cl <- df_atac%>%pull(kk2)
names(atac_cl) <- rownames(df_atac)    

### new clusters    
summ <- as.data.frame(table(atac_cl))%>%arrange(desc(Freq))    
new_cl <- as.character((1:nrow(summ))-1)
names(new_cl) <- summ$atac_cl

atac_cl <- new_cl[atac_cl]
names(atac_cl) <- rownames(df_atac)

    
###
### plot data
res_rna <- as.character(c("0.1", "0.15", "0.2", "0.3", "0.4"))
plotDF <- map_dfr(res_rna, function(ii){
    ##
    fn <- paste(outdir, "1_RNA.res", ii, ".meta.rds", sep="")
    df0 <- read_rds(fn)%>%dplyr::select(predicted.id)
    cellSel <- rownames(df0)
    ## if not same, resort by cellSel
    cat("resolution", ii, identical(names(atac_cl), cellSel), "\n")
    if ( ! identical(names(atac_cl), cellSel)) atac_cl <- atac_cl[cellSel]
    df0$ATAC_cluster <- atac_cl    
    ### 
    df2  <- df0%>%
       group_by(predicted.id, ATAC_cluster)%>%
       summarize(Freq=n(),.groups="drop")%>%
       group_by(ATAC_cluster)%>%
       mutate(nt=sum(Freq), Perc=Freq/sum(Freq))%>%ungroup()
    df2$RNA_res <- paste("res", ii, sep="")
    df2
})


### sort
plotDF <- plotDF%>%
    mutate(rn_cl=paste(RNA_res, predicted.id, sep="_"),
           RNA_res_val=as.numeric(gsub("^res", "", RNA_res)),
           predicted.id_val=as.numeric(predicted.id))

plotDF2 <- plotDF%>%arrange(RNA_res_val, predicted.id_val)


###
###

mat <- plotDF2%>%
    pivot_wider(id_cols=rn_cl, names_from=ATAC_cluster, values_from=Perc, values_fill=NA)%>%
    column_to_rownames(var="rn_cl")

mat <- as.matrix(mat)

mat2 <- t(mat)

## sort
rnSel <- data.frame(rn=rownames(mat2))%>%mutate(rn_val=as.numeric(rn))%>%arrange(rn_val)%>%pull(rn)

mat2 <- mat2[rnSel,]    
    
###
###
x <- str_split(colnames(mat2), "_", simplify=T)

###
ha <- HeatmapAnnotation(rna_gr=anno_block(gp=gpar(fill=1:5),
     labels=c("0.1", "0.15", "0.2",  "0.3", "0.4"),
     labels_gp=gpar(col="white", fontsize=10)))

split <- gsub("res", "", x[,1])     

colnames(mat2) <- x[,2]
     
y2 <- as.numeric(mat2)

### set color
### "#ffffcc", high="#e31a1c"
## mycol <- colorRamp2(seq(0, 1, by=0.1), colorRampPalette(c("#ffffcc", "#e31a1c"))(11)) 
mycol <- colorRamp2(seq(0, 1, length.out=20), colorRampPalette(brewer.pal(n=7, name="YlGnBu"))(20))

###
### figures
  
p0 <- Heatmap(mat2, col=mycol, cluster_rows=F, cluster_columns=F, na_col="white", 
    rect_gp=gpar(col="grey", lwd=0.7),          
    show_row_names=T, row_names_gp=gpar(fontsize=10), 
    show_column_names=T, column_names_gp=gpar(fontsize=9),
    column_split=split,
    column_title="RNA cluster", column_title_gp=gpar(fontsize=10),
    row_title=paste("ATAC cluster(", kk, ")", sep=""), row_title_gp=gpar(fontsize=10), 
    top_annotation=ha, 
    heatmap_legend_param=list(title="Percent", title_gp=gpar(fontsize=9),
        at=seq(0, 1, by=0.25), labels_gp=gpar(fontsize=8),
        grid_width=grid::unit(0.4, "cm"), legend_height=grid::unit(5, "cm")),
    use_raster=T, raster_device="png")

###
### save figures
## png(figfn, width=950, height=320, res=100)
## p0 <- draw(p0)
## dev.off()
pdf(figfn, width=9.5, height=3.2)
p0 <- draw(p0)
dev.off()        

}




#############################
### calculate metrices
#############################



outdir <- "./2_Integrate_output/harmony_default/"
infn_atac <- "./1_results_all/Batch_correct_Harmony/1.0_cluster.rds"
df_atac <- read_rds(infn_atac)

res_ATAC <- c("0.05", "0.06", "0.07", "0.1", "0.12", "0.15")


prop_df <- NULL
for (kk in res_ATAC){

###    
kk2 <- paste("ATAC_cluster_res.", kk, sep="")
atac_cl <- df_atac%>%pull(kk2)
names(atac_cl) <- rownames(df_atac)    

cat("ATAC-res", kk, "\n")
    
res_rna <- as.character(c("0.1", "0.15", "0.2",  "0.3", "0.4"))
prop_df0 <- map_dfr(res_rna, function(ii){
    ##
    fn <- paste(outdir, "1_RNA.res", ii, ".meta.rds", sep="")
    df0 <- read_rds(fn)%>%dplyr::select(predicted.id, predicted.score)
    cellSel <- rownames(df0)
    ## if not same, resort by cellSel
    ## cat("resolution", ii, identical(names(atac_cl), cellSel), "\n")
    if ( ! identical(names(atac_cl), cellSel)) atac_cl <- atac_cl[cellSel]
    df0$ATAC_cluster <- atac_cl    
    ### 
    df2  <- df0%>%
       group_by(predicted.id, ATAC_cluster)%>%
       summarize(Freq=n(), summ_score=sum(predicted.score), .groups="drop")%>%
       group_by(ATAC_cluster)%>%
       mutate(nt=sum(Freq), perc=Freq/sum(Freq))%>%ungroup()

   df3 <- df2%>%group_by(ATAC_cluster)%>%
       slice_max(order_by=summ_score, n=1)%>%ungroup()%>%dplyr::arrange(desc(nt))%>%
        mutate(weight=nt/sum(nt), perc_w=perc*weight)

   summ_df <- data.frame(res_atac=kk, ncluster_atac=length(unique(df0$ATAC_cluster)),
                         res_rna=ii, ncluster_rna=length(unique(df0$predicted.id)),
                         summ_prop=round(sum(df3$perc_w), 3),
                         rangeL=round(min(df3$perc_w), 3),
                         rangeR=round(max(df3$perc_w), 3)) 
   summ_df  
})


###
###
prop_df <- rbind(prop_df, prop_df0)

}

###
### 
opfn <- paste(outdir, "table_summary.metric.xlsx", sep="")
write.xlsx(prop_df, file=opfn)





###################
### UMAP ATAC
##################
## outdir2 <- "./2_Integrate_output/harmony_default/"

## infn_atac <- "./1_results_all/Batch_correct_Harmony/1.0_cluster.rds"
## df_atac <- read_rds(infn_atac)


###
### UMAP RNA

## infn <- "./2_Integrate_output/scHOLD_RNA_all.outs/sc_rna.meta.rds"
## df_rna <- read_rds(infn)
 


            
