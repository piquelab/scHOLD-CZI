###
library(Matrix)
library(tidyverse)
library(data.table)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(Signac)
###library(SeuratWrappers)
library(SeuratObject)
library(ArchR) ##, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
##
library(cowplot)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(openxlsx)

## library(gganimate, lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## install.packages("gganimate", lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
## install.packages("gifski", lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")  
 
rm(list=ls())
 

outdir2 <- "./1_results_all/Batch_correct_Harmony/"
if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

###
dir_ArchR <- "./ArchR_all_output/"
##dir_ArchR <- "./ArchR_all_output/"
proj <- loadArchRProject(path=dir_ArchR)





######################################################
### summary metrics before fastdemuxlet filtering
#######################################################


meta <- read_rds("./1_results_all/Batch_correct_Harmony/0_all.meta.rds")
meta <- meta%>%mutate(NEW_BARCODE=rownames(meta), EXP=gsub("#.*","", NEW_BARCODE))

df0 <-  meta%>%dplyr::group_by(EXP)%>%
    dplyr::summarise(ncell=n(), nFrag_per_cell=mean(nFrags),
         TSS_enrich=mean(TSSEnrichment), Nucleo_ratio=mean(NucleosomeRatio),
         Black_ratio=mean(BlacklistRatio), .groups="drop")%>%dplyr::ungroup()

df0 <- df0%>%dplyr::mutate(Batch_val=as.numeric(gsub("^HOLD|-ATAC-.*", "", EXP)),
           treat=gsub(".*-ATAC-", "", EXP))

###
#### save
opfn <- paste(outdir2, "0.0_all_summary_stat.tsv", sep="")
write_tsv(df0, file=opfn)



fn <- "./1_results_all/Batch_correct_Harmony/0.0_all_summary_stat.tsv"
x <- read_tsv(fn)


## fn <- "./1_results_all/Batch_correct_Harmony/0.0_summary_stat.tsv"
## x <- read_tsv(fn)


###
### 
nnSel <- c("EXP", "Batch_val", "ncell", "nFrag_per_cell", "TSS_enrich", "Nucleo_ratio", "Black_ratio")
plotDF <- x%>%dplyr::select(all_of(nnSel))%>%
    mutate(Batch=gsub("-ATAC-.*", "", EXP),
           Batch_sort=fct_reorder(Batch, Batch_val),
           EXP_sort=fct_reorder(EXP, Batch_val),
           nFrag2 = round(nFrag_per_cell))

###
### #cells per EXP
p1 <- ggplot(plotDF, aes(x=EXP_sort, y=ncell, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("#Barcodes")+
   ggtitle("CellRanger-ATAC: Barcodes per experiment")+ 
   geom_text(aes(label=ncell),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure0.1_qc_barcodes.pdf", sep="")
ggsave(figfn, p1, width=10, height=5)



###
### Fragments per cell for each EXP
p2 <- ggplot(plotDF, aes(x=EXP_sort, y=nFrag2, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("#Fragments")+
   ggtitle("CellRanger-ATAC: Fragments per cell(mean)")+ 
   geom_text(aes(label=nFrag2),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure0.2_qc_Fragments_per_cell.pdf", sep="")
ggsave(figfn, p2, width=10, height=5)



###
### TSS_enrich per cell for each EXP
p2 <- ggplot(plotDF, aes(x=EXP_sort, y=TSS_enrich, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("TSS enrichment score")+
   ggtitle("CellRanger-ATAC: TSS enrichment score per cell(mean)")+ 
   geom_text(aes(label=nFrag2),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure0.3_qc_Fragments_per_cell.pdf", sep="")
ggsave(figfn, p2, width=10, height=5)



###############################################################################
### summary QC statistics for each library after fastdemuxlet filtering 
##############################################################################

outdir2 <- "./1_results_all/Batch_correct_Harmony/"

meta <- read_rds("./1_results_all/Batch_correct_Harmony/1_filtered.meta.rds")
meta <- meta%>%mutate(NEW_BARCODE=rownames(meta), EXP=gsub("#.*","", NEW_BARCODE))

df0 <-  meta%>%dplyr::group_by(EXP)%>%
    dplyr::summarise(ncell=n(), nind=length(unique(sampleID)), nFrag_per_cell=mean(nFrags),
         TSS_enrich=mean(TSSEnrichment), Nucleo_ratio=mean(NucleosomeRatio),
         Black_ratio=mean(BlacklistRatio), .groups="drop")%>%dplyr::ungroup()
 
df2 <- meta%>%dplyr::group_by(EXP, sampleID)%>%
    dplyr::summarise(ncell_ind=n(), nFrag_per_cell=mean(nFrags), .groups="drop")%>%dplyr::ungroup()

###
df2 <- df2%>%dplyr::group_by(EXP)%>%
    dplyr::summarise(ncell_ind_md=round(median(ncell_ind)), 
              ncell_ind_min=round(min(ncell_ind)),
              ncell_ind_max=round(max(ncell_ind)),
              nFrags_per_cell_md=round(median(nFrag_per_cell)), 
              nFrags_per_cell_min=round(min(nFrag_per_cell)),
              nFrags_per_cell_max=round(max(nFrag_per_cell)), .groups="drop")%>%dplyr::ungroup()

summ_df <- df0%>%dplyr::left_join(df2, by="EXP")%>%
    dplyr::mutate(Batch_val=as.numeric(gsub("^HOLD|-ATAC-.*", "", EXP)),
           treat=gsub(".*-ATAC-", "", EXP))

####
###
opfn <- paste(outdir2, "1.0_filterd_summary_stat.tsv", sep="")
write_tsv(summ_df, file=opfn)



####
#### QC plot (1) barplots of basic statistics for each library
fn <- "./1_results_all/Batch_correct_Harmony/1.0_filterd_summary_stat.tsv"
x <- read_tsv(fn)
## x2 <- x%>%arrange(Batch_val)%>%
##     dplyr::select(EXP, nind, ncell, ncell_ind_md, nFrag_per_cell, TSS_enrich, Nucleo_ratio, Black_ratio)
## opfn <- "./888_pub_figures_output/supp/TableS00_summ_statistics.tsv"
## write_tsv(x2, file=opfn)



###
### 
nnSel <- c("EXP", "Batch_val", "ncell", "ncell_ind_md",
           "nFrag_per_cell", "nFrags_per_cell_md", "TSS_enrich", "Nucleo_ratio", "Black_ratio")
plotDF <- x%>%dplyr::select(all_of(nnSel))%>%
    mutate(Batch=gsub("-ATAC-.*", "", EXP),
           Batch_sort=fct_reorder(Batch, Batch_val),
           EXP_sort=fct_reorder(EXP, Batch_val),
           nFrag2 = round(nFrag_per_cell))
           ## nFrag2_md = round(nFrags_per_cell_md)

###
### #cells per EXP
p1 <- ggplot(plotDF, aes(x=EXP_sort, y=ncell, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("#Barcodes")+
   ggtitle("CellRanger-ATAC: Barcodes per experiment")+ 
   geom_text(aes(label=ncell),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure1.1_qc_barcodes.pdf", sep="")
ggsave(figfn, p1, width=10, height=5)


###
### 
p1b <- ggplot(plotDF, aes(x=EXP_sort, y=ncell_ind_md, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("#Barcodes per individual")+
   ggtitle("CellRanger-ATAC: Barcodes/individual (median, experiment)")+ 
   geom_text(aes(label=ncell_ind_md),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure1.1_qc_per_indi_barcodes.pdf", sep="")
ggsave(figfn, p1b, width=10, height=5)


###
### Fragments per cell for each EXP
p2 <- ggplot(plotDF, aes(x=EXP_sort, y=nFrag2, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("#Fragments")+
   ggtitle("CellRanger-ATAC: Fragments per cell(mean)")+ 
   geom_text(aes(label=nFrag2),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure1.2_qc_Fragments_per_cell.pdf", sep="")
ggsave(figfn, p2, width=10, height=5)


###
### TSS enrichment score 
plotDF <- plotDF%>%
    mutate(TSS_enrich=round(TSS_enrich, 2),
          Nucleo_ratio=round(Nucleo_ratio, 2), Black_ratio=round(Black_ratio, 2))

p3 <- ggplot(plotDF, aes(x=EXP_sort, y=TSS_enrich, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("TSS enrichment score")+
   ggtitle("CellRanger-ATAC: TSS enrichment score")+ 
   geom_text(aes(label=TSS_enrich),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure1.3_qc_TSS_enrich.pdf", sep="")
ggsave(figfn, p3, width=10, height=5)


###
### 
p4 <- ggplot(plotDF, aes(x=EXP_sort, y=Nucleo_ratio, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("Nucleosome signal")+
   ggtitle("CellRanger-ATAC: nucleosome signal")+ 
   geom_text(aes(label=Nucleo_ratio),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure1.4_qc_Nucleo_ratio.pdf", sep="")
ggsave(figfn, p4, width=10, height=5)


###
###
p5 <- ggplot(plotDF, aes(x=EXP_sort, y=Black_ratio, fill=factor(Batch_sort)))+
   geom_bar(stat="identity")+
   xlab("")+
   ylab("Blacklist ratio")+
   ggtitle("CellRanger-ATAC: blacklist ratio")+ 
   geom_text(aes(label=Black_ratio),vjust=-0.5, size=2.5)+
   theme_bw()+
   theme(legend.title=element_blank(),
         legend.key.size=grid::unit(0.5, "cm"),
         legend.text=element_text(size=9),
         plot.title=element_text(hjust=0.5, size=12),
         axis.text.x=element_text(angle=90, hjust=1, size=9))

figfn <- paste(outdir2, "Figure1.5_qc_Black_ratio.pdf", sep="")
ggsave(figfn, p5, width=10, height=5)
 



### 
### QC plot (2) Distribution of basic statistics for total libraries 
meta <- read_rds("./1_results_all/Batch_correct_Harmony/1_filtered.meta.rds")

###
feature_nn <- c("TSSEnrichment", "NucleosomeRatio", "BlacklistRatio")                
figs_ls <- lapply(feature_nn, function(ii){
    ##
    df0 <- data.frame(cells=rownames(meta), yy=meta%>%pull(ii))
    p0 <- ggplot(df0, aes(x="", y=yy))+
        geom_violin(fill="grey", color="grey50", linewidth=0.6)+
        geom_boxplot(width=0.2, color="grey50", outlier.color="red", outlier.size=0.6, outlier.stroke=0.35)+
        ggtitle(ii)+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.title=element_text(size=9, hjust=0.5))
    ###
    p0
 })

pcomb <- plot_grid(plotlist=figs_ls, nrow=1, ncol=3) 

###
### save
figfn <- paste(outdir2, "Figure1.6_vlnplot.pdf", sep="")
pdf(figfn, width=8, height=3)
print(pcomb)
dev.off()
###
 
####
### plot (3) scatter plots


outdir2 <- "./1_results_all/Batch_correct_Harmony/"

### ArchR project
dir_ArchR <- "./ArchR_all_output/"
addArchRThreads(12)
proj <- loadArchRProject(path=dir_ArchR)


###
### 
## meta <- read_rds("./1_results_all/Batch_correct_Harmony/1_filtered.meta.rds")
## dd <- meta%>%dplyr::select(TSSEnrichment, nFrags, EXP)%>%
##     mutate(log10_nFrags=log10(nFrags))

## xmax <- quantile(dd$log10_nFrags, probs=0.99)
## ymax <- quantile(dd$TSSEnrichment, probs=0.99)

## p1 <- ggplot(dd, aes(x=log10_nFrags, y=TSSEnrichment))+
##     geom_point(color)+
##     stat_density_2d(aes(fill=after_stat(level)), geom="polygon", contour=T)+
##     geom_vline(xintercept=3, linetype="dashed", color="grey")+
##     geom_hline(yintercept=4, linetype="dashed", color="grey")+
##     scale_fill_viridis_c(option="plasma")+
##     scale_x_continuous(bquote(log[10]~("Unique fragment")), limits=c(log10(500), xmax))+
##     scale_y_continuous("TSS Enrichment", limits=c(0, ))+
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           legend.text=element_text(size=9),
##           legend.key.size=grid::unit(0.7, "cm"))

## figfn <- paste(outdir2, "Figure1.6_qc_TSS_scatter_density.pdf", sep="")
## ggsave(figfn, p1)
    

plotDF <- getCellColData(proj, select=c("log10(nFrags)", "TSSEnrichment")) 
p <- ggPoint(x = plotDF[,1], y = plotDF[,2],
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = c(log10(500), quantile(plotDF[,1], probs = 1)),
        ylim = c(0, quantile(plotDF[,2], probs = 1))
    )+
    geom_hline(yintercept = 4, linetype = "dashed") +
    geom_vline(xintercept = 3, linetype="dashed")

## plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

figfn <- paste(outdir2, "Figure1.6_qc_TSS_scatter_density.pdf", sep="")
pdf(figfn)
p
dev.off()

###
### Fragmentsize distribution                      
p2 <- plotFragmentSizes(proj, returnDF=TRUE)

pp2 <- ggplot(p2, aes(x=fragmentSize, y=fragmentPercent))+
   geom_line(color="brown2", linewidth=0.8)+
   xlab("Size of Fragments(bp)")+
   ylab("Fragments(%)")+     
   theme_bw()+
   theme(axis.title=element_text(size=12),
         axis.text=element_text(size=12)) 
   
figfn <- paste(outdir2, "Figure1.6_qc_FragSize_distribution.pdf", sep="")
pdf(figfn)
pp2
dev.off()
    
###
opfn <- paste(outdir2, "Figure1.6_qc_FragSize_distribution.plotdata.rds", sep="")
write_rds(p2, file=opfn)



###
### QC for filtering data

outdir2 <- "./1_results_all/Batch_correct_Harmony/"

### ArchR project
dir_ArchR <- "./ArchR_all_output/"
proj <- loadArchRProject(path=dir_ArchR)
 
x <- as.data.frame(getCellColData(proj))


### 
colSel <- c("TSSEnrichment", "PromoterRatio", "nFrags")
summ <- map_dfr(sort(unique(x$Cluster2)), function(ii){
   ##
   x2 <- x%>%filter(Cluster2==ii)%>%dplyr::select(all_of(colSel))     
   summ0 <- apply(x2, 2, quantile, c(0, 0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1))
   summ0 <- as.data.frame(summ0)%>%rownames_to_column(var="th")
   summ0$Cluster <- ii
   summ0
})

opfn <- paste(outdir2, "2.0_qc_cluster.tsv", sep="")
write_tsv(summ, file=opfn)


###
### wide table

colSel2 <- c("th", paste("C", 0:10, sep=""))
qc_name <- c("TSSEnrichment", "PromoterRatio", "nFrags")

### TSS
for (qc0 in qc_name){
    
   summ2 <- summ%>%dplyr::select(all_of(c("th", "Cluster", qc0)))
   names(summ2)[3] <- "qc_index" 
   summ22 <- summ2%>%pivot_wider(id_cols=th, names_from=Cluster, values_from=qc_index)%>%
      as.data.frame()
   summ22 <- summ22[, colSel2]

   ###   
   opfn2 <- paste(outdir2, "2.0_qc_", qc0, "_cluster.xlsx", sep="")
   write.xlsx(summ22, file=opfn2)
}

### nFrags



###########################################################################
### Use Heatmap to compare cluster results between current and previous 
############################################################################

outdir2 <- "./1_results_all/Batch_correct_Harmony/"

###
### no DEX results
fn <- paste(outdir2, "1.0_cluster.rds", sep="")
df <- read_rds(fn)
df2 <- df%>%mutate(Barcode=rownames(df))%>%
    dplyr::select(Barcode, cluster=ATAC_cluster_res.0.12,
                  umap_1="Harmony#UMAP_Dimension_1", umap_2="Harmony#UMAP_Dimension_2")

###
summ <- df2%>%group_by(cluster)%>%summarize(ncell=n(),.groups="drop")%>%ungroup()%>%
    arrange(dplyr::desc(ncell))
## opfn <- paste(outdir2, "Table1.0_res0.07.summ.tsv", sep="")
## write_tsv(summ%>%mutate(cluster_new=(1:nrow(summ))-1, prop=round(ncell/sum(ncell), 3)), file=opfn)
          


cl_id <- paste("C", (1:nrow(summ))-1, sep="")
names(cl_id) <- summ$cluster

df2 <- df2%>%mutate(cluster_new=cl_id[cluster])%>%
    dplyr::select(all_of(c("Barcode", "cluster_new", "umap_1", "umap_2")))



###
### umap plots
### p1 
plotDF <- df2%>%
    mutate(cl_val=as.numeric(gsub("C", "", cluster_new)),
           cluster2 = fct_reorder(cluster_new, cl_val))
 
p1 <- ggplot(plotDF, aes(x=umap_1, y=umap_2, colour=as.factor(cluster2)))+
    rasterise(geom_point(size=0.1), dpi=300)+
       ## facet_wrap(~factor(EXP_sort), ncol=8)+
    guides(col=guide_legend(override.aes=list(size=1.2), ncol=1))+
    ggtitle("umap (new)")+
    theme_bw()+
    theme(###strip.text=element_text(size=9),
          legend.title=element_blank(),
          plot.title=element_text(hjust=0.5, size=8),
          ###legend.position.inside=c(0.8, 0.3),
          legend.background=element_blank(),
          legend.box.background=element_blank(),
          legend.key.size=grid::unit(0.4,"cm"),
          legend.text=element_text(size=9),
          axis.title=element_text(size=9),
          axis.text=element_text(size=9))




###
### w/ DEX results
fn0 <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/old/analyses_new_2025_08_07/1.2_ArchR_process/1_results_all/Batch_correct_Harmony/1.0_cluster.rds"
old <- read_rds(fn0)
old2 <- old%>%mutate(Barcode=rownames(old))%>%
    dplyr::select(Barcode, cluster=ATAC_cluster_res.0.07,
                  umap_1="Harmony#UMAP_Dimension_1", umap_2="Harmony#UMAP_Dimension_2")

###
summ <- old2%>%group_by(cluster)%>%summarize(ncell=n(),.groups="drop")%>%ungroup()%>%
    arrange(dplyr::desc(ncell))
cl_id <- paste("C", (1:nrow(summ))-1, sep="")
names(cl_id) <- summ$cluster

old2 <- old2%>%mutate(cluster_old=cl_id[cluster])%>%
    dplyr::select(all_of(c("Barcode", "cluster_old", "umap_1", "umap_2")))



###
### umap plots 
### p2
plotDF2 <- old2%>%filter(Barcode%in%df2$Barcode)%>%
    mutate(cl_val=as.numeric(gsub("C", "", cluster_old)),
           cluster2 = fct_reorder(cluster_old, cl_val))

p2 <- ggplot(plotDF2, aes(x=umap_1, y=umap_2, colour=as.factor(cluster2)))+
    rasterise(geom_point(size=0.1), dpi=300)+
       ## facet_wrap(~factor(EXP_sort), ncol=8)+
    guides(col=guide_legend(override.aes=list(size=1.2), ncol=1))+
    ggtitle("umap (old)")+
    theme_bw()+
    theme(###strip.text=element_text(size=9),
          legend.title=element_blank(),
          plot.title=element_text(hjust=0.5, size=8),
          ###legend.position.inside=c(0.8, 0.3),
          legend.background=element_blank(),
          legend.box.background=element_blank(),
          legend.key.size=grid::unit(0.4,"cm"),
          legend.text=element_text(size=9),
          axis.title=element_text(size=9),
          axis.text=element_text(size=9))




###
### heatmap plots

barSel <- intersect(rownames(df2), rownames(old2))

dfcomb<- df2[, !grepl("umap", colnames(df2))]%>%
    left_join(old2[, !grepl("umap", colnames(old2))], by="Barcode")
              

df_summ <- dfcomb%>%group_by(cluster_new, cluster_old)%>%summarize(ncell=n(), .groups="drop")%>%
    group_by(cluster_new)%>%
    mutate(nt = sum(ncell), Perc = ncell/nt)%>%
    ungroup()
df_summ <- df_summ%>%
    mutate(cluster_new_cl=as.numeric(gsub("C", "", cluster_new)),
           cluster_new_sort=fct_reorder(cluster_new, cluster_new_cl))

p3 <- ggplot(df_summ, aes(x=cluster_new_sort, y=cluster_old, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
       low="#ffffcc", high="#e31a1c", limits=c(0,1), n.breaks=5, na.value=NA,
       guide=guide_legend(keywidth=unit(0.35, "cm"), keyheight=unit(0.8, "cm")))+
   xlab("cluster 0.12 (new)")+ylab("cluster 0.07 (old)")+
   ggtitle("Fraction of cells")+
   theme_bw()+
   theme(legend.title=element_blank(),
         plot.title=element_text(hjust=0.5, size=9),
         axis.text.x=element_text(hjust=0.5, size=9),
         axis.text.y=element_text(size=9),
         axis.title=element_text(size=9))


###
pcomb <- plot_grid(p1, p2, p3, nrow=1, ncol=3, align="h", axis="tb", rel_widths=c(0.95, 0.95, 1))

figfn <- paste(outdir2, "Figure2.00_cluster_newVSold.pdf", sep="")
ggsave(figfn, pcomb, width=12, height=3.5)





################################
### umap show batch effect
################################


outdir2 <- "./1_results_all/Batch_correct_Harmony/"

###
### no DEX results
fn <- paste(outdir2, "1.0_cluster.rds", sep="")
df_cl <- read_rds(fn)%>%pull(ATAC_cluster_res.0.07)


## plot data 
x <- as.data.frame(getCellColData(proj))

identical(rownames(x), rownames(df_cl))

x$Clusters <- df_cl$ATAC_cluster_res.0.07
x <- x%>%dplyr::select(Clusters, EXP, treat)

##
umap <- getEmbedding(proj, embedding="UMAP")
##
identical(rownames(x), rownames(umap))
df0 <- cbind(umap, x)

##
names(df0)[1:2] <- c("UMAP_1", "UMAP_2")



### number of clusters 
summ_cl <- as.data.frame(table(df0$Clusters))
names(summ_cl) <- c("Cluster", "ncell")
summ_cl <- summ_cl%>%dplyr::arrange(dplyr::desc(ncell))
new_cl <- paste("C", 0:(nrow(summ_cl)-1), sep="")
names(new_cl) <- summ_cl$Cluster 
 


###
### sort experiment
plotDF <- df0%>%
    mutate(EXP2=gsub("-ATAC-", "-", EXP),
           EXP_val=as.numeric(gsub("^HOLD|-ATAC-.*", "", EXP)),
           EXP_sort=fct_reorder(EXP2, EXP_val),
           new_cluster=new_cl[Clusters], 
           cl_val=as.numeric(gsub("^C", "", new_cluster)),
           Clusters_sort=fct_reorder(new_cluster, cl_val))

pp <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2, colour=as.factor(cl_val)))+
   rasterise(geom_point(size=0.1), dpi=300)+
   facet_wrap(~factor(EXP_sort), ncol=7)+
   guides(col=guide_legend(override.aes=list(size=1.2), ncol=1))+
   theme_bw()+
   theme(strip.text=element_text(size=9),
         legend.title=element_blank(),
         ###legend.position.inside=c(0.8, 0.3),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.4,"cm"),
         legend.text=element_text(size=9),
         axis.title=element_text(size=9),
         axis.text=element_text(size=9))

figfn <- paste(outdir2, "Figure3.1_EXP.umap.pdf", sep="")
ggsave(figfn, pp, width=11, height=6.2)


###
### facet By treatment
pp2 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2, colour=as.factor(cl_val)))+
   rasterise(geom_point(size=0.1), dpi=300)+
   facet_wrap(~factor(treat), ncol=2)+
   guides(col=guide_legend(override.aes=list(size=1.2), ncol=1))+
   theme_bw()+
   theme(strip.text=element_text(size=9),
         legend.title=element_blank(),
         ###legend.position.inside=c(0.8, 0.3),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.4,"cm"),
         legend.text=element_text(size=9),
         axis.title=element_text(size=9),
         axis.text=element_text(size=9))

figfn <- paste(outdir2, "Figure3.2_treat.umap.pdf", sep="")
ggsave(figfn, pp2, width=6, height=3)




###
## p0 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2, colour=as.factor(cl_val)))+
##    rasterise(geom_point(size=0.1), dpi=300)+
##    ## facet_wrap(~factor(EXP_sort), ncol=8)+
##    guides(col=guide_legend(override.aes=list(size=1.2), ncol=1))+
##    theme_bw()+
##    theme(strip.text=element_text(size=9),
##          legend.title=element_blank(),
##          ###legend.position.inside=c(0.8, 0.3),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=9),
##          axis.title=element_text(size=9),
##          axis.text=element_text(size=9))
 
## figfn <- paste(outdir2, "Figure1.0.umap.png", sep="")
## ggsave(figfn, p0, width=520, height=420, units="px", dpi=120)

## ###
## ### save
## summ <- as.data.frame(table(as.factor(plotDF$cl_val)))
## names(summ) <- c("cluster", "ncell")
## ##
## opfn <- paste(outdir2, "1_summ_cl.tsv", sep="")
## write_tsv(summ, file=opfn)



######################################################
### colored by Library, batch and treatment  
#######################################################


###
### input data 
## x <- as.data.frame(getCellColData(proj))
## x <- x%>%mutate(treat2=ifelse(grepl("DEX", treat), "LPS-DEX", treat))
## umap <- getEmbedding(proj, embedding="UMAP2")

## identical(rownames(x), rownames(umap))

## ###
## x2 <- x%>%dplyr::select(Clusters, sampleID, EXP, treat=treat2)
## df2 <- cbind(umap, x2)
## names(df2)[1:2] <- c("UMAP_1", "UMAP_2")

## ###
## df2 <- df2%>%
##     mutate(EXP2=gsub("-ATAC-", "-", EXP), Batch=gsub("-ATAC-.*", "", EXP))



                       

 
## ###
## ### visulization-option 1

## plotDF <- df2%>%mutate(val=as.numeric(gsub("^HOLD", "", Batch)),
##                        EXP2_sort=fct_reorder(EXP2, val),
##                        Batch_sort=fct_reorder(Batch, val))

 
## p0_1 <- ggplot(plotDF%>%arrange(EXP2), aes(x=UMAP_1, y=UMAP_2, color=EXP2_sort, alpha=factor(treat)))+
##    rasterise(geom_point(size=0.4), dpi=300)+
##    scale_alpha_manual(values=c("CTRL"=0.2, "LPS"=0.5, "LPS-DEX"=0.5), guide="none")+ 
##    guides(col=guide_legend(override.aes=list(size=1.5),ncol=2))+
##    theme_bw()+
##    theme(
##          legend.title=element_blank(),
##          ###legend.position=c(0.85, 0.25),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=10),
##          axis.text=element_text(size=10),
##          axis.title=element_text(size=10))

## figfn <- paste(outdir2, "Figure1.3_umap.lib.png", sep="")
## ggsave(figfn, p0_1, width=680, height=380, units="px", dpi=100)
 


## ###
## ### color by Batch 
## p0_2 <- ggplot(plotDF%>%arrange(Batch), aes(x=UMAP_1, y=UMAP_2, color=Batch_sort, alpha=factor(treat)))+
##    rasterise(geom_point(size=0.4), dpi=300)+
##    scale_alpha_manual(values=c("CTRL"=0.2, "LPS"=0.5, "LPS-DEX"=0.5), guide="none")+ 
##    guides(col=guide_legend(override.aes=list(size=1.5), ncol=1))+
##    theme_bw()+
##    theme(
##          legend.title=element_blank(),
##          ###legend.position=c(0.85, 0.25),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=10),
##          axis.text=element_text(size=10),
##          axis.title=element_text(size=10))

## figfn <- paste(outdir2, "Figure1.3_umap.batch.png", sep="")
## ggsave(figfn, p0_2, width=450, height=380, units="px", dpi=100)



## ###
## ### color by treatment
## ## plotDF%>%arrange(treat)

## p0_3 <- ggplot(plotDF%>%arrange(treat), aes(x=UMAP_1, y=UMAP_2, color=treat, alpha=factor(treat)))+
##    rasterise(geom_point(size=0.4), dpi=300)+
##    scale_alpha_manual(values=c("CTRL"=1, "LPS"=0.2, "LPS-DEX"=0.05), guide="none")+ 
##    guides(col=guide_legend(override.aes=list(size=1.5),ncol=1))+
##    theme_bw()+
##    theme(
##          legend.title=element_blank(),
##          ###legend.position=c(0.85, 0.25),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=10),
##          axis.text=element_text(size=10),
##          axis.title=element_text(size=10))

## figfn <- paste(outdir2, "Figure1.3_umap.treat.png", sep="")
## ggsave(figfn, p0_3, width=430, height=380, units="px", dpi=100)


## ###
## ### save
## pcomb <- plot_grid(p0_1, p0_2, p0_3, ncol=3, rel_widths=c(1.6, 1.05, 1))
 
## ## png(figfn, width=1100, height=300, res=120)
## ## pcomb
## ## dev.off()
## figfn2 <- paste(outdir2, "Figure1.3_umap.zzz.png", sep="")
## ggsave(figfn2, pcomb, width=1400, height=300, units="px", dpi=100) 






##############################################################################################
### visulization-option 2, colored by two classes, in (library, batch or treatment) or not 
##############################################################################################

###
### plot data for library 
## EXPs <- sort(sort(unique(df2$EXP2)))
## EXP_df <- data.frame(EXP_name=EXPs)%>%
##     mutate(EXP_val=as.numeric(gsub("HOLD|-.*", "", EXP_name)))%>%arrange(EXP_val)
## EXP_sort <- EXP_df$EXP_name

## for ( i in 1:4){

## i0 <- (i-1)*8+1
## i1 <- ifelse(i==4, 31, (i-1)*8+8)
## subi <- i0:i1    
## plotDF <- map_dfr(EXP_sort[subi], function(ii){
##    ##
##    df0 <- df2%>%mutate(gr_col=ifelse(EXP2==ii, "HOLD", "bg"))%>%arrange(gr_col) ##dplyr::desc(gr_col))
##    df0$gr_facet <- ii
##    df0
## })    

## cat("Lib", i, "\n")    

## ###
## ###
## plotDF <- plotDF%>%
##     mutate(gr_value=as.numeric(gsub("^HOLD|-.*", "", gr_facet)),
##            gr2_facet=fct_reorder(gr_facet, gr_value))  


## col1 <- c("HOLD"="#F8766D", "bg"="grey")
## ## "LPS+DEX"="#e31a1c",
## ## "PHA"="#a6cee3", "PHA+DEX"="#1f78b4", "bg"="#c7e9c0")
## ### "bg"="#bdbdbd") ##, 
## ##
## p0 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2))+
##    rasterise(geom_point(aes(colour=factor(gr_col), alpha=gr_col), size=0.4), dpi=300)+
##    facet_wrap(~factor(gr2_facet), ncol=4, nrow=2)+
##    scale_colour_manual(values=col1,
##       breaks=c("bg", "HOLD"),
##       labels=c("bg"="Not in EXP", "HOLD"="In EXP"),
##       guide=guide_legend(override.aes=list(size=1.5)))+
##    scale_alpha_manual(values=c("HOLD"=0.3, "bg"=0.3), guide="none")+
##    ##guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
##    theme_bw()+
##    theme(strip.text=element_text(size=12),
##          legend.title=element_blank(),
##          ###legend.position=c(0.85, 0.25),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=10),
##          axis.text=element_text(size=10),
##          axis.title=element_text(size=10))
 
## figfn <- paste(outdir2, "Figure1.2_umap.lib", i, ".png", sep="")
## ggsave(figfn, p0, width=880, height=410, units="px", dpi=100)

## }    
    


 
## ### colored by Batch 
## ########################
 
## batchs <- sort(sort(unique(df2$Batch)))
## plotDF <- map_dfr(batchs, function(ii){
##    ##
##    df0 <- df2%>%mutate(gr_col=ifelse(Batch==ii, "HOLD", "bg"))%>%arrange(gr_col) ##dplyr::desc(gr_col))
##    df0$gr_facet <- ii
##    df0
## })    


## ###
## ###
## plotDF <- plotDF%>%
##     mutate(gr_value=as.numeric(gsub("^HOLD", "", gr_facet)),
##            gr2_facet=fct_reorder(gr_facet, gr_value))  


## col1 <- c("HOLD"="#F8766D", "bg"="#c7e9c0")
## ## "LPS+DEX"="#e31a1c",
## ## "PHA"="#a6cee3", "PHA+DEX"="#1f78b4", "bg"="#c7e9c0")
## ### "bg"="#bdbdbd") ##, 
## ##
## p1 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2))+
##    rasterise(geom_point(aes(colour=factor(gr_col), alpha=gr_col), size=0.4), dpi=300)+
##    facet_wrap(~factor(gr2_facet), ncol=4, nrow=1)+
##    scale_colour_manual(values=col1,
##       breaks=c("bg", "HOLD"),
##       labels=c("bg"="Not in batch", "HOLD"="In batch"),
##       guide=guide_legend(override.aes=list(size=1.5)))+
##    scale_alpha_manual(values=c("HOLD"=1, "bg"=0.4),
##                       guide="none")+
##    ##guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
##    theme_bw()+
##    theme(strip.text=element_text(size=12),
##          legend.title=element_blank(),
##          ###legend.position=c(0.85, 0.25),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=10),
##          axis.text=element_text(size=10),
##          axis.title=element_text(size=10))

## figfn <- paste(outdir2, "Figure1.2_umap.batch.png", sep="")
## ggsave(figfn, p1, width=850, height=250, units="px", dpi=100)




## ###
## ### colored by treatment
## treats <- sort(sort(unique(df2$treat)))
## plotDF <- map_dfr(treats, function(ii){
##    ##
##    df0 <- df2%>%mutate(gr_col=ifelse(treat==ii, treat, "bg"))%>%arrange(gr_col) ##dplyr::desc(gr_col))
##    df0$gr_facet <- ii
##    df0
## })    


## plotDF <- plotDF%>%mutate(gr2_facet=gr_facet)


## ###
## ###
## ## col1 <- c("CTRL"="#828282", "LPS"="#fb9a99", "bg"="#c7e9c0")
## col1 <- c("bg"="grey", "CTRL"="#F8766D", "LPS"="#00BA38", "LPS-DEX"="#619CFF")
## ## "LPS+DEX"="#e31a1c",
## ## "PHA"="#a6cee3", "PHA+DEX"="#1f78b4", "bg"="#c7e9c0")
## ### "bg"="#bdbdbd") ##, 
## ##
 
## p2 <- ggplot(plotDF%>%dplyr::filter(Batch=="HOLD1"), aes(x=UMAP_1, y=UMAP_2))+
##    rasterise(geom_point(aes(colour=factor(gr_col), alpha=gr_col), size=0.4), dpi=300)+
##    facet_wrap(~factor(gr2_facet), ncol=3, nrow=1, dir="h")+
##    scale_colour_manual(values=col1,
##       breaks=c("bg", "CTRL", "LPS", "LPS-DEX"),
##       labels=c("bg"="Not in", "CTRL"="CTRL", "LPS"="LPS", "LPS-DEX"="LPS-DEX"),
##       guide=guide_legend(override.aes=list(size=1.5)))+
##    scale_alpha_manual(values=c("CTRL"=0.1, "LPS"=0.1, "LPS-DEX"=0.1, "bg"=0.8),
##                       guide="none")+
##    ## transition_manual(gr_col)+    
##    ##guides(col=guide_legend(override.aes=list(size=2),ncol=3))+
##    theme_bw()+
##    theme(strip.text=element_text(size=12),
##          legend.title=element_blank(),
##          ###legend.position=c(0.85, 0.25),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=10),
##          axis.text=element_text(size=10),
##          axis.title=element_text(size=10))


## figfn <- paste(outdir2, "Figure1.2_umap.treat.png", sep="")
## ggsave(figfn, p2, width=680, height=250, units="px", dpi=100)


## anim_save(filename=figfn, p2)
###
### animate plots


## mm <- assay(mat)



##############################
### color by clusters
##############################


## plot data 
## x <- as.data.frame(getCellColData(proj))
## x <- x%>%mutate(treat2=ifelse(grepl("DEX", treat), "LPS-DEX", treat))
 
option <- "Batch_correct_Harmony"
outdir2 <- paste("./1_results_all/", option, "/", sep="")

### input data
infn <- paste(outdir2, "1.0_cluster.rds", sep="")
x <- read_rds(infn)
names(x)[8:9] <- c("UMAP_1", "UMAP_2")

## ###
## ### input data
## infn_meta <- paste(outdir2, "1.0_cluster.rds", sep="")
## x <- read_rds(infn_meta)
## ##
## umap <- getEmbedding(proj, embedding="UMAP")
## names(umap) <- c("UMAP_1", "UMAP_2")
 

## cat(identical(rownames(x), rownames(umap)), "\n")
## if ( ! identical(rownames(x), rownames(umap)) ) umap <- umap[rownames(x),]

###
## plot list for each resolution # 
res_atac <- c("0.05", "0.06", "0.07", "0.08", "0.1", "0.12", "0.15")
figs_ls <- lapply(res_atac, function(ii){
    
### clusters     
rnSel <- paste("ATAC_cluster_res.", ii, sep="")
atac_cl <- x%>%pull(rnSel)

    
### number of clusters 
summ <- as.data.frame(table(atac_cl))
names(summ) <- c("Cluster", "ncell")
summ <- summ%>%dplyr::arrange(dplyr::desc(ncell))
###     
new_cl <- as.character((1:nrow(summ))-1)
names(new_cl) <- summ$Cluster

### re-name cluster
atac_cl <- new_cl[atac_cl]    

### plot data    
plotDF <- x%>%dplyr::select(all_of(c("UMAP_1", "UMAP_2")))    
plotDF$Cluster <- atac_cl    

plotDF <- plotDF%>%
    mutate(cl_val=as.numeric(Cluster), Cluster2=paste("C", Cluster, sep=""),
           Cluster2_sort=fct_reorder(Cluster2, cl_val))
    
cat("ATAC resolution", ii, "\n")    

###
p0 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2, colour=Cluster2_sort))+
   rasterise(geom_point(size=0.1), dpi=300)+
   ## facet_wrap(~factor(EXP_sort), ncol=8)+
   guides(col=guide_legend(title=paste("#", ii, sep=""),
          override.aes=list(size=1.2), ncol=1))+
   theme_bw()+
   theme(strip.text=element_text(size=9),
         legend.title=element_text(size=8),
         ###legend.position.inside=c(0.8, 0.3),
         legend.background=element_blank(),
         legend.box.background=element_blank(),
         legend.key.size=grid::unit(0.4,"cm"),
         legend.text=element_text(size=9),
         axis.title=element_text(size=9),
         axis.text=element_text(size=9))
p0
})
names(figs_ls) <- paste("atac", res_atac, sep="")

pcomb <- plot_grid(plotlist=figs_ls, ncol=4)

###
### save for all options 
figfn <- paste(outdir2, "Figure2_comb.umap.pdf", sep="")
ggsave(figfn, pcomb, width=11, height=5)


###
### save for 0.05, 0.07, 0.1, 0.12
pcomb2 <- plot_grid(figs_ls[[1]], figs_ls[[3]], figs_ls[[5]], figs_ls[[6]], figs_ls[[7]], ncol=1)
figfn <- paste(outdir2, "Figure2.0_comb.umap.v1.pdf", sep="")
ggsave(figfn, pcomb2, width=3.5, height=12)


###
### save for 0.12
ii <- "atac0.12"
figfn <- paste(outdir2, "Figure2.1_", ii, "_umap.pdf", sep="")
ggsave(figfn, figs_ls[[ii]], width=4.2, height=4)


### 
## figfn <- paste(outdir2, "Figure2.0_res0.07.umap.png", sep="")
## ggsave(figfn, p0, width=420, height=380, units="px", dpi=100)


##########################################
### colored by cluster for options
##########################################


## outdir2 <- "./1_results_all/option_estLSI_ncell380K/"
## if ( !file.exists(outdir2) ) dir.create(outdir2, showWarnings=F, recursive=T)

## dir_ArchR <- "./ArchR_all_output/"
## proj <- loadArchRProject(path=dir_ArchR)

## ## outdir2 <- "./1_results_all/Batch_correct_Harmony/"

## ### combine x and umap
## ## umap <- getEmbedding(proj, embedding="UMAP3")
## ## names(umap) <- c("UMAP_1", "UMAP_2")

## ## x <- as.data.frame(getCellColData(proj))
## ## x <- x%>%dplyr::select(EXP, treat, Cluster=Cluster_res0.12_sampleCell350K)
## ## ##, select="Cluster_res0.12"))
## ## identical(rownames(umap), rownames(x))
## ## x2 <- cbind(x, umap)

## ## ### save
## ## opfn <- paste(outdir2, "1_cluster.df.rds", sep="")
## ## write_rds(x2, file=opfn)



## ###
## ###
## fn <- paste(outdir2, "1.0_cluster.rds", sep="")
## x2 <- read_rds(fn)
## ## x2 <- x2%>%
## ##     dplyr::select(Cluster_res0.12=ATAC_cluster_res.0.12,
## ##     UMAP_1="Harmony#UMAP_Dimension_1", UMAP_2="Harmony#UMAP_Dimension_2")

## atac_cl <- x2$Cluster 
    
## ### number of clusters 
## summ <- as.data.frame(table(atac_cl))
## names(summ) <- c("Cluster", "ncell")
## summ <- summ%>%dplyr::arrange(dplyr::desc(ncell))
## ###     
## new_cl <- paste("C", as.character((1:nrow(summ))-1), sep="")
## names(new_cl) <- summ$Cluster

## ### re-name cluster
## atac_cl <- new_cl[atac_cl]    


## ### plot data    
## plotDF <- x2%>%dplyr::select(all_of(c("UMAP_1", "UMAP_2")))    
## plotDF$Cluster <- atac_cl    

## plotDF <- plotDF%>%
##     mutate(cl_val=as.numeric(as.character(gsub("C", "", Cluster))), Cluster2=fct_reorder(Cluster, cl_val))
    

## ###
## p0 <- ggplot(plotDF, aes(x=UMAP_1, y=UMAP_2, colour=Cluster2))+
##    rasterise(geom_point(size=0.1), dpi=300)+
##    ## facet_wrap(~factor(EXP_sort), ncol=8)+
##    guides(col=guide_legend(title="res0.12",
##           override.aes=list(size=1.2), ncol=1))+
##    theme_bw()+
##    theme(strip.text=element_text(size=9),
##          legend.title=element_text(size=8),
##          ###legend.position.inside=c(0.8, 0.3),
##          legend.background=element_blank(),
##          legend.box.background=element_blank(),
##          legend.key.size=grid::unit(0.4,"cm"),
##          legend.text=element_text(size=9),
##          axis.title=element_text(size=9),
##          axis.text=element_text(size=9))


## figfn <- paste(outdir2, "Figure1.0_umap.pdf", sep="")
## ggsave(figfn, p0, width=4.5, height=4)



## ###
## ###
## outdir2 <- "./1_results_all/Batch_correct_Harmony/"

## summ <- as.data.frame(table(proj$Cluster2))
## names(summ) <- c("Cluster", "ncell")

## summ <- summ%>%arrange(dplyr::desc(ncell))

## opfn <- paste(outdir2, "1_summary.tsv", sep="")
## write_tsv(summ, file=opfn)
