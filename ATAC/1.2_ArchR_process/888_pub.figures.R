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
###library(ArchR) ###, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
##
library(colorspace)
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

outdir <- "./888_pub_figures_output/main/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


###
### plot data

infn_atac <- "./1_results_all/Batch_correct_Harmony/1.0_cluster.rds"
df <- read_rds(infn_atac)

df2 <- df%>%mutate(Barcode=rownames(df))%>%
    dplyr::select(Barcode, cluster=ATAC_cluster_res.0.12,
                  umap_1="Harmony#UMAP_Dimension_1", umap_2="Harmony#UMAP_Dimension_2")

###
summ <- df2%>%dplyr::group_by(cluster)%>%summarize(ncell=n(),.groups="drop")%>%ungroup()%>%
    arrange(dplyr::desc(ncell))

## opfn <- paste(outdir2, "Table1.0_res0.07.summ.tsv", sep="")
## write_tsv(summ%>%mutate(cluster_new=(1:nrow(summ))-1, prop=round(ncell/sum(ncell), 3)), file=opfn)
          


cl_id <- paste("C", (1:nrow(summ))-1, sep="")
names(cl_id) <- summ$cluster

df2 <- df2%>%mutate(cluster_new=cl_id[cluster])%>%
    dplyr::select(all_of(c("Barcode", "cluster_new", "umap_1", "umap_2")))

opfn <- paste(outdir, "1_meta_final.rds", sep="")
write_rds(df2, file=opfn)

###
### add EXP and treatment
fn <- "./888_pub_figures_output/main/1_meta_final.rds"
x <- read_rds(fn)
x2 <- x%>%dplyr::select(NEW_BARCODE=Barcode)%>%
    mutate(EXP=gsub("#.*", "", NEW_BARCODE), treat=gsub(".*-", "", EXP))

opfn <- "./888_pub_figures_output/main/1_meta_final.rds"
write_rds(x2, file=opfn)


#########################
### atac umap
##########################

fn <- paste(outdir, "1_meta_final.rds", sep="")
df0 <- read_rds(file=fn)


## col_cl <- c("C0"="#E60A0A",  "C1"="#A65192", "C2"="#0570B0", "C3"="#00AF72", "C4"="#87CEFA",
##             "C5"="#F9A602", "C6"="#F9A602", "C7"="#F9A602")

## celltype <- c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B")
## col_cl <- c("C0"="#E69F00","C1"="#984EA3", "C2"="#D55E00", "C3"="#CC79A7", "C4"="#56B4E9",
##             "C5"="#F9A602", "C6"="#F9A602", "C7"="#F9A602")
 
col_cl <- c("C0"="#FF7F00",
   "C3"=darken("#FF7F00", 0.1), "C5"=lighten("#FF7F00", 0.2), "C7"=lighten("#FF7F00",0.5),
   "C8"=lighten("#FF7F00", 0.5), "C9"=lighten("#FF7F00", 0.5),
   "C1"="#E6E600", "C2"="#4DAF4A", "C4"="#984EA3", "C6"="#D97986", 
   "C10"="#D4B9DA")


#C6="#D97986"
## col_cl <- c("C0"="#FF7F00", "C1"="#E6E600", "C2"="#4DAF4A",
##       "C3"="#984EA3", "C4"="#AA4B56",
##        "C5"=darken("#FF7F00", 0.1), "C6"=lighten("#FF7F00", 0.2), "C7"=lighten("#FF7F00",0.5),  
##        "C8"="#D4B9DA", "C9"=lighten("#AA4B56", 0.5))


### p1
pp <- ggplot(df0, aes(x=umap_1, y=umap_2, colour=cluster_new))+
   rasterise(geom_point(size=0.1), dpi=300)+
   scale_color_manual(values=col_cl, guide="none")+  
   xlab("atacUMAP_1")+
   ylab("atacUMAP_2")+ 
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text=element_text(size=10))

###
### save figures
##
figfn <- paste(outdir, "Figure1.0_atac0.12.umap.pdf", sep="")
ggsave(figfn, pp, width=4, height=3.8)





###
### END 








