require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
library(annotables)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(plyr); library(dplyr)
library(parallel)
library("AnnotationHub")
library(ggseurat)
library(cowplot)
library(sva)
library(ggpubr)

library(Matrix)
library(future)
library(readr)
library(gtools)
library(RColorBrewer)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(ggrepel)
library(ggrastr)
library(scales)
library(openxlsx)




rm(list=ls())
setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/umap/")
#outFolderH="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/umap/"
outFolderH="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/umap/"

args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.1) #for testingbase <- args[1]

base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
outFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")

#schold <- read_rds("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/ALL.0.2.11.wavefilt.bticfilt.seurat.RDS")
#schold <- read_rds("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/ALL.0.15.13.bticfilt.seurat.RDS")
#sc <- schold
#metahold <- sc@meta.data
#metahold <- schold@meta.data

metahold$celltype <- "NA"
metahold$celltype[metahold$RNA_snn_res.0.1 == "0"] <- "R0 T CD4+"
metahold$celltype[metahold$RNA_snn_res.0.1 == "1"] <- "R1 T CD8+"
metahold$celltype[metahold$RNA_snn_res.0.1 == "2"] <- "R2 NK"
metahold$celltype[metahold$RNA_snn_res.0.1 == "3"] <- "R3 Monocyte"
metahold$celltype[metahold$RNA_snn_res.0.1 == "4"] <- "R4 B"
metahold$celltype[metahold$RNA_snn_res.0.1 == "5"] <- "R5 DC"
metahold$celltype[metahold$RNA_snn_res.0.1 == "6"] <- "R6 Other"


opfn <- paste0(outFolder,project,".",resset,".",dimset,".bticfilt.seurat.RDS")
schold <- read_rds(opfn)
metahold <- schold@meta.data

dim(schold)
table(metahold$RNA_snn_res.0.1)


outdir=paste0(base,"5b_IdenCelltype_",method,"/noDEX/")
opfn <- paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
sc <- read_rds(opfn)
meta <- sc@meta.data

dim(sc)
dim(meta)
table(meta$RNA_snn_res.0.1)


meta$celltype <- "NA"
meta$celltype[meta$RNA_snn_res.0.1 == "0"] <- "R0 T CD4+"
meta$celltype[meta$RNA_snn_res.0.1 == "1"] <- "R1 T CD8+"
meta$celltype[meta$RNA_snn_res.0.1 == "2"] <- "R2 NK"
meta$celltype[meta$RNA_snn_res.0.1 == "3"] <- "R3 Monocyte"
meta$celltype[meta$RNA_snn_res.0.1 == "4"] <- "R4 B"
meta$celltype[meta$RNA_snn_res.0.1 == "5"] <- "R5 DC"
meta$celltype[meta$RNA_snn_res.0.1 == "6"] <- "R6 Other"


#fname= paste0(outFolderH, "meta_data_scHOLD_seurat.bticfilt.071125.txt")
fname <- paste0(outdir,project,"meta_data_harmony-sctype-",resset,".",dimset,".txt")
write.table(meta, fname, sep = "\t", row.names = FALSE, quote = FALSE) 

# write out the umap coordinates
umap_df <- sc@reductions$umap@cell.embeddings
umap_df <- as.data.frame(umap_df)
umap_df$barcode <- rownames(umap_df)

fname <- paste0(outdir,project,"UMAP_harmony-sctype-",resset,".",dimset,".txt")
write.table(umap_df, fname, sep = "\t", row.names = FALSE, quote = FALSE) 

col_cl <- c(
  "R0 T CD4+"   = "#FF7F00",  
  "R1 T CD8+"   = "#E6E600",  
  "R2 NK"       = "#4DAF4A",  
  "R3 Monocyte" = "#984EA3",  
  "R4 B"        = "#D97986",  #AA4B56
  "R5 DC"       = "#D4B9DA",  
  "R6 Other"   = "#999999"  
)


x <- meta
plotDF <- as.data.frame(sc[["umap"]]@cell.embeddings)
#identical(rownames(plotDF), rownames(x))
plotDF$Clusters <- paste("C", x$RNA_snn_res.0.1, sep="")
plotDF$celltype <- x$celltype
plotDF$treat <- x$treat
plotDF$BATCH <- x$BATCH


### p2
p2 <- ggplot(plotDF, aes(x=umap_1, y=umap_2, colour=celltype))+
   rasterise(geom_point(size=0.05), dpi=300)+
   scale_color_manual(values=col_cl)+#, guide="none")+  
      guides(colour = guide_legend(override.aes = list(size = 3))) +
   #geom_text(data = label_positions, aes(label = celltype), 
   #         color = "black", size = 3, fontface = "bold") +
   xlab("rnaUMAP_1")+
   ylab("rnaUMAP_2")+ 
   theme_bw()+
   theme(axis.title=element_text(size=10),
         axis.text=element_text(size=10), 
        legend.title = element_text(size = 10),
      legend.text = element_text(size = 9), 
      legend.position = "right")

###
### save figures
##
figfn <- paste(outFolderH, "Figure1.1_rna0.2.umap.labled.v4.png", sep="")
png(figfn,width=1300,height=1000, res=240)
print(p2)
dev.off()

figfn <- paste(outFolderH, "Figure1.1_rna0.2.umap.labled.v4.pdf", sep="")
pdf(figfn, width=5, height=4)
print(p2)
dev.off()



############################# facet wrap by treatment
p2 <- ggplot(plotDF, aes(x = umap_1, y = umap_2, colour = celltype)) +
  rasterise(geom_point(size = 0.1), dpi = 600) +
  #geom_point(size = 0.1) +
  scale_color_manual(values = col_cl) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ treat, nrow = 1) +  
  xlab("rnaUMAP_1") +
  ylab("rnaUMAP_2") + 
  #theme_bw(base_size = 20) +
  theme_bw(base_size = 15) +

  theme(
    #axis.title = element_text(size = 10),
    #axis.text = element_text(size = 10), 
    #legend.title = element_text(size = 10),
    #legend.text = element_text(size = 9), 
    legend.position = "right"
  )

#figfn <- paste(outFolderH, "Figure1.2_rna0.2.umap.facet.treat.png", sep="")
figfn <- paste(outFolderH, "Figure1.3_rna0.2.umap.facet.treat.highres.png", sep="")
#png(figfn,width=3000,height=1000, res=240)
png(figfn,width=2000,height=1000, res=240)
#png(figfn,width=6000,height=3000, res=360)
print(p2)
dev.off()


figfn <- paste(outFolderH, "Figure1.3_rna0.2.umap.facet.treat.highres.pdf", sep="")
pdf(figfn, width=14, height=7)
print(p2)
dev.off()


############################ facet wrap by batch 
plotDF$BATCH <- factor(plotDF$BATCH, levels= c("HOLD1", "HOLD2", "HOLD3", "HOLD4", "HOLD5", "HOLD6", "HOLD7", "HOLD8", "HOLD9", "HOLD10", "HOLD11", "HOLD12", "HOLD13", "HOLD14"))

p2 <- ggplot(plotDF, aes(x = umap_1, y = umap_2, colour = celltype)) +
  #rasterise(geom_point(size = 0.1), dpi = 600) +
  geom_point(size = 0.1) +
  scale_color_manual(values = col_cl) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ BATCH, nrow = 2) +  
  xlab("rnaUMAP_1") +
  ylab("rnaUMAP_2") + 
  theme_bw(base_size = 20) +
  theme(
    #axis.title = element_text(size = 10),
    #axis.text = element_text(size = 10), 
    #legend.title = element_text(size = 10),
    #legend.text = element_text(size = 9), 
    legend.position = "right"
  )

#figfn <- paste(outFolderH, "Figure1.4_rna0.2.umap.facet.batch.png", sep="")
figfn <- paste(outFolderH, "Figure1.5_rna0.2.umap.facet.batch.highres.png", sep="")
png(figfn,width=7000,height=3000, res=360)
print(p2)
dev.off()


figfn <- paste(outFolderH, "Figure1.5_rna0.2.umap.facet.batch.highres.pdf", sep="")
pdf(figfn, width=12, height=5)
print(p2)
dev.off()




################################################################################################## STATs and plots ########################################################################
# stats will be reported prior to >1000 cells and >20comb filtering for psudobulking, but after the >10%mt filtering, >200 geanefeature filtering. 
# this is the seurat object used: 
#outdir=paste0(base,"5b_IdenCelltype_",method,"/noDEX/")
#opfn <- paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
#sc <- read_rds(opfn)
#meta <- sc@meta.data

library(dplyr)
# library stats
#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta %>% group_by(Library)%>%
             summarise(ncell=n(),
                       reads=mean(nCount_RNA),
                       ngene=mean(nFeature_RNA),
                       percent.mt=mean(percent.mt)#,
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       #.groups="drop"
                       )
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))


# library stats
sum(dd$ncell) # 575885
median(dd$ncell) # 21684
median(dd$reads) # 8193.465
median(dd$ngene) # 2645.313
dim(dd) # [1] 28  


# treatment stats
treatment_stats <- meta %>% 
  group_by(treats) %>% 
  summarise(
    ncell = n(),
    median_cell = median(ncell),
    median_reads = median(nCount_RNA),
    median_ngene = median(nFeature_RNA)
  )


# wide format for number of cells
# Assuming `metahold@meta.data` is your metadata
metasc2 <- meta
metasc2 <- metasc2 %>%
  mutate(BATCH = factor(BATCH, levels = mixedsort(unique(BATCH))))

metadata <- meta

# Count unique Sample_IDs for each batch and treatment
sample_counts <- metasc2 %>%
  group_by(Sample_ID, treats) %>%
  summarise(#nSamples = n_distinct(Sample_ID),
            nCells = n(), .
            #groups = "drop"
            )

sample_counts_wide <- sample_counts %>%
  pivot_wider(names_from = treats, values_from = nCells, values_fill = 0)

sample_counts_wide <- as.data.frame(sample_counts_wide)

fname <- paste0(outFolderH, "/10.1_cells_per_sample_per_treatment.txt")
write.table(sample_counts_wide, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


stats <- metadata %>% 
  group_by(treats) %>% 
  summarise(
    n_unique_samples = n_distinct(Sample_ID),
    median_cells_per_sample = median(table(Sample_ID)),
    median_reads_per_cell = median(nCount_RNA),
    median_genes_per_cell = median(nFeature_RNA)
  ) %>% as.data.frame()

fname <- paste0(outFolderH, "10.2_stats_cells_UMI_genes.txt")
write.table(stats, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

metadata$celltype <- factor(metadata$celltype, levels= c("R0 T CD4+", "R1 T CD8+", "R2 NK", "R3 Monocyte", "R4 B", "R5 DC", "R6 Other"))
metadata$BATCH <- factor(metadata$BATCH, levels= c("HOLD1", "HOLD2", "HOLD3", "HOLD4", "HOLD5", "HOLD6", "HOLD7", "HOLD8", "HOLD9", "HOLD10", "HOLD11", "HOLD12", "HOLD13", "HOLD14"))

cell_medians_treatment <- metadata %>%
  group_by(treats, Sample_ID, celltype) %>% 
  summarise(cell_count = n(), .groups = "drop") %>%  # Count cells per sample per cell type
  group_by(treats, celltype) %>% 
  summarise(median_cells_per_sample = median(cell_count), .groups = "drop") %>% # Median per treatment-celltype
  pivot_wider(names_from = celltype, values_from = median_cells_per_sample, values_fill = 0) %>% as.data.frame()

fname <- paste0(outFolderH, "/10.3_stats_ncells_per_celltype_per_sample_per_treatmet.txt")
write.table(cell_medians_treatment, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

################################################################ bar plots
# bar plots for number of cells and barcodes and genes, etc

dd$ident <- factor(dd$ident, levels= c(
                "HOLD1-RNA-CTRL", "HOLD1-RNA-LPS" ,"HOLD2-RNA-CTRL", "HOLD2-RNA-LPS", "HOLD3-RNA-CTRL", "HOLD3-RNA-LPS", "HOLD4-RNA-CTRL", "HOLD4-RNA-LPS", 
                "HOLD5-RNA-CTRL", "HOLD5-RNA-LPS", "HOLD6-RNA-CTRL", "HOLD6-RNA-LPS", "HOLD7-RNA-CTRL", "HOLD7-RNA-LPS", "HOLD8-RNA-CTRL", 
                "HOLD8-RNA-LPS","HOLD9-RNA-CTRL",  "HOLD9-RNA-LPS", "HOLD10-RNA-CTRL", "HOLD10-RNA-LPS", "HOLD11-RNA-CTRL",
                "HOLD11-RNA-LPS", "HOLD12-RNA-CTRL", "HOLD12-RNA-LPS", "HOLD13-RNA-CTRL", "HOLD13-RNA-LPS", 
                "HOLD14-RNA-CTRL", "HOLD14-RNA-LPS"
                ))
dd$batch <- factor(dd$batch, levels= c("HOLD1", "HOLD2", "HOLD3", "HOLD4", "HOLD5", "HOLD6", "HOLD7", "HOLD8", "HOLD9", "HOLD10", "HOLD11", "HOLD12", "HOLD13", "HOLD14"))



###(1), barcodes for each experiment           
fig0 <- ggplot(dd,aes(x=ident, y=ncell, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("number of cells per library")+
        geom_text(aes(label=ncell),position = position_dodge(0.9), vjust = 0.5, size=2.5, angle=90, hjust=0)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  
fname=paste0(outFolderH, "Figure7.1_barcodes_post_filt",".png")
png(fname, width=1000, height=600, res=120)
print(fig0)
dev.off()

figfn <- paste(outFolderH, "Figure7.1_barcodes_post_filt",".pdf", sep="")
pdf(figfn, width=10, height=5)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads)
dd2 <- dd%>%dplyr::select(ident,ngene,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="#UMIs per cell", "2"="#Genes per cell"))
fig0 <- ggplot(ddnew, aes(x=ident, y=y, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        facet_wrap(~factor(stats), nrow=2, scales="free_y", labeller=stats)+
        geom_text(aes(label=round(y)),vjust=-0.7, size=3.0)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=10),
              axis.text.y=element_text(size=10), 
              strip.text = element_text(size = 20),
              strip.background=element_blank())
##   
fname=paste0(outFolderH, "/Figure7.2_genes_post_filt.png")     
png(fname, width=3000, height=2500, res=240)
print(fig0)
dev.off() 

figfn <- paste(outFolderH, "Figure7.2_genes_post_filt",".pdf", sep="")
pdf(figfn, width=12, height=10)
print(fig0)
dev.off()


## separate plots for reads and gene numbers
fig0 <- ggplot(dd,aes(x=ident, y=reads, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("#UMI per cell")+
        geom_text(aes(label=round(reads)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

fname=paste0(outFolderH, "/Figure7.3_UMI_numb_post_merge_post_filt.png")     
png(fname, width=1000, height=600, res=120)
print(fig0)
dev.off()

fig0 <- ggplot(dd,aes(x=ident, y=ngene, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("#Gene per cell")+
        geom_text(aes(label=round(ngene)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

fname=paste0(outFolderH, "/Figure7.4_Gene_numb_post_filt.png")     
png(fname, width=1000, height=600, res=120)
print(fig0)
dev.off()


# mitochondrial content
fig0 <- ggplot(dd,aes(x=ident, y=percent.mt, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("percent mitochondria (mean)")+
        geom_text(aes(label=round(percent.mt)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  
fname=paste0(outFolderH, "/Figure7.5_percent_mt_merge_post_filt.png")     
png(fname, width=1000, height=600, res=120)
print(fig0)
dev.off()



# Count cells per cluster
tab <- table(meta$celltype)

# Convert to dataframe with percentages
celltype_df <- data.frame(
  celltype = names(tab),
  count = as.vector(tab)
) %>%
  mutate(percentage = round(100 * count / sum(count), 2))

celltype_df

fname <- paste0(outFolderH, "/10.4_celltype_count_percentage.txt")
write.table(celltype_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


df_summary <- dd %>%
  separate(ident, into = c("BATCH", "RNA", "Treat"), sep = "-") %>%
  select(BATCH, Treat, ncell) %>%
  pivot_wider(names_from = Treat, values_from = ncell, names_prefix = "RNA-") %>%
  arrange(BATCH)

df_summary$BATCH <- factor(df_summary$BATCH, levels= c("HOLD1", "HOLD2", "HOLD3", "HOLD4", "HOLD5", "HOLD6", "HOLD7", "HOLD8", "HOLD9", "HOLD10", "HOLD11", "HOLD12", "HOLD13", "HOLD14"))
df_summary <- df_summary[order(df_summary$BATCH), ]

df_summary


fname <- paste0(outFolderH, "/10.5_celltype_scRNA_counts_per_batch_per_treat.txt")
write.table(df_summary, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



df_summary <- meta %>%
  group_by(BATCH, treats) %>%
  summarise(
    n_individuals = n_distinct(Sample_ID),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = treats,
    values_from = n_individuals,
    values_fill = 0
  )
df_summary$BATCH <- factor(df_summary$BATCH, levels= c("HOLD1", "HOLD2", "HOLD3", "HOLD4", "HOLD5", "HOLD6", "HOLD7", "HOLD8", "HOLD9", "HOLD10", "HOLD11", "HOLD12", "HOLD13", "HOLD14"))
df_summary <- as.data.frame(df_summary[order(df_summary$BATCH), ])
df_summary

fname <- paste0(outFolderH, "/10.6_celltype_scRNA_nsamples_per_batch_per_treat.txt")
write.table(df_summary, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)





############################################# scatac UMAP ###########################################

#atacdir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/1.2_ArchR_process"
#/2_Integrate_output/new/

#infn_meta <- paste0(atacdir, "/2_Integrate_output/new/1_RNA.res0.2.meta.rds")
#df0 <- as.data.frame(read_rds(infn_meta))
#df0$rowname <- rownames(df0)

#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/1.2_ArchR_process/888_pub_figures_output/1_meta_final.rds"
#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/1.2_ArchR_process/1_results_all/Batch_correct_Harmony/1.0_cluster.rds"
#df1 <- read_rds(fname)

#fname = "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_new_2025_08_07/1.2_ArchR_process/888_pub_figures_output/1_meta_final.rds"
#df1 <- read_rds(fname)
#"/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/1.2_ArchR_process"


outFolderH="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/umap/"

outdir <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/1.2_ArchR_process/888_pub_figures_output/main/"


infn_atac <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/analyses_correct_cellRanger_2025_11_14/1.2_ArchR_process/1_results_all/Batch_correct_Harmony/1.0_cluster.rds"
df <- read_rds(infn_atac)

fn <- paste(outdir, "1_meta_final.rds", sep="")
df1 <- read_rds(file=fn)



#/1.2_ArchR_process/1_results_all/Batch_correct_Harmony/1.0_cluster.rds
rn <- rownames(df1)
df1$BATCH <- sub("-.*", "", rn)
df1$treat <- sub(".*-ATAC-([^#]+)#.*", "\\1", rn)


df3 <- df1


table(df3$cluster_new)

# assign atac metadata: 
df3$celltype <- "NA"
df3$celltype[df3$cluster_new == "C0"] <- "A0 T CD4+"
df3$celltype[df3$cluster_new == "C1"] <- "A1 T CD8+"
df3$celltype[df3$cluster_new == "C2"] <- "A2 NK"
df3$celltype[df3$cluster_new == "C3"] <- "A3 T CD4+"
df3$celltype[df3$cluster_new == "C4"] <- "A4 Monocyte"
df3$celltype[df3$cluster_new == "C5"] <- "A5 T CD4+"
df3$celltype[df3$cluster_new == "C6"] <- "A6 B"
df3$celltype[df3$cluster_new == "C7"] <- "A7 T CD4+"
df3$celltype[df3$cluster_new == "C8"] <- "A8 T CD4+"
df3$celltype[df3$cluster_new == "C9"] <- "A9 T CD4+"
df3$celltype[df3$cluster_new == "C10"] <- "A10 DC"

#df0$celltype[df0$Clusters == "7"] <- "R7 T CD4+"



#col_cl <- c("C0"="#E60A0A",  "C1"="#A65192", "C2"="#0570B0", "C3"="#00AF72", "C4"="#87CEFA",
#                        "C5"="#F9A602", "C6"="#F9A602", "C7"="#F9A602", "C8"="#993F50")


#col_cl <- c("R0 T CD4+"="#E60A0A", "R1 T CD4+"="#F9A602", "R2 T CD8+"="#A65192", "R3 NK"="#0570B0", "R4 Monocyte"="#00AF72", "R5 B"="#87CEFA",
#                         "R6 DC"="#993F50", "R7 T CD4+"="#FF7F0E")
col_cl <- c(
  "A0 T CD4+"   = "#FF7F00",  # orange
  "A1 T CD8+"   = "#E6E600",  # blue
  "A2 NK"      = "#4DAF4A",  # reddish orange
  "A4 Monocyte" = "#984EA3",  # purple/magenta
  "A6 B"       = "#D97986",  # light blue
  "A10 DC"      = "#D4B9DA",  # yellow
  "A3 T CD4+"  = "#FF7F00",   #  teal green
  "A5 T CD4+"   = "#FF7F00",   #  teal green
  "A7 T CD4+"  = "#FF7F00",   #  teal green
  "A8 T CD4+"  = "#FF7F00",   #  teal green
  "A9 T CD4+" = "#FF7F00"   #  teal green
)



## julongs colors: 
#col_cl <- c("C0"="#FF7F00",
#   "C3"=darken("#FF7F00", 0.1), "C5"=lighten("#FF7F00", 0.2), "C7"=lighten("#FF7F00",0.5),
#   "C8"=lighten("#FF7F00", 0.5), "C9"=lighten("#FF7F00", 0.5),
#   "C1"="#E6E600", "C2"="#4DAF4A", "C4"="#984EA3", "C6"="#D97986", 
#   "C10"="#D4B9DA")

col_cl <- c(
  "A0 T CD4+"   = "#FF7F00",  # orange
  "A1 T CD8+"   = "#E6E600",  # blue
  "A2 NK"      = "#4DAF4A",  # reddish orange
  "A4 Monocyte" = "#984EA3",  # purple/magenta
  "A6 B"       = "#D97986",  # light blue
  "A10 DC"      = "#D4B9DA",  # yellow
  "A3 T CD4+"  = darken("#FF7F00", 0.1),   #  teal green
  "A5 T CD4+"   = lighten("#FF7F00", 0.2),   #  teal green
  "A7 T CD4+"  = lighten("#FF7F00",0.5),   #  teal green
  "A8 T CD4+"  = lighten("#FF7F00", 0.5),   #  teal green
  "A9 T CD4+" = lighten("#FF7F00", 0.5)   #  teal green
)



#   "A0 T CD4+"   = "#FF7F00",  # orange
#  "A5 T CD4+"   = "#FF7F00",   #  teal green
#  "A6 T CD4+"   = "#FF7F00",   #  teal green
#  "A7 T CD4+"   = "#FF7F00",   #  teal green
 df3$celltype <- factor(df3$celltype, levels= c("A0 T CD4+",
  "A1 T CD8+",
  "A2 NK", 
  "A3 T CD4+",
  "A4 Monocyte",
  "A5 T CD4+",
  "A6 B",
  "A7 T CD4+",
  "A8 T CD4+",
  "A9 T CD4+", 
  "A10 DC"
))

### p1
p2 <- ggplot(df3, aes(x=umap_1, y=umap_2, colour=celltype))+
   rasterise(geom_point(size=0.1), dpi=300)+
   scale_color_manual(values=col_cl)+#, guide="none")+  
     guides(colour = guide_legend(override.aes = list(size = 3))) +
   xlab("atacUMAP_1")+
   ylab("atacUMAP_2")+ 
   theme_bw()+
   theme(axis.title=element_text(size=10),
         axis.text=element_text(size=10))

###
### save figures
##

figfn <- paste(outFolderH, "Figure2.0_atac0.umap.png", sep="")
#figfn <- paste(outFolderH, "Figure1.5_rna0.2.umap.facet.batch.highres.png", sep="")
png(figfn,width=1300,height=1000, res=240)
print(p2)
dev.off()


############################# facet wrap by treatment
p2 <- ggplot(df3, aes(x = umap_1, y = umap_2, colour = celltype)) +
  rasterise(geom_point(size = 0.1), dpi = 600) +
  #geom_point(size = 0.1) +
  scale_color_manual(values = col_cl) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ treat, nrow = 1) +  
  xlab("atacUMAP_1") +
  ylab("atacUMAP_2") + 
  theme_bw(base_size = 15) +
  theme(
    #axis.title = element_text(size = 10),
    #axis.text = element_text(size = 10), 
    #legend.title = element_text(size = 10),
    #legend.text = element_text(size = 9), 
    legend.position = "right"
  )

figfn <- paste(outFolderH, "Figure2.2_atac0.07.umap.facet.treat.png", sep="")
#figfn <- paste(outFolderH, "Figure2.3_atac0.07.umap.facet.treat.highres.png", sep="")
png(figfn,width=2000,height=1000, res=240)
#png(figfn,width=7000,height=3000, res=360)
print(p2)
dev.off()



############################ facet wrap by batch 
df3$BATCH <- factor(df3$BATCH, levels= c("HOLD1", "HOLD2", "HOLD3", "HOLD4", "HOLD5", "HOLD6", "HOLD7", "HOLD8", "HOLD9", "HOLD10", "HOLD11", "HOLD12", "HOLD13", "HOLD14"))

p2 <- ggplot(df3, aes(x = umap_1, y = umap_2, colour = celltype)) +
  rasterise(geom_point(size = 0.1), dpi = 600) +
  #geom_point(size = 0.1) +
  scale_color_manual(values = col_cl) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(~ BATCH, nrow = 2) +  
  xlab("atacUMAP_1") +
  ylab("atacUMAP_2") + 
  theme_bw(base_size = 15) +
  theme(
    #axis.title = element_text(size = 10),
    #axis.text = element_text(size = 10), 
    #legend.title = element_text(size = 10),
    #legend.text = element_text(size = 9), 
    legend.position = "right"
  )

figfn <- paste(outFolderH, "Figure2.4_atac0.07.umap.facet.batch.png", sep="")
#figfn <- paste(outFolderH, "Figure2.5_atac0.07.umap.facet.batch.highres.png", sep="")
png(figfn,width=7000,height=3000, res=360)
print(p2)
dev.off()




meta <- df3

# Count cells per cluster
tab <- table(meta$celltype)

# Convert to dataframe with percentages
celltype_df <- data.frame(
  celltype = names(tab),
  count = as.vector(tab)
) %>%
  mutate(percentage = round(100 * count / sum(count), 2))

celltype_df

fname <- paste0(outFolderH, "/10.4_atac_celltype_count_percentage.txt")
write.table(celltype_df, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

dd <- meta %>% group_by(BATCH, treat, celltype)%>%
             summarise(ncell=n()
                       #reads=mean(nCount_RNA),
                       #ngene=mean(nFeature_RNA),
                       #percent.mt=mean(percent.mt)#,
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       #.groups="drop"
                       ) %>% as.data.frame()

dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))


df_summary <- dd %>%
  separate(ident, into = c("BATCH", "treat"), sep = "-") %>%
  select(BATCH, treat, ncell) %>%
  pivot_wider(names_from = treat, values_from = ncell, names_prefix = "RNA-") %>%
  arrange(BATCH)

df_summary$BATCH <- factor(df_summary$BATCH, levels= c("HOLD1", "HOLD2", "HOLD3", "HOLD4", "HOLD5", "HOLD6", "HOLD7", "HOLD8", "HOLD9", "HOLD10", "HOLD11", "HOLD12", "HOLD13", "HOLD14"))
df_summary <- df_summary[order(df_summary$BATCH), ]

df_summary


fname <- paste0(outFolderH, "/10.5_celltype_scRNA_counts_per_batch_per_treat.txt")
write.table(df_summary, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



df_summary <- meta %>%
  group_by(BATCH, treat) %>%
  summarise(
    n_individuals = n_distinct(Sample_ID),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = treat,
    values_from = n_individuals,
    values_fill = 0
  )
df_summary$BATCH <- factor(df_summary$BATCH, levels= c("HOLD1", "HOLD2", "HOLD3", "HOLD4", "HOLD5", "HOLD6", "HOLD7", "HOLD8", "HOLD9", "HOLD10", "HOLD11", "HOLD12", "HOLD13", "HOLD14"))
df_summary <- as.data.frame(df_summary[order(df_summary$BATCH), ])
df_summary

fname <- paste0(outFolderH, "/10.6_celltype_scRNA_nsamples_per_batch_per_treat.txt")
write.table(df_summary, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



