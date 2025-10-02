##############################################
### cell type classification using SingleR ###
### reference: parturition elife           ### 
##############################################

library(Seurat)
library(Matrix)
library(tidyverse)
library(harmony)
library(pheatmap)
library(symphony)
library(ggrepel)
#################
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/",0.1,"ALL","demux",50) 
future::plan(strategy = 'multicore', workers = 7)
options(future.globals.maxSize = 100 * 1024 ^ 3)

base <- args[1]
resset <- as.numeric(args[2])
project <- args[3]
method <- args[4]
dimset=args[5]
cat("resolution=",resset,"\nproject=",project,"\n","method=",method,"\n")
#filter <- "noDEX"
filter <- "CTRLonly" #ALOFT used

indir=paste0(base,"5b_IdenCelltype_",method,"/",filter,"/")

outdir1=paste0(base,"5b_IdenCelltype_",method,"symphonyTBRU/")
if (!file.exists(outdir1)) dir.create(outdir1, showWarnings=F)

outdir=paste0(base,"5b_IdenCelltype_",method,"symphonyTBRU/",filter,"/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outdir,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#anno <- read_rds("/wsu/home/groups/prbgenomics/preeclampsia/preeclampsia_sc/3_runHarmony/anno-after-mt-filtering.2023-05-07.rds")

#cell.type.anno <- read_tsv("cell.type.anno.res0p5_2023-05-07.tsv")
#cell.type.anno$seurat_clusters <- factor(cell.type.anno$seurat_clusters)

opfn <- paste0(indir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
sc2 <- readRDS(opfn)

ref <- read_rds("/wsu/el7/groups/piquelab/refData/symphony/tbru_ref.rds")
ref$save_uwot_path <- "/wsu/el7/groups/piquelab/refData/symphony/tbru_uwot_model"

rmd <- read_tsv("/wsu/el7/groups/piquelab/refData/symphony/tbru_metadata.txt.gz")
ref$meta_data <- rmd

md <- sc2@meta.data

# Map query
query = mapQuery(sc2@assays$RNA@counts,             # query gene expression (genes x cells)
                 md,        # query metadata (cells x attributes)
                 ref,             # Symphony reference object
                 vars = c("EXP"),           # Query batch variables to harmonize over (NULL treats query as one batch)
                 verbose = "TRUE",
                 do_normalize = TRUE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                 do_umap = TRUE)        # project query cells into reference UMAP


p = plotReference(ref,
                  as.density = TRUE,      # plot density or individual cells
                  bins = 14,              # if density, nbins parameter for stat_density_2d
                  bandwidth = 1,        # if density, bandwidth parameter for stat_density_2d
                  title = "Symphony Reference_tbru ref",    # Plot title
                  color.by = 'cluster_name', # metadata column name for cell type labels
           ##       celltype.colors = pbmc_colors, # custom color palette
                  show.legend = TRUE,     # Show cell type legend
                  show.labels = TRUE,     # Show cell type labels
                  show.centroids = FALSE) # Plot soft cluster centroid locations)
png(width = 9, height = 8, file=paste0(figuredir,"plotReference.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()


tlabels <-  ref$meta_data$cluster_name

tlabels[is.na(tlabels)] <- "Unk"

# Predict query cell types
query = knnPredict(query, ref, 
                   train_labels = tlabels,
                   k = 30, 
                   confidence = TRUE) # calculate prediction confidence

#####

md2 <- cbind(query$meta_data,query$umap)

p2 = md2 %>%
  sample_frac(1L) %>% # permute rows randomly
  ggplot(aes(x = UMAP1, y = UMAP2, col = cell_type_pred_knn)) +
  geom_point(size = 0.4, stroke = 0.2, shape = 16) +
  theme_bw() +
  labs(title = 'Query after Symphony: Colored by cell type', color = '') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="right", legend.text = element_text(size=9), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 3))) 


centroids <- md2 %>%
  group_by(cell_type_pred_knn) %>%
  summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

p2_centroids = md2 %>%
  sample_frac(1L) %>% # permute rows randomly
  ggplot(aes(x = UMAP1, y = UMAP2, col = cell_type_pred_knn)) +
  geom_point(size = 0.4, stroke = 0.2, shape = 16) +
  stat_ellipse(aes(group = cell_type_pred_knn, color = cell_type_pred_knn), 
               type = "norm", level = 0.68, linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = 'Query after Symphony: Colored by cell type', color = '') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.text = element_text(size = 9), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 14)) +
  guides(color = guide_legend(override.aes = list(size = 3))) + # Increase legend symbol size
  geom_text(data = centroids, aes(x = UMAP1, y = UMAP2, label = cell_type_pred_knn), 
            color = "black", size = 3)


png(width = 9, height = 8, file=paste0(figuredir,"queryaftersymphony_celltypecolor.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
print(p2)
dev.off()

png(width = 9, height = 8, file=paste0(figuredir,"queryaftersymphony_celltypecolor_centroids.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
print(p2_centroids)
dev.off()


#not sure what this does???
#p / p2


write.csv(md2,paste0(outdir,"sc2_metadata_symphony_tbru.csv"))

tt_knn <- table(md2$cell_type_pred_knn,md2$seurat_clusters)  #abundance of cell types in each cluster
tt_knn
write.csv(tt_knn,paste0(outdir,"symphony_tbru_celltype_counts_per_clusters.csv"))

saveRDS(query,paste0(outdir,"SymphonyQuery.TBRU.rds"))


png(width = 8, height = 8, file=paste0(figuredir,"heatmap_celltype_Symphony_TBRU.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
p <- pheatmap(tt_knn,cluster_rows=TRUE,cluster_cols=FALSE,scale = "column")
print(p)
dev.off()

##########################################################################

#from not tbru script
outdir1=paste0(base,"5b_IdenCelltype_",method,"symphonyPBMC/")
if (!file.exists(outdir1)) dir.create(outdir1, showWarnings=F)

outdir=paste0(base,"5b_IdenCelltype_",method,"symphonyPBMC/",filter,"/")
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outdir,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

ref <- read_rds("/wsu/el7/groups/piquelab/refData/symphony/pbmcs_10x_reference.rds")
ref$save_uwot_path <- "/wsu/el7/groups/piquelab/refData/symphony/pbmcs_10x_uwot_model"

# Map query
query = mapQuery(sc2@assays$RNA@counts,             # query gene expression (genes x cells)
                 md,        # query metadata (cells x attributes)
                 ref,             # Symphony reference object
                 vars = c("EXP"),           # Query batch variables to harmonize over (NULL treats query as one batch)
                 do_normalize = TRUE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                 do_umap = TRUE)        # project query cells into reference UMAP


p = plotReference(ref,
                  as.density = TRUE,      # plot density or individual cells
                  bins = 14,              # if density, nbins parameter for stat_density_2d
                  bandwidth = 1,        # if density, bandwidth parameter for stat_density_2d
                  title = "Symphony Reference PBMC",    # Plot title
                  color.by = 'pub_cell_subtype', # metadata column name for cell type labels
           ##       celltype.colors = pbmc_colors, # custom color palette
                  show.legend = TRUE,     # Show cell type legend
                  show.labels = TRUE,     # Show cell type labels
                  show.centroids = FALSE) # Plot soft cluster centroid locations)

png(width = 9, height = 8, file=paste0(figuredir,"plotReference.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

# Predict query cell types
query = knnPredict(query, ref, 
                   train_labels = ref$meta_data$pub_cell_subtype,
                   k = 30, 
                   confidence = TRUE) # calculate prediction confidence

md2 <- cbind(query$meta_data,query$umap)

p2 = md2 %>%
  sample_frac(1L) %>% # permute rows randomly
  ggplot(aes(x = UMAP1, y = UMAP2, col = cell_type_pred_knn)) +
  geom_point(size = 0.3, stroke = 0.2, shape = 16) +
  theme_bw() +
  labs(title = 'Query after Symphony: Colored by cell type', color = '') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none", legend.text = element_text(size=11), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 14)) 

png(width = 9, height = 8, file=paste0(figuredir,"queryaftersymphony_celltypecolor.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
print(p2)
dev.off()

png(width = 9, height = 8, file=paste0(figuredir,"queryaftersymphony_celltypecolor_centroids.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
print(p2_centroids)
dev.off()

write.csv(md2,paste0(outdir,"sc2_metadata_symphony_pbmc.csv"))

tt <- table(md2$cell_type_pred_knn,md2$seurat_clusters)  #abundance of cell types in each cluster
tt
write.csv(tt,paste0(outdir,"SymphonyPBMC_classification_SeuratCluster_.csv"))

png(width = 8, height = 8, file=paste0(figuredir,"heatmap_celltype_Symphony_PBMC.png"), 
  pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 600)
p <- pheatmap(tt,cluster_rows=TRUE,cluster_cols=FALSE,scale = "column")
print(p)
dev.off()

saveRDS(query,paste0(outdir,"SymphonyQuery.PBMC.rds"))

#updating celltype labels on UMAP given types from sctype and symphony
#changing C1 to CD4+T (CD4+ CD27+), C7 to CD4+ Tcell (either naive or lncrna)
sc2@meta.data <- transform(sc2@meta.data, final_celltype=ifelse(seurat_clusters=="1", "CD4+ CD27+ T cells", ifelse(seurat_clusters=="7","CD4+ T cells",ifelse(customclassif=="Macrophages","Monocytes",customclassif))))

 [1] "Naive CD4+ T cells"    "Natural killer  cells" "CD4+ CD27+ T cells"
 [4] "CD8+ NKT-like cells"   "Pre-B cells"           "CD4+ T cells"
 [7] "Classical Monocytes"   "Memory CD4+ T cells"   "Monocytes"
[10] "γδ-T cells"            "Plasma B cells"

suboutFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/",filter,"/")
figuredir=paste0(suboutFolder,"figures/")

my_cols <- c('Naive CD4+ T cells'='#31C53F','Natural killer  cells'='#F68282','CD4+ CD27+ T cells'='#1FA195',
  'CD8+ NKT-like cells'='#ff9a36','Pre-B cells'='#E6C122', 'CD4+ T cells'='#25aff5','Classical Monocytes'='#B95FBB',
  'Memory CD4+ T cells'='midnightblue','Monocytes'='purple4','γδ-T cells'='darkgreen',
  'Plasma B cells'='magenta4')

  png(width = 9, height = 8, file=paste0(figuredir,project,".",resset,".",dimset,".harmony_umap_QC_finalcelltype.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
  p <- DimPlot(sc2, reduction = "umap", label = TRUE, repel=TRUE, group.by = 'final_celltype', pt.size = .8,
      cols = my_cols)+ labs(title = NULL)+
  theme(panel.background = element_rect(fill="white",colour = "black"))
  print(p)
  dev.off()

my_cols <- c('0'='#31C53F','1'='#1FA195','2'='#ff9a36',
  '3'='#F68282','4'='#E6C122', '5'='purple4','6'='darkgreen',
  '7'='#25aff5','8'='#B95FBB','9'='midnightblue',
  '10'='magenta4')
fname=paste0(figuredir,project,".Figure5.4_UMAP_Harmony-res",resset,".",dimset,"_group_seurat_cluster_with_names",".png");
png(width = 8, height = 8, file=fname, res=600, pointsize=12,bg = "transparent", units = "in")
fig1 <- DimPlot(sc2, reduction = "umap", label=T, group.by="seurat_clusters", label.size=10,pt.size=0.8,cols = my_cols)+
  theme(legend.position = "none",panel.background = element_rect(fill="white",colour = "black"))+ labs(title = NULL)
print(fig1)
dev.off()


cat("saving \n")
opfn <- paste0(indir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
write_rds(sc2, opfn)

    df <- unique(data.frame(seurat_clusters=unique(sc2@meta.data$seurat_clusters),cluster=unique(paste0("C",sc2@meta.data$seurat_clusters)),cell_type=unique(sc2@meta.data$final_celltype)))
df <- df[order(df$seurat_clusters),-1]
fwrite(df, sep='\t', quote=F, row.names=F, col.names=T, paste0(base,method,"_pseudobulk_ctrl/",filter,"/",project,".",resset,".",dimset,".cluster_celltype.txt"))
