###
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(Signac)
library(SeuratWrappers)
library(SeuratData)

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.1000genomes.hs37d5, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
#library(ChIPseeker, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(ChIPseeker,lib="/wsu/home/hd/hd37/hd3768/Bin/Rpackages")

###
library(ggplot2)
#library(cowplot, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(ComplexHeatmap)
library(circlize)
#library(ggtext, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(glue)
theme_set(theme_grey())


# set the out put
outdir <- "/wsu/home/groups/piquelab/SCAIP_2022/Ali/19-atac-motif-degs"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



    
#######################################
### TFs enrichment analysis for DEG ###
#######################################


### fisher test
cal.fisher <- function(df){
   ###
   resfisher <- map_dfr(1:nrow(df),function(i){
      dmat <- matrix(as.numeric(df[i,]),2,2)
      colnames(dmat) <- c("interest", "not.interest")
      rownames(dmat) <- c("in.motif", "not.motif")
      res <- fisher.test(dmat, alternative="greater")
      res2 <- data.frame(odds=res$estimate, pval.fisher=res$p.value)
      res2
   })
   resfisher
}



###################################
###  motifs are enriched in DEG ###
###################################
##

### Results Of Deg
fn <- "/nfs/rprdata/julong/SCAIP/analyses/SCAIP-B1-6_2020.03.23/6_DEG.CelltypeNew_output/Filter2/2_meta.rds"
res <- read_rds(fn)
res <- res%>%as.data.frame()%>%drop_na(beta, qval)%>%
   mutate(direction=ifelse(beta>0, "Up", "Down"),
          comb=paste(MCls, contrast, direction, sep="_"),
          is_sig=ifelse(qval<0.1&abs(beta)>0.5, 1, 0))


#### load results of gene expression analysis

resNK <- read.table("/wsu/home/groups/piquelab/SCAIP_2022/Ali/09-ALOFT-SCAIP-batch-removed/00-ALOFT-SCAIP-pseudobulk-CTRL/out_data_SCAIP_DESeq_measured_variables_4-27-22/pnsi/stats/SCAIP_DESeq_measured_variables_4-27-22_pnsi_NKcells_DESeq_results.txt" )

resNK$MCls <- "NKcell"

resNK <- resNK %>% 
            drop_na(padj)%>%
            mutate(is_sig=ifelse(padj<0.1, 1, 0),
                    gene=gsub("\\..*", "", rownames(.)))

# in all cell types there many p-values with NAs. about 8k padj=NA

# per cell type
rn <- gsub("\\..*", "", rownames(resNK))
resNK$gene <- rn
sigNK <- resNK %>% dplyr::filter(padj<0.1) 


# load motif-gene annotation

# motif-gene list
fn <- "/nfs/rprdata/julong/sc-atac/analyses.2021-02-05/3_motif/4.2_motif.enrich.outs/Motif_gene_DAR/Motif_genes.txt"
genesTF <- read.table(fn, header=T)
dim(genesTF)
head(genesTF)

tfNK <- genesTF %>% filter(MCls=="NKcell")

anno <- tfNK
res <- resNK

comb <- "NKcell"


DF_enrich <- NULL
for (ii in comb){
   ###
   ###
   #oneMCl <- gsub("_.*", "", ii)
   oneMCl <- ii
   BG <- res%>%dplyr::filter(MCls==oneMCl)%>%pull(gene)%>%unique()
   BG <- intersect(BG, anno$gene) 

    
   DEG <- res%>%dplyr::filter(MCls==oneMCl, is_sig==1)%>%pull(gene)%>%unique() 
   DEG <- intersect(DEG, BG) #gives 17 DEGs that match (out of the 38)
   notDEG <- setdiff(BG, DEG)
    
   anno2 <- anno%>%dplyr::filter(MCls==oneMCl)
   #motifList <- unique(anno2$motif_ID)
    
    # get unique motifs that have more than 100 genes
    motif_counts <- anno2 %>%
      group_by(motif_ID) %>%
      summarise(num_genes = n_distinct(gene)) %>%
      filter(num_genes > 100) %>%
      arrange(desc(num_genes))
    motifList <-  motif_counts$motif_ID #623 

   ##

   DF2 <- map_dfr(motifList, function(ID){
       ##
       anno3 <- anno2%>%dplyr::filter(motif_ID==ID, gene%in%BG)
       TF_gene <- anno3$gene
       interest_in.motif <- length(intersect(TF_gene, DEG))
       interest_not.motif <- length(DEG)-interest_in.motif
       not.interest_in.motif <- length(intersect(notDEG, TF_gene)) 
       not.interest_not.motif <- length(notDEG)-not.interest_in.motif
       ##
       dvec <- c(interest_in.motif, interest_not.motif, not.interest_in.motif, not.interest_not.motif)
       dmat <- matrix(dvec, 2, 2)
       colnames(dmat) <- c("interest", "not.interest")
       rownames(dmat) <- c("in.motif", "not.motif")
      
       test <- fisher.test(dmat, alternative="greater", conf.int=F)
       test0 <- fisher.test(dmat)
       CI <- test0$conf.int

       ##
       df <- data.frame(comb=ii, motif_ID=ID,
           "interest_in.motif"=interest_in.motif, "interest_not.motif"=interest_not.motif,
           "not.interest_in.motif"=not.interest_in.motif, "not.interest_not.motif"=not.interest_not.motif,
           odds=test$estimate, pval=test$p.value, CI_lower=CI[1], CI_upper=CI[2])
       df
   })

   DF_enrich <- rbind(DF_enrich, DF2) 
    #DF_enrich$qvalue <- p.adjust(DF_enrich$pval, "BH")

   cat(ii, "\n") 
}    

DF_enrich2 <- DF_enrich%>%group_by(comb)%>%mutate(qval=p.adjust(pval, "BH"))%>%ungroup()

#opfn <- paste(outdir, "1_enrich_DEG.rds", sep="")
#write_rds(DF_enrich2, opfn)

sum(DF_enrich2$qval<0.1) # 0
sum(DF_enrich2$pval<0.01) # 5



#enrich <- enrich.motif(atac, anno2, feature=feature0, resSig, resDP)
#opfn <- paste(outdir, "1.", i, "_motif.DEG.direction.", feature0, ".rds", sep="")
#write_rds(enrich, opfn)










#################
### dot plots ###
#################


## read data for dot plots
feature0 <- "Allpeaks"
i <- 3
fn <- paste(outdir, "1.", i, "_motif.DEG.direction.", feature0, ".rds", sep="")
enrich <- read_rds(fn)%>%
   mutate(MCls=gsub("_.*", "", comb),
          contrast=gsub(".*[le]_|_[UD].*", "", comb),
          direction=gsub(".*_", "", comb),
          cluster=paste(contrast, MCls, direction, sep="."))%>%as.data.frame()


cl <- paste( rep(c("LPS", "PHA", "LPS+DEX", "PHA+DEX"),each=4),
   rep(c("Bcell", "Monocyte", "NKcell", "Tcell"), times=4), sep=".")
cl <- paste(rep(cl,times=2), rep(c("Up","Down"), each=16), sep=".")

cluster2 <- 1:32
names(cluster2) <- cl


## ii <- paste(rep(c("A","B","C","D"),each=4),rep(1:4,times=4), sep="")
## cl2 <- paste(rep(c("X","Y"),each=16), rep(ii,times=2), sep=".")
## cluster2 <- setNames(cl2, cl)
## lab2 <- setNames(gsub("-","+",cl),cl2)
col1 <- c("LPS"="#fb9a99", "LPS+DEX"="#e31a1c", "PHA"="#a6cee3", "PHA+DEX"="#1f78b4")
col2 <- c("Bcell"="#4daf4a", "Monocyte"="#984ea3", "NKcell"="#aa4b56", "Tcell"="#ffaa00")


enrich <- enrich%>%
  mutate(cluster=gsub("-", "+", cluster),
         contrast=gsub("-", "+", contrast),
         col.contrast=col1[contrast],
         col.MCls=col2[MCls],
         ClusterValue=as.numeric(cluster2[cluster]),
         ClusterNew=glue("<i style='color:{col.contrast}'>{contrast}.<i style='color:{col.MCls}'>{MCls}.<i style='color:black'>{direction}"),
         ClusterNew=fct_reorder(ClusterNew, ClusterValue))




## top 5 motifs
enrich2 <- enrich%>%
   dplyr::filter(fold.enrichment>1, qvalue.hyper<0.1)%>%
   group_by(ClusterNew)%>%
   top_n(n=5, wt=fold.enrichment)%>%ungroup()

topmotif <- enrich2%>%dplyr::pull(motif)
topmotif <- unique(topmotif)

enrich3 <- enrich%>%
   dplyr::filter(motif%in%topmotif, fold.enrichment>1, qvalue.hyper<0.1) 

###
p <- ggplot(enrich3, aes(x=ClusterNew, y=motif.name))+
   geom_point(aes(size=fold.enrichment, colour=qvalue.hyper))+
   scale_colour_gradient(name="FDR",
      low="blue", high="red", na.value=NA, trans="reverse", n.breaks=5,
      guide=guide_colourbar(order=1))+
   scale_size_binned("fold enrichment",
      guide=guide_bins(show.limits=T, axis=T,
      axis.show=arrow(length=unit(1.5,"mm"), ends="both"), order=2),
      n.breaks=4)+
   ## scale_x_discrete(labels=lab2)+
   theme_bw()+
   theme(axis.title=element_blank(),
         axis.text.x=element_markdown(angle=45, hjust=1, size=10),
         axis.text.y=element_text(size=8))
###
figfn <- paste(outdir, "Figure1.", i, "_dotplot.direction.", feature0, ".png", sep="")
png(figfn, width=1000, height=1000, res=120)
print(p)
dev.off()
 

