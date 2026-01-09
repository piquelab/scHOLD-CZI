##
library(Matrix)
library(tidyverse)
library(SummarizedExperiment)
###

###
### Generate pseudubulk motif activity exclude X chrosome 
option <- "option_nFeature15K_cluster_res0.12"
outdir <- paste("./2_motif.outs/", option, "/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T) 


###
### 



YtX_all <- NULL
dd_all <- NULL
  

###
### read data
infn <- "../1.2_ArchR_process/4b_motif_output/Motif_jaspar2022.activity.mat.rds"
mat <- read_rds(infn)

    
###
### column
meta <- as.data.frame(colData(mat))
meta0 <- meta%>%
    mutate(bti=paste(Cluster2, treat, sampleID, sep="_"))%>%
    rownames_to_column(var="NEW_BARCODE")%>%
    dplyr::select(NEW_BARCODE, bti)
 
## C5, C6 and C7 combine new cluster 9
## meta2 <- meta%>%dplyr::filter(Cluster2%in%c("C5", "C6", "C7"))%>%
##     mutate(Cluster2="C9", bti=paste(Cluster2, treat, sampleID, sep="_"))%>%
##     rownames_to_column(var="NEW_BARCODE")%>%
##     dplyr::select(NEW_BARCODE, bti)
    
meta <- meta0 ## rbind(meta0, meta2)    
 
dd <- meta%>%group_by(bti)%>%summarise(ncell=n(), .groups="drop")    

    
###
### features
## infn_feature <- "../1.2_ArchR_process/4_motif_output/Motif_jaspar2022.motifinfor.S4.rds"
## motif_DF <- as.data.frame(read_rds(infn_feature))

## motif_DF <- as.data.frame(motif_DF)
## motif_DF <- motif_DF%>%mutate(motif_name=gsub("\\.\\.", "::", gsub("_.*", "", name)))

## infn_feature2 <- "../1.2_ArchR_process/4_motif_output/Motif_jaspar2022.motifinfor.df.txt"
## motif_DF2 <- read.table(infn_feature2, header=T)
## motif_DF2 <- motif_DF2[,1:2]

    
###    
### Assays and obtain pseudobulk data

Y <- assays(mat)$z
Y <- Y[, meta$NEW_BARCODE]

identical(meta$NEW_BARCODE, colnames(Y))

### pseudobulk data     

X <- model.matrix(~0+meta$bti)
YtX <- Y%*%X

 
## rownames(YtX) <- peak
colnames(YtX) <- gsub(".*bti", "", colnames(YtX))    

   

###
### save data 
opfn <- paste(outdir, "0_ncell.rds", sep="")
write_rds(dd, file=opfn)

###
opfn2 <- paste(outdir, "1_YtX_comb.rds", sep="")
write_rds(YtX, file=opfn2)


###
### average cells 
dd_fn <- paste(outdir, "0_ncell.rds", sep="")
dd <- read_rds(dd_fn)

fn2 <- paste(outdir, "1_YtX_comb.rds", sep="")
YtX <- read_rds(fn2)


###
### number of cell and keep the name of ncell the same to the colnames of the YtX 
ncell <- dd$ncell
names(ncell) <- dd$bti

ncell <- ncell[colnames(YtX)]
 
identical(names(ncell), colnames(YtX))

###
YtX_ave <- sweep(YtX, 2, 1/ncell, "*")
opfn3 <- paste(outdir, "1_YtX_ave.rds", sep="")
write_rds(YtX_ave, file=opfn3)
 


###
### 




## ###
## ### summary 

## ncell_fn <- paste(outdir, "0_ncell.rds", sep="")
## dd <- read_rds(ncell_fn)

## x <- str_split(dd$bti, "_", simplify=T)

## dd <- dd%>%mutate(Cluster=factor(gsub("^C", "", x[,1])), treat=x[,2])

## ###
## ### plots
## p0 <- ggplot(dd, aes(x=ncell))+
##    geom_histogram(fill="white", color="grey50")+
##    xlab("#cells")+ylab("#combinations")+ 
##    facet_grid(treat~Cluster, scales="free")+
##    theme_bw()+
##    theme(axis.title=element_text(size=9),
##          axis.text=element_text(size=6))
         

## ### save
## figfn <- paste(outdir, "Figure0_ncell.hist.png", sep="")
## ggsave(figfn, p0, width=950, height=300, units="px", dpi=120)




##############################
### get clean data
##############################

###
ncell_fn <- paste(outdir, "0_ncell.rds", sep="")
dd <- read_rds(ncell_fn)

### peak count data 
fn2 <- paste(outdir, "1_YtX_ave.rds", sep="")
count <- read_rds(fn2)


### auto peaks
## EXPs <- read.table("HOLD.lib.txt")$V1
## peakfn <- paste("../1.2_ArchR_process/3_reCallPeaks_output/Peak_matrix/", EXPs[1], "_peakMatrix.rds", sep="")
## mat <- read_rds(peakfn)
## peak_df <- rowRanges(mat)
## anno <- data.frame(peak=paste(seqnames(peak_df), start(peak_df), end(peak_df), sep="_"),
##                    chr=as.character(gsub("^chr", "", seqnames(peak_df))))
## anno$rnz <- rowSums(count)
## auto <- as.character(1:22)
## peakSel <- anno%>%filter(rnz>0, chr%in%auto)%>%pull(peak)


##
for ( num in c(0, 20, 100)){
###    
bti2 <- dd%>%filter(ncell>num)%>%pull(bti)
count2 <- count[, bti2]

### save out
opfn <- paste(outdir, "1.3_YtX_ave.th", num, ".clean.rds", sep="")
write_rds(count2, file=opfn)

cat("Filter:", num, " motifs:", nrow(count2), " combs:", ncol(count2), "\n")

}

## fn <- "/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/analyses_2024-09-13/2_Differential_analysis/2_motif.outs/Cluster_res0.07/1.3_YtX_ave.th20.clean.rds"
## x <- read_rds(fn)

###
### END






