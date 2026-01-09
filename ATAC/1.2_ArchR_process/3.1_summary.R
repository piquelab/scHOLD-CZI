##
library(Matrix)
library(tidyverse)
library(data.table)
library(SummarizedExperiment)
library(openxlsx)

###

option <- "option_nFeature15K"
outdir <- paste("./3_reCallPeaks_output/", option, "/Peak_matrix_cluster/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=F) 



#####################################################################
### get peak infor, what proportion of cells peaks are measured
#####################################################################

conditions <- read.table("cluster.txt")$V1

for ( ii in conditions){    
#### check if row names is the same across library 
peakfn <- paste(outdir, ii, "_peakMatrix.rds", sep="")
mat <- read_rds(peakfn)
peak_df <- rowRanges(mat)
peak0 <- paste(as.character(seqnames(peak_df)), as.character(start(peak_df)), as.character(end(peak_df)), sep="_")
    
Y <- assays(mat)$PeakMatrix

cat(ii, "peaks", nrow(Y), "cells", ncol(Y), "\n")    
    
df0 <- data.frame(Cluster = ii, peak_name = peak0, chr = as.character(seqnames(peak_df)), 
      prop_greater0 = rowMeans(Y>0), prop_greater1 = rowMeans(Y>1),
      prop_greater2 = rowMeans(Y>2), prop_greater3 = rowMeans(Y>3)) 

opfn <- paste(outdir, "summary_", ii, ".peak.txt.gz", sep="")
fwrite(df0, file = opfn, quote = F, sep = "\t")

}


 
### 5 no cells greater than 5




###
### summary  
conditions <- read.table("cluster.txt")$V1
autosome <- paste("chr", 1:22, sep = "")

for ( ii in conditions){
##ii <- conditions[1]
fn_summ  <- paste(outdir, "summary_", ii, ".peak.txt.gz", sep = "")
x <- fread(fn_summ, sep = "\t", data.table = F)
x <- x%>%filter(chr%in%autosome)

## rn <- names(x)
## rn2 <- rn[grepl("prop", rn)]

## x2 <- x%>%dplyr::select(all_of(rn2))
## colSums(x2>0)

dd <- x%>%
    pivot_longer(cols=starts_with("prop"), names_to = "thresh", names_prefix="prop_", values_to = "prop")

###
### plots
p0 <- ggplot(dd, aes(x=prop))+
   geom_histogram(fill="white", color="grey50")+
   xlab("proportion of cells")+ylab("#peaks")+
   ggtitle(ii)+ 
   facet_wrap(~thresh, ncol = 4, scales="free")+
   theme_bw()+
   theme(axis.title = element_text(size=9),
         axis.text = element_text(size=9),
         plot.title = element_text(size=12, hjust=0.5))
         

### save
figfn <- paste(outdir, "Figure_", ii, ".peaks.hist.pdf", sep="")
ggsave(figfn, p0, width=9.50, height=3)
}


###
###
summ_df <- NULL
for ( ii in conditions){
##ii <- conditions[1]
fn_summ  <- paste(outdir, "summary_", ii, ".peak.txt.gz", sep = "")
x <- fread(fn_summ, sep = "\t", data.table = F)
x <- x%>%filter(chr%in%autosome)

rn <- names(x)
rn2 <- rn[grepl("prop", rn)]

x2 <- x%>%dplyr::select(all_of(rn2))
summ0 <- apply(x2, 2, median)    
summ_df <- rbind(summ_df, summ0)
##
}    

summ_df <- as.data.frame(summ_df)
summ_df$cluster <- conditions

opfn <- paste(outdir, "summary_prop.median.xlsx", sep="")
write.xlsx(summ_df, file=opfn)


###
### calculate the proportion of peaks  greater than 0 across cells
### select the peaks with the proportion greater 0.01
summ_df <- NULL
th0 <- 0.01
for ( ii in conditions){
##ii <- conditions[1]
fn_summ  <- paste(outdir, "summary_", ii, ".peak.txt.gz", sep = "")
x <- fread(fn_summ, sep = "\t", data.table = F)
x <- x%>%filter(chr%in%autosome)

###    
rn <- names(x)
rn2 <- "prop_greater0"  ##rn[grepl("prop", rn)]
peak_val <- x%>%pull(rn2)
df0 <- data.frame("cluster"=ii, prop=mean(peak_val>th0), npeaks=sum(peak_val>th0))

summ_df <- rbind(summ_df, df0)
##
}    


opfn <- paste(outdir, "summary_prop.th0.01.xlsx", sep="")
write.xlsx(summ_df, file=opfn)
   

#####################################################
### get the cell type peaks used for analysis 
####################################################

conditions <- read.table("cluster.txt")$V1
autosome <- paste("chr", 1:22, sep = "")
th0 <- 0.01

###
###
peak_df <- map_dfr(conditions, function(ii){
   ##ii
   fn_summ  <- paste(outdir, "summary_", ii, ".peak.txt.gz", sep = "")
   x <- fread(fn_summ, sep = "\t", data.table = F)
   x2 <- x%>%filter(chr%in%autosome, prop_greater0 > th0)
   ### 
   cat(ii, nrow(x), sum(x$chr%in%autosome), nrow(x2), "\n")        
   x2 
})


###
###
opfn <- paste(outdir, "1_celltype_peak.txt.gz", sep="")
fwrite(peak_df, file = opfn, sep = "\t", quote=F)


###
### cell type peaks
fn <- paste(outdir, "1_celltype_peak.txt.gz", sep="")
x <- fread(fn, header=T)

summ <- x%>%group_by(Cluster)%>%summarise(npeak=n(), .groups = "drop")%>%ungroup()

summ2 <- summ%>%mutate(cl_val=as.numeric(gsub("C", "", Cluster)))%>%arrange(cl_val)

###
### save
opfn <- paste(outdir, "1_summ.xlsx", sep="")
write.xlsx(summ, file = opfn)



##############################
### get clean data
##############################

###
## ncell_fn <- paste(outdir, "0_ncell.rds", sep="")
## dd <- read_rds(ncell_fn)

## ### peak count data 
## peak_fn <- paste(outdir, "1_YtX_comb.rds", sep="")
## count <- read_rds(peak_fn)


## ### auto peaks
## EXPs <- read.table("HOLD.lib.txt")$V1
## peakfn <- paste("../1.2_ArchR_process/3_reCallPeaks_output/Peak_matrix/", EXPs[1], "_peakMatrix.rds", sep="")
## mat <- read_rds(peakfn)
## peak_df <- rowRanges(mat)
## anno <- data.frame(peak=paste(seqnames(peak_df), start(peak_df), end(peak_df), sep="_"),
##                    chr=as.character(gsub("^chr", "", seqnames(peak_df))))
## anno$rnz <- rowSums(count)
## auto <- as.character(1:22)
## peakSel <- anno%>%filter(rnz>0, chr%in%auto)%>%pull(peak)


## ##
## for ( num in c(0, 20, 100)){
## ###    
## bti2 <- dd%>%filter(ncell>num)%>%pull(bti)
## count2 <- count[peakSel, bti2]

## ### save out
## opfn <- paste(outdir, "1.2_YtX_comb.th", num, ".clean.rds", sep="")
## write_rds(count2, file=opfn)

## cat("Filter", num, nrow(count2), "peaks", ncol(count2), "combs", "\n")

## }

###
### END


## num <- 20
## fn <- paste(outdir, "1.2_YtX_comb.th", num, ".clean.rds", sep="")
## x <- read_rds(file=fn)

## fn <- paste(outdir, "0_ncell.rds", sep="")
## cvt <- read_rds(fn)%>%filter(ncell>20)

## x <- str_split(cvt$bti, "_", simplify=T)

## cvt$cluster <- x[,1]
## cvt$condition <- x[,2]
## cvt$sampleID <- x[,3]


## fn <- "/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/HOLD-CZI_geno_pc_nokin.txt"
## pca <- read.table(fn, header=T)
## pca <- pca%>%filter(Batch!="HOLD7")


## cvt2 <- cvt%>%filter(sampleID%in%pca$Sample_ID, condition=="CTRL")
 





