##
library(Matrix)
library(tidyverse)
library(SummarizedExperiment)
library(argparse)
rm(list=ls())

time0 <- Sys.time()

###
### passing arguments 
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                        help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
                        dest="verbose", help="Print little output")
parser$add_argument("-cond", "--condition", type="character", default="HOLD1-ATAC-CTRL",
                        help="Library, cluster or cell type",
                        metavar="character")
parser$add_argument("-opt", "--option", type="character", default="option_nFeature15k",
                        help="Different approach to define cluster",
                        metavar="character")

###
### passing arguments
args <- parser$parse_args()
lib_ii <- args$condition
option <- args$option


###
### input and output path
## option <- "option_nFeature15K"
outdir <- paste("./0_pseudobulk.outs/", option, "_cluster_res0.12/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T) 

dir_peak <- paste("../1.2_ArchR_process/3_reCallPeaks_output/", option, "/Peak_matrix_lib/", sep="")



###
### check if row names is the same across library 
EXPs <- read.table("HOLD.lib.txt")$V1

peakfn <- paste(dir_peak, EXPs[1], "_peakMatrix.rds", sep="")
mat <- read_rds(peakfn)
peak_df <- rowRanges(mat)
peak0 <- paste(as.character(seqnames(peak_df)), as.character(start(peak_df)), as.character(end(peak_df)), sep="_")



###
### read data    
peakfn <- paste(dir_peak, lib_ii, "_peakMatrix.rds", sep="")
mat <- read_rds(peakfn)

    
###
### column
meta <- as.data.frame(colData(mat))
meta <- meta%>%
    mutate(bti=paste(Cluster2, treat, sampleID, sep="_"))%>%
    rownames_to_column(var="NEW_BARCODE")%>%
    dplyr::select(NEW_BARCODE, bti)

dd <- meta%>%group_by(bti)%>%summarise(ncell=n(), .groups="drop")    

    
###
### features
peak_df <- rowRanges(mat)
peak <- paste(as.character(seqnames(peak_df)), as.character(start(peak_df)), as.character(end(peak_df)), sep="_")

cat(lib_ii, "the same peak name", identical(peak, peak0), "\n")



###    
### Assays and obtain pseudobulk data
Y <- assays(mat)$PeakMatrix

cat(lib_ii, "#Cells", ncol(Y), "#Peaks", nrow(Y), "\n")

### pseudobulk data     
bti <- factor(meta$bti)
X <- model.matrix(~0+bti)
YtX <- Y%*%X

rownames(YtX) <- peak
colnames(YtX) <- gsub("^bti", "", colnames(YtX))    


###
### save data 
opfn <- paste(outdir, "0_", lib_ii, ".ncell.rds", sep="")
write_rds(dd, file=opfn)

###
opfn2 <- paste(outdir, "0_", lib_ii, ".YtX_comb.rds", sep="")
write_rds(YtX, file=opfn2)


###
### END

 





