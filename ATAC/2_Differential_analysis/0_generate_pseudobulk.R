##
library(Matrix)
library(tidyverse)
library(SummarizedExperiment)
###

option <- "option_nFeature15K"
outdir <- paste("./0_pseudobulk.outs/", option, "_cluster_res0.12/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T) 



####
### generate pseudobulk data 

### get pseudobulk data for each library
### 0_generate_pseudobulk.R and 0_submit_perLib.sh
  
## ###
## is_same <- sapply(EXPs, function(ii){
## peakfn <- paste("../1.2_ArchR_process/3_reCallPeaks_output/Peak_matrix/", ii, "_peakMatrix.rds", sep="")
## mat <- read_rds(peakfn)
## peak_df <- rowRanges(mat)
## peak2 <- paste(as.character(seqnames(peak_df)), as.character(start(peak_df)), as.character(end(peak_df)), sep="_")

## cat(ii, "\n")
## identical(peak0, peak2)
## })    

## sum(is_same)
## length(EXPs)

###
### combine pseudobulk counts
YtX_all <- lapply(EXPs, function(ii){
   ##
   fn <- paste(outdir, "0_", ii, ".YtX_comb.rds", sep="")
   mat <- read_rds(fn)
   cat(ii, dim(mat), "\n") 
   mat
})
### save
YtX_all <- do.call(cbind, YtX_all)
opfn <- paste(outdir, "1_YtX_comb.rds", sep="")
write_rds(YtX_all, file=opfn)

###
### rbind data frame of ncells 
dd_all <- map_dfr(EXPs, function(ii){
    ##
    cat(ii, "\n")
    fn <- paste(outdir, "0_", ii, ".ncell.rds", sep="")
    df0 <- read_rds(fn)
    df0
})
### save
opfn <- paste(outdir, "1_dd.ncell.rds", sep="")
write_rds(dd_all, file=opfn)




###
### summary 
ncell_fn <- paste(outdir, "1_dd.ncell.rds", sep="")
dd <- read_rds(ncell_fn)

x <- str_split(dd$bti, "_", simplify=T)

dd <- dd%>%
    mutate(Cluster=x[,1], cl_val=as.numeric(as.character(gsub("^C", "", x[,1]))),
           treat=x[,2])%>%
    mutate(Cluster_sort=fct_reorder(Cluster, cl_val))


### 
### plots
p0 <- ggplot(dd, aes(x=ncell))+
   geom_histogram(fill="white", color="grey50")+
   xlab("#cells")+ylab("#combinations")+ 
   facet_grid(treat~Cluster_sort, scales="free")+
   theme_bw()+
   theme(axis.title=element_text(size=9),
         axis.text.x=element_text(size=6),
         axis.text.y=element_text(size=9))
         

### save
figfn <- paste(outdir, "Figure0_ncell.hist.pdf", sep="")
ggsave(figfn, p0, width=10, height=3)




##############################
### get clean data
##############################

option <- "option_nFeature15K"
outdir <- paste("./0_pseudobulk.outs/", option, "_cluster_res0.12/", sep="")
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)  

## dir_peak <- paste("../1.2_ArchR_process/3_reCallPeaks_output/", option, "/Peak_matrix_lib/", sep="")

###
### ncells per combination
ncell_fn <- paste(outdir, "1_dd.ncell.rds", sep="")
dd <- read_rds(ncell_fn)

###
### peak count data 
peak_fn <- paste(outdir, "1_YtX_comb.rds", sep="")
count <- read_rds(peak_fn)

### select peaks in auotosome 
anno <- data.frame(peak=rownames(count), chr=gsub("_.*", "", rownames(count)))
auto <- paste("chr", 1:22, sep="")
peakSel <- anno%>%filter(chr%in%auto)%>%pull(peak)


##
for ( num in c(0, 20, 100)){
###    
bti2 <- dd%>%filter(ncell>num)%>%pull(bti)
count2 <- count[peakSel, bti2]
rnz <- rowSums(count2)
count2 <- count2[rnz>0, ]    
### save out
opfn <- paste(outdir, "2_YtX_comb.th", num, ".clean.rds", sep="")
write_rds(count2, file=opfn)

cat("Filter", num, nrow(count2), "peaks", ncol(count2), "combs", "\n")

}



###
### summary number of individual for each cluster and treatment

fn <- paste(outdir, "2_YtX_comb.th20.clean.rds", sep="")
count <- read_rds(fn)
x <- as.data.frame(str_split(colnames(count), "_", simplify=TRUE)) 
names(x) <- c("Cluster", "treat", "sampleID")
x$bti <- colnames(count)

####
fn <- paste(outdir, "1_dd.ncell.rds", sep="")
dd <- read_rds(fn)


###
x2 <- x%>%left_join(dd, by="bti")
summ <- x2%>%group_by(Cluster, treat)%>%
    summarize(nind=n(), ncell_md=median(ncell), .groups="drop")%>%
    ungroup()


summ2 <- summ%>%mutate(cl_val=as.numeric(gsub("C", "", Cluster)))%>%arrange(cl_val) 
summ2 <- summ2%>%dplyr::select(-cl_val)
 
opfn <- paste(outdir, "2_summ.infor.tsv", sep="")
write_tsv(summ2, file=opfn)



###
###
fn <- paste(outdir, "2_summ.infor.tsv", sep="")
summ <- read_tsv(fn)

opfn <- paste("./Cluster/", option, "_cluster_res0.12.cl.txt", sep="")
write.table(unique(summ$Cluster), file=opfn, col.names=FALSE, row.names=FALSE, quote=FALSE)
###
### END



