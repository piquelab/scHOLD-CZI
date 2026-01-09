######################################
### demux2 for downstream analysis ###
######################################
library(Matrix)
library(tidyverse)
library(parallel)
library(data.table)



rm(list=ls())

outdir <- "./0_demux.outs/"
if (!file.exists(outdir)) dir.create(outdir, showWarnings=F)

###
### script for demux results
### filtering criteria
### SNG-955,996 cells
### samples in covariates (165)-937,717 cells
### correct batch_individual-798,363 cells 



############################
### 1, read demuxlet data ###
############################


### demux file
basefolder <- "/rs/rs_grp_schold/CZI/ATAC/counts_cellranger_atac/fastdemux/fdout/"
###
demuxfn <- list.files(basefolder, pattern="*.fdout.raw.info.txt.gz")
expNames <- gsub(".fdout.raw.info.txt.gz", "", demuxfn) 
names(demuxfn) <- expNames

##
demux_ls <- mclapply(expNames,function(ii){

  ### demux results  
  fn <- paste0(basefolder, demuxfn[ii], sep="")
  dd <- fread(fn, header=T, data.table=F)
  ###
    
  ## ### filtered barcodes    
  ## fn2 <- paste(base2folder, ii, "/filtered_peak_bc_matrix/barcodes.tsv", sep="")
  ## barSel <- fread(fn2, header=F, data.table=F)$V1  
  ###   
  ## dd2 <- dd%>%filter(BARCODE%in%barSel)%>%
  ##     mutate(NEW_BARCODE=paste(ii, "_", gsub("-1", "", BARCODE), sep=""), EXP=ii)
    
  dd2 <- dd%>%mutate(NEW_BARCODE=paste(ii, "_", gsub("-1", "", BARCODE), sep=""), EXP=ii)  
  ###  
  cat(ii, nrow(dd2), sum(dd2$DropType==1), "\n")  
  dd2 
},mc.cores=2)

# merging all the experiments  (31 experiments)
demux <- do.call(rbind, demux_ls)


## saving the final demux file obtained from all experiments into "1_demux_New.ALL.rds" file
### output
opfn <- paste(outdir, "1_demux_all.rds", sep="")
write_rds(demux, opfn)

###
### DROPLET.TYPE: "SNG" filtering  
opfn2 <- paste(outdir, "1_demux.SNG.rds", sep="")
demux2 <- demux %>% filter(DropType==1)
write_rds(demux2, file=opfn2)


## ########################################################
## ### 2. visulization data quality for filtered cells
## #########################################################


## outdir <- "./0_demux.outs/"

## fn <- paste(outdir, "1_demux_filtered.rds", sep="")
## dd <- read_rds(fn)%>%mutate(DropType2=ifelse(DropType>1, "DBL", "SNG"))



## ### droplet type 
## summ2 <- dd%>%
##    group_by(EXP, DropType2)%>%
##    summarise(ncell=n(),.groups="drop")%>%ungroup()

## x <- summ2%>%group_by(EXP)%>%mutate(percent=ncell/sum(ncell))%>%ungroup()
## x0 <- x%>%filter(DropType2=="DBL")

## summ2 <- summ2%>%
##     mutate(EXP_val=as.numeric(gsub("^HOLD|-ATAC-.*", "", EXP)),
##            EXP2=gsub("-ATAC", "", EXP),             
##            EXP_sort=fct_reorder(EXP2, EXP_val))
 

## p <- ggplot(summ2)+
##    geom_bar(stat="identity", position=position_fill(reverse=F),
##         aes(x=EXP_sort, y=ncell, fill=DropType2))+
##    theme_bw()+
##    theme(legend.title=element_blank(),
##          legend.key.size=grid::unit(0.4, "cm"),
##          axis.title=element_blank(),
##          axis.text.x=element_text(angle=45, hjust=1, size=7))

## figfn <- paste(outdir, "Figure0.1_droplet.png", sep="")
## ggsave(figfn, p, width=700, height=400, units="px", dpi=120)



## ### Best score distribution

## p2 <- ggplot(dd, aes(BestScore, fill=DropType2, color=DropType2))+
##     geom_density(alpha=0.1)+
##     xlab("Best score")+ylab("density")+
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           legend.key.size=grid::unit(0.4, "cm"))

## ###
## figfn <- paste(outdir, "Figure0.2_bestscore.density.png", sep="")
## ggsave(figfn, p2, width=450, height=380, units="px", dpi=120)
 


##############################
### 3. summary SNG cells  
##############################



###
### covariate file individual batch infor
fn <- "/rs/rs_grp_schold/covariates/dbgap/HOLD_library_metadata_prelim_n165_dbgapIDs_batch_09_17_2024.txt"
cv <- read.table(fn, header=T) ## 165 individuals
cv2 <- cv%>%
    dplyr::mutate(sampleID_batch=paste(dbgap.ID, Batch, sep="_"))%>%
    dplyr::select(sampleID=dbgap.ID, sampleID_batch)

###
### 183 individuals 
fn <- "./0_demux.outs/1_demux.SNG.rds"
dd <- read_rds(fn) ## 955,996 cells
## setdiff(unique(dd$BestSample), cv$dbgap.ID)
dd <- dd%>%filter(BestSample%in%cv2$sampleID) ##937,717 cells


summ <- dd%>%group_by(EXP, BestSample)%>%summarize(ncell=n(), .groups="drop")%>%ungroup()

summ2 <- summ%>%group_by(BestSample)%>%mutate(percent_ind=ncell/sum(ncell))%>%ungroup()
summ2 <- summ2%>%group_by(EXP)%>%mutate(percent_EXP=ncell/sum(ncell))%>%ungroup()
summ2 <- summ2%>%
    dplyr::rename("sampleID"="BestSample")%>%
    mutate(EXP2=gsub("-ATAC", "", EXP), Batch=gsub("-ATAC-.*", "", EXP),
           EXP_val=as.numeric(gsub("^HOLD|-ATAC.*", "", EXP)))
           


summ2 <- summ2%>%left_join(cv2, by="sampleID")



summ3 <- summ2%>%
    mutate(sampleID_val=as.numeric(gsub(".*_HOLD", "", sampleID_batch)))%>%
    mutate(EXP2_sort=fct_reorder(EXP2, EXP_val),
           sampleID_batch_sort=fct_reorder(sampleID_batch, sampleID_val))

summ3 <- as.data.frame(summ3)


###
### percent for each individual
p <- ggplot(summ3, aes(x=sampleID_batch_sort, y=EXP2_sort, fill=percent_ind))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=7),
         axis.title=element_blank())

###
figfn <- paste(outdir, "Figure1.1_ind.heatmap.png", sep="")
ggsave(figfn, p, width=800, height=350, units="px", dpi=100)

###
figfn <- paste(outdir, "Figure1.1_ind.heatmap.pdf", sep="")
ggsave(figfn, p, width=8, height=3.5)



###
### percent for each experiment 
p2 <- ggplot(summ3, aes(x=sampleID_batch_sort, y=EXP2_sort, fill=percent_EXP))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=7),
         axis.title=element_blank())

###
figfn <- paste(outdir, "Figure1.2_EXP.heatmap.png", sep="")
ggsave(figfn, p2, width=800, height=350, units="px", dpi=100)


figfn <- paste(outdir, "Figure1.2_EXP.heatmap.pdf", sep="")
ggsave(figfn, p2, width=8, height=3.5)



##############################################################################################
### 4. Summary clean data, focusing on filtered cells, SNG and  removing missmatch cells 
##############################################################################################



fn <- "/rs/rs_grp_schold/covariates/dbgap/HOLD_library_metadata_prelim_n165_dbgapIDs_batch_09_17_2024.txt"
cv <- read.table(fn, header=T) ## 165 individuals

## fn <- "/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_HOLD01-HOLD14_dbgap.ID_cziexp_fixed_05_03_2024.txt"
## cv <- read.table(fn, header=T)
cv2 <- cv%>%
    dplyr::mutate(sampleID_batch=paste(dbgap.ID, Batch, sep="_"))%>%
    dplyr::select(sampleID=dbgap.ID, sampleID_batch)


### 
fn <- "./0_demux.outs/1_demux.SNG.rds"
dd <- read_rds(fn)
## setdiff(unique(dd$BestSample), cv$dbgap.ID)
dd <- dd%>%filter(BestSample%in%cv2$sampleID)%>%
   mutate(Batch=gsub("-ATAC-.*", "", EXP),
          id=gsub("^HOLD", "", Batch),
          id2=ifelse(nchar(id)==1, paste0("0", id), id), Batch2=paste0("HOLD", id2),
          sampleID_batch=paste(BestSample, Batch2, sep="_"))



dd2 <- dd%>%filter(sampleID_batch%in%cv2$sampleID_batch)  ### 798,363 cells

### 462,673 cells
opfn <- paste(outdir, "2_demux.SNG.correctBatch.rds", sep="")
write_rds(dd2, file=opfn)

print(object.size(dd2), units="GB")


## x  <- read_rds("./torm/0_demux.outs/2_demux.SNG.correctBatch.rds")  ## 691,125 cells for old demux 





##################################################
### visulization after filtering mismatching  
####################################################

fn <- "./0_demux.outs/2_demux.SNG.correctBatch.rds"
dd <- read_rds(fn)   #### 798,363 cells


summ <- dd%>%group_by(EXP, sampleID_batch)%>%summarize(ncell=n(), .groups="drop")%>%ungroup()

summ2 <- summ%>%group_by(sampleID_batch)%>%mutate(percent_ind=ncell/sum(ncell))%>%ungroup()
summ2 <- summ2%>%group_by(EXP)%>%mutate(percent_EXP=ncell/sum(ncell))%>%ungroup()
summ2 <- summ2%>%
    mutate(EXP2=gsub("-ATAC", "", EXP),
           EXP_val=as.numeric(gsub("^HOLD|-ATAC.*", "", EXP)),
           Batch=gsub(".*_", "", sampleID_batch),
           Batch_val=as.numeric(gsub("^HOLD", "", Batch)))
           

summ3 <- summ2%>%
    mutate(EXP2_sort=fct_reorder(EXP2, EXP_val),
           sampleID_batch_sort=fct_reorder(sampleID_batch, Batch_val))%>%
    as.data.frame()


###
### percent for each individual
p <- ggplot(summ3, aes(x=sampleID_batch_sort, y=EXP2_sort, fill=percent_ind))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=7),
         axis.title=element_blank())

### png
figfn <- paste(outdir, "Figure2.1_ind.heatmap.png", sep="")
ggsave(figfn, p, width=800, height=350, units="px", dpi=100)

### pdf
figfn <- paste(outdir, "Figure2.1_ind.heatmap.pdf", sep="")
ggsave(figfn, p, width=8, height=3.5)



###
### percent for each experiment 
p2 <- ggplot(summ3, aes(x=sampleID_batch_sort, y=EXP2_sort, fill=percent_EXP))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=7),
         axis.title=element_blank())

###
figfn <- paste(outdir, "Figure2.2_EXP.heatmap.png", sep="")
ggsave(figfn, p2, width=800, height=350, units="px", dpi=100)
 
### pdf
figfn <- paste(outdir, "Figure2.2_EXP.heatmap.pdf", sep="")
ggsave(figfn, p2, width=8, height=3.5)

###
### END

###
###
## fn <- "./0_demux.outs/2_demux_filtered.SNG.correctBatch.rds"
## dd <- read_rds(fn)


## summ <- dd%>%group_by(EXP, sampleID_batch)%>%summarize(ncell=n(), .groups="drop")%>%ungroup()
## summ <- summ%>%mutate(treat=gsub(".*-ATAC-", "", EXP), Batch=gsub("-ATAC-.*", "", EXP))


## b1 <- summ%>%filter(treat=="CTRL")%>%pull(Batch)%>%unique()%>%sort()
## b2 <- summ%>%filter(treat=="LPS")%>%pull(Batch)%>%unique()%>%sort() 
## b3 <- summ%>%filter(treat=="LPS-DEX")%>%pull(Batch)%>%unique()%>%sort()

###
### END




###############################
### summary demulet results ###
###############################

## demux <- read_rds("./1_demux_output/1_demux_New.SNG.rds")
## aa <- demux %>% dplyr::filter(NUM.READS>20,NUM.SNPS>20,NUM.READS<20000) %>%
##     select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,Sample_ID=SNG.BEST.GUESS)

## ##
## aa2 <- aa %>% group_by(EXP,Sample_ID) %>%
##     summarize(n=n(), .groups="drop")%>%ungroup()%>%
##     pivot_wider(id_cols=EXP, names_from=Sample_ID, values_from=n, values_fill=0)

## ###
## opfn <- paste(outdir, "cell.counts.matrix.tsv", sep="")
## write_tsv(aa2, opfn)



## aa3 <- aa %>% group_by(EXP,Sample_ID) %>%
##     summarize(n=n(), .groups="drop")%>%ungroup()
## opfn <- paste(outdir, "cell.counts.tsv", sep="")
## write_tsv(aa3, opfn)












#
