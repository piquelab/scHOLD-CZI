#
library(tidyverse)
library(parallel)
library(viridis)
library(pheatmap)

##library(data.table)

### 
args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_HOLD01-HOLD14_dbgap.ID_cziexp_02_16_2024.txt","CZ1_group.txt") #for testing
outFolder=paste0(args[1],"1_demux_output/")
cov_file=args[2]
# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"umatched_removed/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#read in samples file (just list of samples to run, each sample on newline)
if(!is.na(args[3])){
samples=read.table(args[3],header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
}

opfn <- paste0(outFolder,"1_demux_New.SNG.rds")
demux <- read_rds(opfn)

opfn <- paste0(outFolder,"1_demux_New.ALL.rds")
tdemux <- read_rds(opfn)

###
demux <- demux%>%
         mutate(
                #BEST.GUESS=gsub(",.*", "", BEST.GUESS),
                #NEXT.GUESS=gsub(",.*", "", NEXT.GUESS),
                #BATCH=substring(EXP,1,6),
                BATCH=gsub("-.*", "", EXP),
                treats=gsub(".*[0-9].{,2}-","",EXP),
                #chemi=grepl("SCAIP5V3|SCAIP6",BATCH2),
                #chem=ifelse(chemi,"V3", "V2")
                ) 

#head(demux)

##This step not necessary for current data as only control and LPS - may need to uncomment
# correct the EtOH ETOH discrepency in treat column 
#demux$treats <- gsub("EtOH", "ETOH", demux$treats)

# filter 
aa <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID=SNG.BEST.GUESS) 

# check numbers when grouped by SampleID, treats, Batches, separatly
xx <- aa %>% group_by(treats) %>% summarize(n=n())
cat("numbers grouped by treat= ")
print(xx)
xx <- aa %>% group_by(BATCH) %>% summarize(n=n())
cat("numbers grouped by batch= ")
print(xx)
xx <- aa %>% group_by(Sample_ID) %>% summarize(n=n())
cat("numbers grouped by sample= ")
print(xx)

xx <- xx %>% filter(n>100) 
cat("sum sample numbers after 100 filter= ",sum(xx$n))
cat("unique sample numbers after 100 filter= ",length(unique(xx$Sample_ID)))

cell.counts <- aa %>% group_by(EXP,Sample_ID, BATCH, treats) %>% summarize(n=n()) 
cat("sum cell counts= ",sum(cell.counts$n))
cat("unique cell counts samples= ",length(unique(cell.counts$Sample_ID)))

cell.counts <- cell.counts%>%group_by(Sample_ID)%>%
   mutate(Perc=n/sum(n), comb=paste(EXP, Sample_ID, sep="_"))%>%as.data.frame()

cell.counts.filt <- cell.counts %>% filter(n>100) 
n100toremove <- cell.counts %>% filter(!n>100) 

cat("sum cell counts filtered n>100= ",sum(cell.counts.filt$n))
cat("unique cell counts filtered n>100 samples= ",length(unique(cell.counts.filt$Sample_ID)))

cell.counts.filt <- cell.counts.filt %>% mutate(ID_BATCH=paste(Sample_ID,BATCH, sep="_"))

write_tsv(cell.counts.filt,paste0(outFolder,"cell.count.filter.100min.v1.tsv"))


##################################################################
############ experimental data sheet to find unmatched ###########
##################################################################


# load experimental data sheet
#covariate file numbers are in format 01,02,etc which does not match file naming labels of 1,2,etc
exp <- read.table(cov_file, row.names=NULL,header=T)
exp$Batch <- gsub("HOLD0","HOLD",exp$Batch)
if(!is.na(args[3])){
  exp <- exp %>% filter(Batch %in% samples$Batch)
}
dim(exp)
table(exp$Batch)

# find what IDs are missing from the filtered data
missing <- setdiff(exp$dbgap.ID, cell.counts.filt$Sample_ID)
IDmiss <- exp %>% filter(dbgap.ID %in% missing) %>% mutate(Sample_ID=dbgap.ID) %>% select(dbgap.ID, Batch, Sample_ID) 
print(IDmiss)
write.csv(IDmiss, paste0(outFolder,"missing_samples_in_demultiplexing_results.csv"), row.names=F)

# merge experimental data and filtered data
expsub <- exp %>% mutate(Sample_ID=dbgap.ID) %>% select(dbgap.ID, Batch, Sample_ID) #removed Participant.ID as not a variable
merge <- left_join(cell.counts.filt, expsub, by = "Sample_ID")

#there were a bunch of summary items here, kept only QC
cat("checking NA ids and mismatches of ids")
sum(is.na(merge$Sample_ID))
sum(is.na(merge$dbgap.ID))
sum(merge$Sample_ID != merge$dbgap.ID)
sum(merge$Batch != merge$BATCH)

# samples where IDs dont match
unmatch_ID <- merge %>%
  filter(Sample_ID != dbgap.ID)
print(unmatch_ID)

# samples where ID_batch dont m

# samples where demultiplexing BATCH does not match to the experimental batch for a sampleIDs
unmatch_batch <- merge %>% filter(BATCH != Batch) # -- at least in current data seems to only be instances of NA in exp batch
merge <-  transform(merge, ID_Batch= paste0(dbgap.ID, "_", Batch))# , sep="_" ))
unmatch_comb <- subset(merge, ID_BATCH != ID_Batch)
print(unmatch_comb)

write.csv(unmatch_comb, paste0(outFolder,"umatched_batches_samples_IDs_compared_to_exp_file.csv"), row.names=F)

### remove the unmatched from the >100 filtered data. 
dim(cell.counts.filt)
unmatch_comb$comb

# how many cells are unmatched
cell.unm <- cell.counts.filt %>% filter(comb %in% unmatch_comb$comb)
dim(cell.unm)
sum(cell.unm$n)

# how many cells remain
rem.unm <- cell.counts.filt %>% filter(!comb %in% unmatch_comb$comb)
dim(rem.unm)
sum(rem.unm$n)

write_tsv(rem.unm,paste0(outFolder,"cell.count.filter.100min.v2_unmatched_cells_removed.tsv"))


################################
############ Figures ###########
################################

cell.counts.filt <- rem.unm

fig0 <- ggplot(cell.counts.filt,aes(x=EXP, y=n, fill=factor(BATCH)))+
        geom_bar(stat="identity", position = position_dodge())+ #, position = position_dodge())+
        ylab("Number Barcodes")+xlab("")+
        #geom_text(aes(label=n),color="black",position = position_dodge(0.9), vjust = 0, size=3)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=15),
              axis.text.y=element_text(hjust=1, size=15),
              axis.title.y=element_text(size=20),
              legend.text=element_text(size=15))  
figfn <- paste(figuredir, "Figure01_barcodes.png", sep="")
png(figfn, width=8500, height=4000, res=380)
fig0
dev.off()


fig0 <- cell.counts.filt %>% ggplot(aes(y = Sample_ID, x = EXP, fill = n)) + # can also try x=BATCH
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") 
figfn <- paste(figuredir, "Figure02_tileplot_100min.png", sep="")
png(figfn, width=1300, height=1000, res=120)
fig0
dev.off()
## should make hierarchical, so all the samples in one batch/exp are plotted together
top10m <- reshape2::dcast(EXP  ~ Sample_ID, data=cell.counts.filt, value.var="n",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$EXP  
top10m[is.na(top10m)] <- 0
breaksList = c(0,seq(0.1, 4, by = 0.5))

figfn <- paste(figuredir, "Figure02_heatmap_100min.png", sep="")
png(figfn, width=1300, height=1000, res=120)
pheatmap(top10m[,-1], color = c("grey",inferno(length(breaksList))[-1]), breaks = breaksList, scale="row", cluster_cols = T, cluster_rows = FALSE, clustering_distance_rows = 'correlation',show_colnames = T)
dev.off()

### Heat map of fractions per individuals in each library
p <- ggplot(cell.counts.filt, aes(x=comb, y=EXP, fill=Perc))+
   geom_tile()+
   scale_fill_gradient("Fraction of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())
###
figfn <- paste(figuredir, "Figure03.2_heatmap.png", sep="")
png(figfn, width=1300, height=1000, res=120)
print(p)
dev.off()


### Heat map of number of cells per individuals in each library
p <- ggplot(cell.counts.filt, aes(x=comb, y=EXP, fill=n))+
   geom_tile()+
   scale_fill_gradient("Number of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())
###
figfn <- paste(figuredir, "Figure03.3_heatmap_n.png", sep="")
png(figfn, width=1300, height=1000, res=120)
print(p)
dev.off()


### Heat map of number of cells per individuals in each library
p <- ggplot(cell.counts.filt, aes(x=Sample_ID, y=EXP, fill=n))+
   geom_tile()+
   scale_fill_gradient("Number of cells",
      low="#ffffc8", high="#7d0025", na.value=NA)+
   theme_bw()+
   theme(axis.text.x=element_blank(), #element_text(hjust=1, vjust=0.5, angle=90, size=6),
         axis.text.y=element_text(size=8),
         axis.title=element_blank())
###
figfn <- paste(figuredir, "Figure03.4_heatmap_n.png", sep="")
png(figfn, width=1300, height=1000, res=120)
print(p)
dev.off()

# heatmap
cellmatexp <- cell.counts.filt %>% ungroup %>% select(EXP,Sample_ID,n) %>% 
  pivot_wider(names_from=EXP,values_from=n,values_fill=0) %>% 
  column_to_rownames("Sample_ID") %>%
  as.matrix()

library(pheatmap)

figfn <- paste(figuredir, "Figure03_heatmap_100min.png", sep="")
png(figfn, width=1000, height=1800, res=120)
#pdf("heatmap.pdf",height=24,width=10)
pheatmap(cellmatexp)
dev.off()

figfn <- paste(figuredir, "Figure03_heatmap_100min_exp_on_y.png", sep="")
png(figfn, width=2200, height=1000, res=120)
#pdf("heatmap.pdf",height=24,width=10)
pheatmap(t(cellmatexp), 
          cluster_rows = TRUE, 
          cluster_cols = TRUE)
dev.off()


####### load aa again, remove the umatched cells, and the <100. then group to plot per exp and batch
aa <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID=SNG.BEST.GUESS) 

bb <- aa %>% mutate(comb=paste(EXP, Sample_ID, sep="_"))
bb <- bb %>% filter(!comb %in% unmatch_comb$comb) %>% filter(!comb %in% n100toremove$comb)

dd <- bb %>% group_by(EXP) %>% summarize(n=n()) 
sum(dd$n)

# make a batch column
dd <- dd %>% separate(EXP,c('BATCH','x', 'y'),'-',remove=FALSE)  

fig0 <- ggplot(dd,aes(x=EXP, y=n, fill=factor(BATCH)))+
        geom_bar(stat="identity", position = position_dodge())+ #, position = position_dodge())+
        ylab("Number Barcodes")+xlab("")+
        geom_text(aes(label=n),color="black",position = position_dodge(0.9), vjust = 0, size=4)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=15),
              axis.text.y=element_text(hjust=1, size=15),
              axis.title.y=element_text(size=20),
              legend.text=element_text(size=15))  
figfn <- paste(figuredir, "Figure01.2_barcodes.png", sep="")
png(figfn, width=8500, height=4000, res=380)
fig0
dev.off()


################################

#removed section here that did not seem to apply to current data
#files not made: cc.wider.v1.tsv, cell.counts.tsv, cc.wider.tsv
