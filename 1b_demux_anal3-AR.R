#
library(tidyverse)
library(parallel)
##library(data.table)

### 
outFolder="/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/1_demux_output/"

opfn <- "./1_demux_output/1_demux_New.SNG.rds"
demux <- read_rds(opfn)

opfn <- "./1_demux_output/1_demux_New.ALL.rds"
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

head(demux)

# correct the EtOH ETOH discrepency in treat column 
demux$treats <- gsub("EtOH", "ETOH", demux$treats)

# filter 
aa <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID=SNG.BEST.GUESS) 

# check numbers when grouped by SampleID, treats, Batches, separatly
xx <- aa %>% group_by(Sample_ID) %>% summarize(n=n())
xx <- aa %>% group_by(treats) %>% summarize(n=n())
xx <- aa %>% group_by(Batches) %>% summarize(n=n())

xx <- xx %>% filter(n>100) 
sum(xx$n)
length(unique(xx$Sample_ID)) 


cell.counts <- aa %>% group_by(EXP,Sample_ID, BATCH, treats) %>% summarize(n=n()) 
sum(cell.counts$n)
length(unique(cell.counts$Sample_ID)) # all 251 ALOFT samples are present

cell.counts <- cell.counts%>%group_by(Sample_ID)%>%
   mutate(Perc=n/sum(n), comb=paste(EXP, Sample_ID, sep="_"))%>%as.data.frame()

 #mutate(comb=paste(EXP, sampleID, sep="_"))

cell.counts.filt <- cell.counts %>% filter(n>100) 
#cell.counts.filt <- cell.counts %>% filter(n>500) 
n100toremove <- cell.counts %>% filter(!n>100) 

sum(cell.counts.filt$n)
length(unique(cell.counts.filt$Sample_ID)) # 135, with n>100 and 131 with n>500


cell.counts.filt <- cell.counts.filt %>% mutate(ID_BATCH=paste(Sample_ID,BATCH, sep="_"))

write_tsv(cell.counts.filt,paste0(outFolder,"cell.count.filter.100min.v1.tsv"))


##################################################################
############ experimental data sheet to find unmatched ###########
##################################################################


# load experimental data sheet
exp <- read.csv("/wsu/home/groups/piquelab/SCAIP_2022/Ali/20-scALOFT/scALOFT2_Batches_Final_11-21-22 - scALOFT2_Batches_Final_11-21-22.csv", row.names=NULL)
exp <- exp %>% filter(!Batch %in% c('SCAIP19', 'SCAIP20'))
dim(exp)
table(exp$Batch)

# find what IDs are missing from the filtered data
missing <- setdiff(exp$dbgap.ID, cell.counts.filt$Sample_ID)
IDmiss <- exp %>% filter(dbgap.ID %in% missing) %>% mutate(Sample_ID=dbgap.ID) %>% select(dbgap.ID, Batch, Sample_ID) 
write.csv(IDmiss, "1_demux_output/missing_samples_in_demultiplexing_results.csv", row.names=F)

# merge experimental data and filtered data
expsub <- exp %>% mutate(Sample_ID=dbgap.ID) %>% select(dbgap.ID, Batch, Sample_ID, Participant.ID)
merge <- left_join(cell.counts.filt, expsub, by = "Sample_ID")

#summary(merge$Sample_ID)
#summary(merge$dbgap.ID)
#sum(is.na(merge$Sample_ID))
#sum(is.na(merge$dbgap.ID))
#sum(merge$Sample_ID != merge$dbgap.ID)
#sum(merge$Batch != merge$BATCH)

# samples where demultiplexing BATCH does not match to the experimental batch for a sampleIDs
unmatch_batch <- merge %>%
  filter(BATCH != Batch)
unmatch_batch

# samples where IDs dont match
unmatch_ID <- merge %>%
  filter(Sample_ID != dbgap.ID)
unmatch_ID

# samples where ID_batch dont match 
merge <-  merge %>% mutate(ID_Batch= paste0(dbgap.ID, "_", Batch))# , sep="_" ))
unmatch_comb <- merge %>%
  filter(ID_BATCH != ID_Batch)
unmatch_comb

write.csv(unmatch_comb, "1_demux_output/umatched_batches_samples_IDs_compared_to_exp_file.csv", row.names=F)


### remove the unmatched from the >100 filtered data. 13 unmatched combinations
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

# set new output dir for filtered out unmatched figures
outdir="/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/1_demux_output/umatched_removed/"


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
figfn <- paste(outdir, "Figure01_barcodes.png", sep="")
png(figfn, width=8500, height=4000, res=380)
fig0
dev.off()


fig0 <- cell.counts.filt %>% ggplot(aes(y = Sample_ID, x = EXP, fill = n)) + # can also try x=BATCH
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") 
figfn <- paste(outdir, "Figure02_tileplot_100min.png", sep="")
png(figfn, width=1300, height=1000, res=120)
fig0
dev.off()
## should make hierarchical, so all the samples in one batch/exp are plotted together



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
figfn <- paste(outdir, "Figure03.2_heatmap.png", sep="")
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
figfn <- paste(outdir, "Figure03.3_heatmap_n.png", sep="")
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
figfn <- paste(outdir, "Figure03.4_heatmap_n.png", sep="")
png(figfn, width=1300, height=1000, res=120)
print(p)
dev.off()



# heatmap
cellmatexp <- cell.counts.filt %>% ungroup %>% select(EXP,Sample_ID,n) %>% 
  pivot_wider(names_from=EXP,values_from=n,values_fill=0) %>% 
  column_to_rownames("Sample_ID") %>%
  as.matrix()

library(pheatmap)

figfn <- paste(outdir, "Figure03_heatmap_100min.png", sep="")
png(figfn, width=1000, height=1800, res=120)
#pdf("heatmap.pdf",height=24,width=10)
pheatmap(cellmatexp)
dev.off()

figfn <- paste(outdir, "Figure03_heatmap_100min_exp_on_y.png", sep="")
png(figfn, width=2200, height=1000, res=120)
#pdf("heatmap.pdf",height=24,width=10)
pheatmap(t(cellmatexp), 
          cluster_rows = TRUE, 
          cluster_cols = TRUE)
dev.off()



## figure out whats goin on with 14-PHA-DEX
cell.counts.filt %>% filter(BATCH=="SCAIP14")

# demuxlet did not run. gave an error. look at the slurm.out file
# /nfs/rprdata/SCAIP-ALOFT2/counts_cellranger_2023-05-27/demuxlet/slurm.SCAIP14-PHA-DEX.out

#ddx <- aa%>%group_by(Sample_ID,BATCH)%>%summarise(ncell=n())%>%arrange(BATCH)
#dd0 <- ddx[-17,]
#ddx <- aa%>%group_by(Sample_ID)%>%summarise(ncell=n())



####### load aa again, remove the umatched cells, and the <100. then group to plot per exp and batch
aa <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID=SNG.BEST.GUESS) 

# remove unmatched and remove <100
bb <- aa %>% mutate(comb=paste(EXP, Sample_ID, sep="_"))
head(n100toremove)
head(rem.unm)

bb <- bb %>% filter(!comb %in% unmatch_comb$comb) %>% filter(!comb %in% n100toremove$comb)

dd <- bb %>% group_by(EXP) %>% summarize(n=n()) 
sum(dd$n)
#length(unique(cell.counts$Sample_ID)) # all 251 ALOFT samples are present

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
figfn <- paste(outdir, "Figure01.2_barcodes.png", sep="")
png(figfn, width=8500, height=4000, res=380)
fig0
dev.off()





################################

cc.wider <- cell.counts %>% ungroup() %>% select(EXP,Pregnancy_ID,Origin2,n,Sample_ID) %>%
  filter(!is.na(Origin2)) %>% 
  group_by(EXP,Pregnancy_ID,Origin2) %>%  
  pivot_wider(names_from=Origin2,values_from=c(n,Sample_ID),values_fill=list(n=0,Sample_ID="UNK")) %>%
  mutate(n_T=n_M+n_F) %>% 
  filter(n_T>500)  %>%
  arrange(EXP,-n_T)


cc.wider <- cell.counts %>% ungroup() %>% select(EXP,n,Sample_ID,BATCH, treats) %>%
  #filter(!is.na(Origin2)) %>% 
  group_by(EXP,Sample_ID,BATCH, treats) %>%  
  pivot_wider(names_from=Origin2,values_from=c(n,Sample_ID),values_fill=list(n=0,Sample_ID="UNK")) %>%
  mutate(n_T=n_M+n_F) %>% 
  filter(n_T>500)  %>%
  arrange(EXP,-n_T)


write_tsv(cc.wider,paste0(outFolder,"cc.wider.v1.tsv"))

cv = read_tsv("./SampleListAndCV.txt")

aa <- left_join(aa,cv)

table(aa$Origin,aa$Origin2)
table(aa$Library_ID1,aa$EXP)


cell.counts <- aa %>% group_by(EXP,Sample_ID,Pregnancy_ID,Origin2,Library_ID1,SampleNum,Set,Preeclampsia,Combined) %>%
  summarize(n=n()) 

write_tsv(cell.counts,paste0(outFolder,"cell.counts.tsv"))


cc.wider <- cell.counts %>% ungroup() %>% select(EXP,Pregnancy_ID,Origin2,n,Sample_ID,Library_ID1,SampleNum,Set,Preeclampsia,Combined) %>%
  filter(!is.na(Origin2)) %>% 
  group_by(EXP,Pregnancy_ID,Origin2) %>%  
  pivot_wider(names_from=Origin2,values_from=c(n,Sample_ID),values_fill=list(n=0,Sample_ID="UNK")) %>%
  mutate(n_T=n_M+n_F) %>% 
  filter(n_T>300)  %>%
  arrange(SampleNum)

write_tsv(cc.wider,paste0(outFolder,"cc.wider.tsv"))

