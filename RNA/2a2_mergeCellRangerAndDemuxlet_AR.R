library(Seurat)
library(Matrix)
library(future)
library(readr)
library(tidyverse)
library(harmony)


#####################################################################
### 07/25/2023, Ali R                                          ###### 
###  SCAIP7-18 CellRanger alignment and demultiplexing         ######
###  comparison of CellRanger and Kallisto                     ######
###  include post renaming EtOH -> CTRL barcod reruns          ######
###  modified from pre-eclampsia paper:                        ######
###                    2_mergeCellRangerAndDemuxlet.R          ######
#####################################################################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","fastdemux") #for testing "CZI2_group.txt"
#args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/", "demux")
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","fastdemux")

base <- args[1]
method <- args[2]
#read in samples file (just list of samples to run, each sample on newline)
if(!is.na(args[3])){
samples=read.table(paste0(base,args[3]),header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
project=sapply(strsplit(args[3],"_"),function(y)y[1])
cat("samplefile= ",args[3])
}else{
   project="ALL"
}
if(method=="demux"){
   outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
   demux_in <- paste0(base,"1_demux_output/")
   kallisto_in <- paste0("2b_mergeKallistoAnd",method,"/")
   opfn <- paste0(demux_in,project,".1_demux_filt.SNG.rds")
   demux <- read_rds(opfn)
} else {
   outFolder=paste0(base,"2.1_mergeCellRangerAnd",method,"/")
   demux_in <- paste0(base,"1_demux_alt_output/")
   opfn <- paste0(demux_in,project,".1_demux_alt_filt.SNG.rds")
   demux <- read_rds(opfn)   
   kallisto_in <- paste0("2b_mergeKallistoAnd",method,"/")
}

if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_cellranger_2024-04-19/",base)  ##

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

libList <- scan(paste0(basefolder,"libList.txt"),what=character(0)) 
if(!is.na(args[3])){
  libList <- libList[libList %in% samples$V1]
}


##Removig DEX libraries here here:
libList <- libList[grep("DEX",libList,invert=TRUE)]


## find matching barcodes demuxlet and the sc object. 

opfn <- paste0(outFolder,project,".seuratObj-merge.","2025-11-10",".rds")

sc_list <- read_rds(opfn)

sc <- merge(sc_list[[1]],sc_list[-1],add.cell.ids = libList, project=paste0("cellranger-CZI.",project))

rm(sc_list)

opfn <- paste0(outFolder,project,".seuratObj-all-ulnist-prior-to-demux.",Sys.Date(),".rds") 
sc <- read_rds(opfn)

### RPR 2025-11-10  Error: vector::reserve
##sc[["RNA"]] <- JoinLayers(sc[["RNA"]])

opfn <- paste0(outFolder,project,".seuratObj-all-ulnist-prior-to-demux.",Sys.Date(),".rds") 
write_rds(sc, opfn)

############# find matching barcodes demuxlet and the sc object. 

cat("find matching barcodes demuxlet and the sc object")

###

# filter 

demux$NEW_BARCODE = paste0(demux$NEW_BARCODE,"-1")

matching_bc = intersect(colnames(sc),demux$NEW_BARCODE)

length(matching_bc) 
length(colnames(sc)) 
length(rownames(demux)) 

sc <- subset(sc, cells = matching_bc)

md <- sc@meta.data %>% rownames_to_column("NEW_BARCODE")
md <- left_join(md,demux)
#md <- left_join(md,cv)

dim(md)

##sc <- subset(sc, subset = nFeature_RNA > 200) # change to: #>200  #<20000 , #doubletfinder #soupx  #diem

md <- md %>% column_to_rownames("NEW_BARCODE")
md$NEW_BARCODE = rownames(md)

stopifnot(identical(rownames(sc@meta.data),md$NEW_BARCODE))

sc@meta.data=md

opfn <- paste0(outFolder,project,".seuratObj-merge.md.post-merge-demux.",Sys.Date(),".rds") 
write_rds(sc, opfn)


 ###########################################################################
 ##############  cellranger stats after merging with demuxlet ##############
 ###########################################################################

cat("cellranger stats after merging with demuxlet")

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-merge.md.post-merge-demux."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)

## RPR 2025-11-10 I removed the statistics before merging with demuxlet, but I need to add back this here. 
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")


## RPR 2025-11-10  Maybe not needed so I'm trying ithout it. 
###count <- sc[["RNA"]]$counts

#anno <- data.frame(rn=rownames(count))%>%
#        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
#               uns=grepl("S-",rn),rnz=rowSums(count))

meta <- sc@meta.data
#meta$nCount_spliced <- nCount_spliced
#meta$nFeature_spliced <- nFeature_spliced

#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell=n(),
                       reads=mean(nCount_RNA),
                       ngene=mean(nFeature_RNA),
                       percent.mt=mean(percent.mt),
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       #.groups="drop"
                       )
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))

              
###(1), barcodes for each experiment           
fig0 <- ggplot(dd,aes(x=ident, y=ncell, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #Barcodes per experiment - merged")+
        geom_text(aes(label=ncell),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,project,".Figure4.1_barcodes_post_merge.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads)
dd2 <- dd%>%dplyr::select(ident,ngene,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="CellRanger: #UMIs per cell - merged", "2"="CellRanger: #Genes per cell - merged"))
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
png(paste0(figuredir,project,".Figure4.2_genes_post_merge.png"), width=3000, height=2500, res=240)
print(fig0)
dev.off() 

## separate plots for reads and gene numbers
fig0 <- ggplot(dd,aes(x=ident, y=reads, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #UMI per cell - merged")+
        geom_text(aes(label=round(reads)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,project,".Figure4.3_UMI_numb_post_merge.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

fig0 <- ggplot(dd,aes(x=ident, y=ngene, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #Gene per cell - merged")+
        geom_text(aes(label=round(ngene)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,project,".Figure4.4_Gene_numb_post_merge.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

fig0 <- ggplot(dd,aes(x=ident, y=percent.mt, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: percent mitochondria (mean) - merge")+
        geom_text(aes(label=round(percent.mt)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,project,".Figure4.5_percent_mt_merge.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

# cell ranger stats post merge before mt filter
sum(dd$ncell)
median(dd$ncell)
median(dd$reads)
median(dd$ngene)
dim(dd)

#################################################################
#filters post demux merge
cat("post demux merge")

## RPR commenting out bc already done above. 
##sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

mean(sc[["percent.mt"]]<10) # 
mean(sc[["percent.mt"]]<15) # 
mean(sc[["percent.mt"]]<20) # 

sc[["nFeature_RNA"]] %>% summary()
mean(sc[["nFeature_RNA"]]>200) # 
sc[["nCount_RNA"]] %>% summary()
mean(sc[["nCount_RNA"]] < 20000) # 

## RPR verify methods on this. or is this not used at all? AR: going over paper methods and other downstream scripts, this output is not used again. 
scsub <- subset(sc, subset = percent.mt < 10 & nFeature_RNA > 10000) 
opfn <- paste0(outFolder,project,".seuratObj-postmerge-greater10kfeature.",Sys.Date(),".rds") 
write_rds(scsub,opfn)
rm(scsub)
gc()

sc <- subset(sc, subset = nFeature_RNA > 200 & percent.mt < 10) #& nFeature_RNA < 20000

opfn <- paste0(outFolder,project,".seuratObj-postmerge-after-mt-filtering.",Sys.Date(),".rds") 
write_rds(sc, opfn)


 ###########################################################################
 ##############  cellranger stats after removing high mithoch ##############
 ###########################################################################

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-postmerge-after-mt-filtering."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)


cat("cellranger stats after filters")

## RPR 2025-11-10 probably not needed. 
##count <- sc[["RNA"]]$counts

#anno <- data.frame(rn=rownames(count))%>%
#        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
#               uns=grepl("S-",rn),rnz=rowSums(count))

meta <- sc@meta.data
#meta$nCount_spliced <- nCount_spliced
#meta$nFeature_spliced <- nFeature_spliced

#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell=n(),
                       reads=mean(nCount_RNA),
                       ngene=mean(nFeature_RNA),
                       percent.mt=mean(percent.mt),
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       #.groups="drop"
                       )
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))

# cell ranger stats post merge after mt filter
sum(dd$ncell)
median(dd$ncell)
median(dd$reads)
median(dd$ngene)
dim(dd)

###(1), barcodes for each experiment           
fig0 <- ggplot(dd,aes(x=ident, y=ncell, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #Barcodes per experiment - merged")+
        geom_text(aes(label=ncell),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,project,".Figure6.1_barcodes_post_merge_post_filt.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads)
dd2 <- dd%>%dplyr::select(ident,ngene,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="CellRanger: #UMIs per cell - filtr", "2"="CellRanger: #Genes per cell - filtr"))
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
png(paste0(figuredir,project,".Figure6.2_genes_post_merge_post_filt.png"), width=3000, height=2500, res=240)
print(fig0)
dev.off() 


## separate plots for reads and gene numbers
fig0 <- ggplot(dd,aes(x=ident, y=reads, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #UMI per cell - filtr")+
        geom_text(aes(label=round(reads)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,project,".Figure6.3_UMI_numb_post_merge_post_filt.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

fig0 <- ggplot(dd,aes(x=ident, y=ngene, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #Gene per cell - filtr")+
        geom_text(aes(label=round(ngene)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,project,".Figure6.4_Gene_numb_post_merge_post_filt.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()


# mitochondrial content
fig0 <- ggplot(dd,aes(x=ident, y=percent.mt, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: percent mitochondria (mean) - filtr")+
        geom_text(aes(label=round(percent.mt)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,project,".Figure6.5_percent_mt_merge_post_filt.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()


