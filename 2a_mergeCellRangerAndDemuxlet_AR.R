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
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","CZ1_group.txt") #for testing
base <- args[1]
outFolder=paste0(base,"2.1_mergeCellRangerAndDemuxlet_renamed/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

#read in samples file (just list of samples to run, each sample on newline)
if(!is.na(args[2])){
samples=read.table(paste0(base,args[2]),header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
}

basefolder=gsub("analysis/","counts_cellranger_hg38/",base)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

libList <- scan(paste0(basefolder,"libList.txt"),what=character(0))
if(!is.na(args[2])){
  libList <- libList[libList %in% samples$V1]
}

cat("creating seurat object with ",libList)

sc_list<-sapply(libList, function(x){
##  x<-libList[1]
  cat("#Processing: ",x,"\n")
  gp.data<- Read10X(data.dir = paste0(basefolder,x,"/filtered_feature_bc_matrix"))
  #################################################################################
  # creating seurat object
  #################################################################################
  sc <- CreateSeuratObject(counts = gp.data, project = "cellranger-CZI",min.cells = 3, min.features=200)
  sc@meta.data$Library<-rep(x,nrow(sc@meta.data))
  sc
} )

## find matching barcodes demuxlet and the sc object. 

opfn <- paste0(outFolder,"seuratObj-merge.",Sys.Date(),".rds") 
write_rds(sc_list, opfn)

# read it again if needed
#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-merge."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc_list <- read_rds(opfn)

sc <- merge(sc_list[[1]],sc_list[-1],add.cell.ids = libList, project="cellranger-CZI")

opfn <- paste0(outFolder,"seuratObj-all-ulnist-prior-to-demux.",Sys.Date(),".rds") 
write_rds(sc, opfn)


 ###################################################################################
 ############## get raw cellranger stats before merging with demuxlet ##############
 ###################################################################################

cat("get raw cellranger stats before merging with demuxlet")

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-all-ulnist-prior-to-demux."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

count <- sc@assays$RNA@counts

anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
               uns=grepl("S-",rn),rnz=rowSums(count))

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

fig0 <- ggplot(dd,aes(x=ident, y=percent.mt, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: percent mitochondria (mean)")+
        geom_text(aes(label=round(percent.mt)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure2.1_percent_mt.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

###(1), barcodes for each experiment           
fig0 <- ggplot(dd,aes(x=ident, y=ncell, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #Barcodes per experiment")+
        geom_text(aes(label=ncell),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure1.1_barcodes.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads)
dd2 <- dd%>%dplyr::select(ident,ngene,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="CellRanger: #UMIs per cell", "2"="CellRanger: #Genes per cell"))
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
png(paste0(figuredir,"Figure1.2_genes.png"), width=3000, height=2500, res=240)
print(fig0)
dev.off() 

## separate plots for reads and gene numbers
fig0 <- ggplot(dd,aes(x=ident, y=reads, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #UMI per cell")+
        geom_text(aes(label=round(reads)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure1.3_UMI_numb.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

fig0 <- ggplot(dd,aes(x=ident, y=ngene, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("CellRanger: #Gene per cell")+
        geom_text(aes(label=round(ngene)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure1.4_Gene_numb.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()


fig0 <- VlnPlot(sc, features = "percent.mt", ncol = 1, group.by='Library',pt.size = FALSE)+ scale_y_continuous(limits = c(0,50))
png(paste0(figuredir,"Figure0.1_violin_percent_mt.png"), width=4000, height=1000, res=120)
print(fig0)
dev.off()


fig0 <- VlnPlot(sc, features = "nFeature_RNA", ncol = 1, group.by='Library',pt.size = FALSE)
png(paste0(figuredir,"Figure0.2_violin_nfeatures_genes.png"), width=4000, height=1000, res=120)
print(fig0)
dev.off()


fig0 <- VlnPlot(sc, features = "nCount_RNA", ncol = 1, group.by='Library',pt.size = FALSE)
png(paste0(figuredir,"Figure0.3_violin_ncount_umi.png"), width=4000, height=1000, res=120)
print(fig0)
dev.off()

#############################################
#  Scatter plots of kallisto vs cell ranger #
#############################################

cat("making scatter plots of kallisto vs cell ranger")
### rename dd

#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell_CR=n(),
                       reads_CR=mean(nCount_RNA),
                       ngene_CR=mean(nFeature_RNA),
                       percent.mt_CR=mean(percent.mt),
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       #.groups="drop"
                       )
dd <- dd%>%dplyr::rename(ident=Library) %>%
           mutate(batch=gsub("-.*","",ident))

### load the kallisto dd to merge and do scatter plots
ddk <- read.csv(paste0(base,"2b_mergeKallistoAndDemuxlet/raw-stats-kallisto-lib.csv"))#, row.names=F, quote=FALSE)

mdd <- left_join(as.data.frame(dd), ddk, by="ident")

fig0 <- ggplot(mdd, aes(x=ncell_CR, y=ncell_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("number of cell per library - CellRanger vs kallisto")) +
            theme_minimal() +
            geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of cells", y="Kallisto # of cells")
png(paste0(figuredir,"Figure3.0-scatter-numb-cell-CR-KL.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

fig0 <- ggplot(mdd, aes(x=ncell_CR, y=ncell_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("number of cell per library - CellRanger vs kallisto")) +
            theme_minimal() +
            geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of cells", y="Kallisto # of cells")
png(paste0(figuredir,"Figure3.0-scatter-numb-cell-CR-KL-labled.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### UMI per cell scatter
fig0 <- ggplot(mdd, aes(x=reads_CR, y=reads_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("number of UMI per cell - CellRanger vs kallisto")) +
            theme_minimal() +
            geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of UMI", y="Kallisto # of UMI")
png(paste0(figuredir,"Figure3.1-scatter-numb-readUMI-per-cell-CR-KL-labled.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### genes per cell scatter
fig0 <- ggplot(mdd, aes(x=ngene_CR, y=ngene_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("number of genes per cell - CellRanger vs kallisto")) +
            theme_minimal() +
            geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of genes", y="Kallisto # of genes")
png(paste0(figuredir,"Figure3.2-scatter-numb-ngenes-per-cell-CR-KL-labled.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### percent mitochondria
fig0 <- ggplot(mdd, aes(x=percent.mt_CR, y=percent.mt_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("Percent Mitochondria (mean) - CellRanger vs kallisto")) +
            theme_minimal() +
            #geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger percent mitochondria", y="Kallisto percent mitochondria")
png(paste0(figuredir,"Figure3.3-scatter-percent-mitochondria-CR-KL-labled.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

############# find matching barcodes demuxlet and the sc object. 

cat("find matching barcodes demuxlet and the sc object")

opfn <- paste0(base,"/1_demux_output/1_demux_New.SNG.rds")
demux <- read_rds(opfn)
###

# filter 
demux <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID) 

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

opfn <- paste0(outFolder,"seuratObj-merge.md.post-merge-demux.",Sys.Date(),".rds") 
write_rds(sc, opfn)


 ###########################################################################
 ##############  cellranger stats after merging with demuxlet ##############
 ###########################################################################

cat("cellranger stats after merging with demuxlet")

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-merge.md.post-merge-demux."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)

count <- sc@assays$RNA@counts

anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
               uns=grepl("S-",rn),rnz=rowSums(count))

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

png(paste0(figuredir,"Figure4.1_barcodes_post_merge.png"), width=1000, height=600, res=120)
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
png(paste0(figuredir,"Figure4.2_genes_post_merge.png"), width=3000, height=2500, res=240)
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

png(paste0(figuredir,"Figure4.3_UMI_numb_post_merge.png"), width=1000, height=600, res=120)
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

png(paste0(figuredir,"Figure4.4_Gene_numb_post_merge.png"), width=1000, height=600, res=120)
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

png(paste0(figuredir,"Figure4.5_percent_mt_merge.png"), width=1000, height=600, res=120)
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

sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

mean(sc[["percent.mt"]]<10) # 
mean(sc[["percent.mt"]]<15) # 
mean(sc[["percent.mt"]]<20) # 

sc[["nFeature_RNA"]] %>% summary()
mean(sc[["nFeature_RNA"]]>200) # 
sc[["nCount_RNA"]] %>% summary()
mean(sc[["nCount_RNA"]] < 20000) # 

scsub <- subset(sc, subset = percent.mt < 10 & nFeature_RNA > 10000) 
opfn <- paste0(outFolder,"seuratObj-postmerge-greater10kfeature.",Sys.Date(),".rds") 
write_rds(scsub,opfn)
rm(scsub)
gc()

sc <- subset(sc, subset = nFeature_RNA > 200 & percent.mt < 10) #& nFeature_RNA < 20000

opfn <- paste0(outFolder,"seuratObj-postmerge-after-mt-filtering.",Sys.Date(),".rds") 
write_rds(sc, opfn)


 ###########################################################################
 ##############  cellranger stats after removing high mithoch ##############
 ###########################################################################

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-postmerge-after-mt-filtering."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)


cat("cellranger stats after filters")

count <- sc@assays$RNA@counts

anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
               uns=grepl("S-",rn),rnz=rowSums(count))

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

png(paste0(figuredir,"Figure6.1_barcodes_post_merge_post_filt.png"), width=1000, height=600, res=120)
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
png(paste0(figuredir,"Figure6.2_genes_post_merge_post_filt.png"), width=3000, height=2500, res=240)
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

png(paste0(figuredir,"Figure6.3_UMI_numb_post_merge_post_filt.png"), width=1000, height=600, res=120)
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

png(paste0(figuredir,"Figure6.4_Gene_numb_post_merge_post_filt.png"), width=1000, height=600, res=120)
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

png(paste0(figuredir,"Figure6.5_percent_mt_merge_post_filt.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()


#############################################
#  Scatter plots of kallisto vs cell ranger post merge
#############################################

cat("Scatter plots of kallisto vs cell ranger post merge")

### rename dd
#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell_CR=n(),
                       reads_CR=mean(nCount_RNA),
                       ngene_CR=mean(nFeature_RNA),
                       percent.mt_CR=mean(percent.mt),
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       #.groups="drop"
                       )
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))

### load the kallisto dd to merge and do scatter plots
ddk <- read.csv(paste0(base,"2b_mergeKallistoAndDemuxlet/post-merge-pos-mt-filter-kallisto-lib.csv"))#, row.names=F, quote=FALSE)

mdd <- left_join(as.data.frame(dd), ddk, by="ident")

##scatter plots

fig0 <- ggplot(mdd, aes(x=ncell_CR, y=ncell_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("number of cell per library - CellRanger vs kallisto")) +
            theme_minimal() +
            geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of cells", y="Kallisto # of cells")
png(paste0(figuredir,"Figure7.0-scatter-numb-cell-CR-KL.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

fig0 <- ggplot(mdd, aes(x=ncell_CR, y=ncell_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("number of cell per library - CellRanger vs kallisto")) +
            theme_minimal() +
            geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of cells", y="Kallisto # of cells")
png(paste0(figuredir,"Figure7.0-scatter-numb-cell-CR-KL-labled.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()


### UMI per cell scatter
fig0 <- ggplot(mdd, aes(x=reads_CR, y=reads_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("number of UMI per cell - CellRanger vs kallisto")) +
            theme_minimal() +
            #geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of UMI", y="Kallisto # of UMI")
png(paste0(figuredir,"Figure7.1-scatter-numb-readUMI-per-cell-CR-KL.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()


### genes per cell scatter
fig0 <- ggplot(mdd, aes(x=ngene_CR, y=ngene_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("number of genes per cell - CellRanger vs kallisto")) +
            theme_minimal() +
            #geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of genes", y="Kallisto # of genes")
png(paste0(figuredir,"Figure7.2-scatter-numb-ngenes-per-cell-CR-KL.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()


### percent mitochondria
fig0 <- ggplot(mdd, aes(x=percent.mt_CR, y=percent.mt_KL, color=batch.x))+ 
            geom_point() + 
            ggtitle(paste0("Percent Mitochondria (mean) - CellRanger vs kallisto")) +
            theme_minimal() +
            #geom_abline(intercept = 0, slope = 1)+
            #geom_smooth(method='lm')+
            #labs(fill = "LPS sig genes") +
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            #geom_text(aes(label=ident),  position = position_dodge(width = 1))+
            #geom_text_repel(aes(label=ident),vjust=-0.7, size=2.5)+
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger percent mitochondria", y="Kallisto percent mitochondria")
png(paste0(figuredir,"Figure7.3-scatter-percent-mitochondria-CR-KL.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()
