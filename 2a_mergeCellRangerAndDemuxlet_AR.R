library(Seurat)
library(Matrix)
library(future)
library(readr)
library(tidyverse)


#####################################################################
### 07/25/2023, Ali R                                          ###### 
###  SCAIP7-18 CellRanger alignment and demultiplexing         ######
###  comparison of CellRanger and Kallisto                     ######
###  include post renaming EtOH -> CTRL barcod reruns          ######
###  modified from pre-eclampsia paper:                        ######
###                    2_mergeCellRangerAndDemuxlet.R          ######
#####################################################################


setwd("/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/")
outFolder="./2.1_mergeCellRangerAndDemuxlet_renamed/"
system(paste0("mkdir -p ", outFolder))

future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 3)

basefolder <- "/nfs/rprdata/SCAIP-ALOFT2/counts_cellranger_2023-05-27/"

libList <- scan("/nfs/rprdata/SCAIP-ALOFT2/counts_cellranger_2023-05-27/libList.txt",what=character(0))

#cv <- read_tsv("SampleCV.txt") %>% select(EXP,Pregnancy_ID, Library_ID1, Set, Preeclampsia, Combined) #?
#"./1_demux_output/1_demux_New.ALL.rds"

##"1_demux_output/1_demux_New.ALL.rds"
opfn <- "./1_demux_output/1_demux_New.SNG.rds"
demux <- read_rds(opfn)
###
demux <- demux %>% mutate(BATCH=gsub("-.*", "", EXP), treats=gsub(".*[0-9].{,2}-","",EXP)) 
head(demux)

# correct the EtOH ETOH discrepency in treat column 
#demux$treats <- gsub("EtOH", "ETOH", demux$treats)

# filter 
demux <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID=SNG.BEST.GUESS) 

## filter for right combinations. 
#expected_combos=paste0(cv$EXP,"_",cv$Pregnancy_ID)
#demux <- demux %>% filter(paste0(EXP,"_",Pregnancy_ID) %in% expected_combos)


sc_list<-sapply(libList, function(x){
##  x<-samples[1]
  cat("#Processing: ",x,"\n")
  gp.data<- Read10X(data.dir = paste0(basefolder,x,"/filtered_feature_bc_matrix"))
  #################################################################################
  # creating seurat object
  #################################################################################
  sc <- CreateSeuratObject(counts = gp.data, project = "scALOFT-SCAIP7-18",min.cells = 3, min.features=200)
  sc@meta.data$Library<-rep(x,nrow(sc@meta.data))
  sc
} )

#?
#sc <- merge(sc_list[[1]],sc_list[-1],add.cell.ids = libList, project="scALOFT-SCAIP7-18")


## find matching barcodes demuxlet and the sc object. 
sc <- sc_list

opfn <- paste0(outFolder,"seuratObj-merge.",Sys.Date(),".rds") 
write_rds(sc, opfn)

# read it again if needed
#opfn <- paste0(outFolder,"seuratObj-merge.",Sys.Date(),".rds") 
opfn <- "./2.1_mergeCellRangerAndDemuxlet_renamed/seuratObj-merge.2023-07-24.rds"
sc_list <- read_rds(opfn)

#sc <- merge(sc_list[[1]],sc_list[2],add.cell.ids = libList[1:2], project="scALOFT-SCAIP7-18")

sc <- merge(sc_list[[1]],sc_list[-1],add.cell.ids = libList, project="scALOFT-SCAIP7-18")

opfn <- paste0(outFolder,"seuratObj-all-ulnist-prior-to-demux.",Sys.Date(),".rds") 
write_rds(sc, opfn)





 ###################################################################################
 ############## get raw cellranger stats before merging with demuxlet ##############
 ###################################################################################

#test <- read_rds("./2_mergeCellRangerAndDemuxlet/seuratObj-all-ulnist-prior-to-demux.2023-06-23.rds")
#sc <- test

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
        ggtitle("CellRanger: #Barcodes per experiment")+
        geom_text(aes(label=ncell),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure1.1_barcodes.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure1.2_genes.png", width=3000, height=2500, res=240)
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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure1.3_UMI_numb.png", width=1000, height=600, res=120)
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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure1.4_Gene_numb.png", width=1000, height=600, res=120)
print(fig0)
dev.off()


#############################################
#  Scatter plots of kallisto vs cell ranger #
#############################################

### rename ddd

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
ddk <- read.csv("./2b_mergeKallistoAndDemuxlet/raw-stats-kallisto-48lib.csv")#, row.names=F, quote=FALSE)

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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure3.0-scatter-numb-cell-CR-KL.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure3.0-scatter-numb-cell-CR-KL-labled.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure3.1-scatter-numb-readUMI-per-cell-CR-KL-labled.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure3.2-scatter-numb-ngenes-per-cell-CR-KL-labled.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure3.3-scatter-percent-mitochondria-CR-KL-labled.png", width=1000, height=600, res=120)
print(fig0)
dev.off()





####################################################################################
# percentage of reads mapping to mitochondrial - before merging with demulet
####################################################################################

refFolder="/wsu/home/groups/piquelab/data/refGenome10x/refdata-gex-GRCh38-2020-A/"

cmd <- paste0("cat ",refFolder,"/genes/genes.gtf",
              " | awk '$3~/gene/'",
              " | sed 's/gene_id //;s/;.* gene_name /\t/;s/;.*transcript_biotype/\t/;s/;.*//'")
cat(cmd,"\n")


aux <- read_tsv(pipe(cmd),col_names = FALSE) %>% mutate(TSS=ifelse(X7=="+",X4,X5)) %>%
  dplyr::select(Chr=X1,Min=X4,Max=X5,kbid=X9,TSS,Strand=X7,gene_name=X10) 

anno <- tibble(gene_name=rownames(sc),rs=rowSums(sc@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>% filter(!is.na(Chr))

table(is.na(anno$Chr))
table(anno$Chr)

head(anno)
head(aux)

sc <- sc[anno$gene_name,]

sc[["percent.mt"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrM",]$gene_name)

sc[["percent.mt"]] %>% summary()

# prorportion of cells with < xx percent mitochondria for CellRanger
mean(sc[["percent.mt"]]<10) # 0.9717731
mean(sc[["percent.mt"]]<15) # 0.9869789
mean(sc[["percent.mt"]]<20) # 0.9916916

sc[["nFeature_RNA"]] %>% summary()
sc[["nCount_RNA"]] %>% summary()


fig0 <- VlnPlot(sc, features = "percent.mt", ncol = 1, group.by='Library')
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure0.1_violin_percent_mt.png", width=4000, height=1000, res=120)
print(fig0)
dev.off()


fig0 <- VlnPlot(sc, features = "nFeature_RNA", ncol = 1, group.by='Library')
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure0.2_violin_nfeatures_genes.png", width=4000, height=1000, res=120)
print(fig0)
dev.off()


fig0 <- VlnPlot(sc, features = "nCount_RNA", ncol = 1, group.by='Library')
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure0.3_violin_ncount_umi.png", width=4000, height=1000, res=120)
print(fig0)
dev.off()



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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure2.1_percent_mt.png", width=1000, height=600, res=120)
print(fig0)
dev.off()


#fig0 <- VlnPlot(sc, features =  "percent.mt", ncol = 28, fill=Library)
#png("./2_mergeCellRangerAndDemuxlet/Figure2.2_percent_mt_violin_plot.png", width=1000, height=600, res=120)
#print(fig0)
#dev.off()


############# find matching barcodes demuxlet and the sc object. 

##"1_demux_output/1_demux_New.ALL.rds"
opfn <- "./1_demux_output/1_demux_New.SNG.rds"
demux <- read_rds(opfn)
###
demux <- demux %>% mutate(BATCH=gsub("-.*", "", EXP), treats=gsub(".*[0-9].{,2}-","",EXP)) 
head(demux)

# correct the EtOH ETOH discrepency in treat column 
#demux$treats <- gsub("EtOH", "ETOH", demux$treats)

# filter 
demux <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID=SNG.BEST.GUESS) 


demux$NEW_BARCODE = paste0(demux$NEW_BARCODE,"-1")

matching_bc = intersect(colnames(sc),demux$NEW_BARCODE)


length(matching_bc) # 570427
length(colnames(sc)) # 710848
length(rownames(demux)) # 647853

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

#sc <- read_rds("/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/2_mergeCellRangerAndDemuxlet/seuratObj-merge.md.post-merge-demux.2023-06-23.rds")
#test <- read_rds("/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/2_mergeCellRangerAndDemuxlet/seuratObj-before-clustering.2023-06-26.rds")
#sc <- read_rds("/wsu/home/groups/piquelab/SCAIP_2022/Ali/scALOFT/2_mergeCellRangerAndDemuxlet/seuratObj-after-mt-filtering.2023-06-26.rds")



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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure4.1_barcodes_post_merge.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure4.2_genes_post_merge.png", width=3000, height=2500, res=240)
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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure4.3_UMI_numb_post_merge.png", width=1000, height=600, res=120)
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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure4.4_Gene_numb_post_merge.png", width=1000, height=600, res=120)
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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure4.5_percent_mt_merge.png", width=1000, height=600, res=120)
print(fig0)
dev.off()

# cell ranger stats post merge before mt filter
sum(dd$ncell)
median(dd$ncell)
median(dd$reads)
median(dd$ngene)
dim(dd)



#############################################
#  Scatter plots of kallisto vs cell ranger post merge
#############################################

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
ddk <- read.csv("./2b_mergeKallistoAndDemuxlet/post-merge-kallisto-48lib.csv")#, row.names=F, quote=FALSE)

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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure5.0-scatter-numb-cell-CR-KL.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure5.0-scatter-numb-cell-CR-KL-labled.png", width=1000, height=600, res=120)
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
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of UMI", y="Kallisto # of UMI")
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure5.1-scatter-numb-readUMI-per-cell-CR-KL.png", width=1000, height=600, res=120)
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
            #geom_text(aes(label=ident),vjust=-1.2, size=2.5,  position = position_dodge(width = 1))+
            theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), plot.subtitle=element_text(hjust=0.5), ) +
            #scale_fill_discrete(name = "Differentially Expressed Genes") +
            #guides(fill=guide_legend(title="Differentially Expressed Genes")) +
            #labs(color = "DEG", subtitle = "all participants")
            #labs(subtitle = paste("correlation for ", nrow(mergedmale), " genes", sep="")) +
            labs(x="CellRanger # of genes", y="Kallisto # of genes")
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure5.2-scatter-numb-ngenes-per-cell-CR-KL.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure5.3-scatter-percent-mitochondria-CR-KL.png", width=1000, height=600, res=120)
print(fig0)
dev.off()


#################################################################
# Quality control
#################################################################

####################################################
# percentage of reads mapping to mitochondrial
####################################################

refFolder="/wsu/home/groups/piquelab/data/refGenome10x/refdata-gex-GRCh38-2020-A/"

cmd <- paste0("cat ",refFolder,"/genes/genes.gtf",
              " | awk '$3~/gene/'",
              " | sed 's/gene_id //;s/;.* gene_name /\t/;s/;.*transcript_biotype/\t/;s/;.*//'")
cat(cmd,"\n")


aux <- read_tsv(pipe(cmd),col_names = FALSE) %>% mutate(TSS=ifelse(X7=="+",X4,X5)) %>%
  dplyr::select(Chr=X1,Min=X4,Max=X5,kbid=X9,TSS,Strand=X7,gene_name=X10) 


anno <- tibble(gene_name=rownames(sc),rs=rowSums(sc@assays$RNA@data)) %>% filter(rs>0) %>% left_join(aux) %>% filter(!is.na(Chr))

table(is.na(anno$Chr))

table(anno$Chr)

##table(anno$Type)

head(anno)

head(aux)


sc <- sc[anno$gene_name,]


sc[["percent.mt"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrM",]$gene_name)

sc[["percent.mt"]] %>% summary()

mean(sc[["percent.mt"]]<10) # 0.9717731
mean(sc[["percent.mt"]]<15) # 0.9869789
mean(sc[["percent.mt"]]<20) # 0.9916916


sc[["nFeature_RNA"]] %>% summary()
mean(sc[["nFeature_RNA"]]>200) # 0.9997327


sc[["nCount_RNA"]] %>% summary()
mean(sc[["nCount_RNA"]] < 20000) # 0.9341406


dim(sc)

sc <- subset(sc, subset = nFeature_RNA > 200 & percent.mt < 20) #& nFeature_RNA < 20000

dim(sc) #567386

opfn <- paste0(outFolder,"seuratObj-after-mt-filtering.",Sys.Date(),".rds") 
write_rds(sc, opfn)

##sc<-read_rds("4_harmony_cellClass_doubletfinder_chrM/seuratObj-after-mt-filtering.2022-12-29.rds")


 ###########################################################################
 ##############  cellranger stats after removing high mithoch ##############
 ###########################################################################

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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure6.1_barcodes_post_merge_post_filt.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure6.2_genes_post_merge_post_filt.png", width=3000, height=2500, res=240)
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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure6.3_UMI_numb_post_merge_post_filt.png", width=1000, height=600, res=120)
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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure6.4_Gene_numb_post_merge_post_filt.png", width=1000, height=600, res=120)
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

png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure6.5_percent_mt_merge_post_filt.png", width=1000, height=600, res=120)
print(fig0)
dev.off()




#############################################
#  Scatter plots of kallisto vs cell ranger post merge
#############################################

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
ddk <- read.csv("./2b_mergeKallistoAndDemuxlet/post-merge-pos-mt-filter-kallisto-48lib.csv")#, row.names=F, quote=FALSE)

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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure7.0-scatter-numb-cell-CR-KL.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure7.0-scatter-numb-cell-CR-KL-labled.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure7.1-scatter-numb-readUMI-per-cell-CR-KL.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure7.2-scatter-numb-ngenes-per-cell-CR-KL.png", width=1000, height=600, res=120)
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
png("./2.1_mergeCellRangerAndDemuxlet_renamed/Figure7.3-scatter-percent-mitochondria-CR-KL.png", width=1000, height=600, res=120)
print(fig0)
dev.off()









############################################################

## Clustering
############################################################

### load the seurat object

opfn <- paste0(outFolder,"seuratObj-after-mt-filtering.2023-07-25.rds") 
sc <- read_rds(opfn)

##sc<-read_rds("4_harmony_cellClass_doubletfinder_chrM/seuratObj-after-mt-filtering.2022-12-29.rds")


future::plan(strategy = 'multicore', workers = 16)
options(future.globals.maxSize = 30 * 1024 ^ 6)


sc <- NormalizeData(sc, verbose=TRUE)

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 3000)

sc <- ScaleData(sc, verbose = TRUE)

sc <- RunPCA(sc,pc.genes = sc@var.genes, npcs = 100, verbose = TRUE)

library(harmony)
sc <- RunHarmony(sc,c("Library"),reduction="pca")

sc <- RunUMAP(sc,reduction = "harmony", dims = 1:30)

sc <- FindNeighbors(sc, reduction = "harmony", dims = 1:30, verbose = TRUE)

#before clustering
opfn <- paste0(outFolder,"seuratObj-before-clustering.",Sys.Date(),".rds") 
write_rds(sc, opfn)


sc <- FindClusters(sc, verbose = TRUE,resolution=0.5)

opfn <- paste0(outFolder,"seuratObj-post-clustering-res0.5.",Sys.Date(),".rds")
write_rds(sc, opfn)


### read the seurat object
opfn <- paste0(outFolder,"seuratObj-post-clustering-res0.5.2023-06-27.rds")
sc <- read_rds(opfn)



fname=paste0(outFolder,"UMAP_Harmony-res0.5.",Sys.Date(),".png");
ggdp = DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 6) #+ NoLegend()
ggdp
dev.off()


# make initial umap group by treatment
fname=paste0(outFolder,"Figure5.1_UMAP_Harmony-res0.5_group_treatment",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="treats")+ #, cols=col0)+
        theme_bw()+
        theme(legend.position="right",
              axis.title=element_text(size=25),
              axis.text=element_text(size=25),
              legend.text=element_text(size=25), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(), 
              plot.title=element_text(size=35),           
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()


# make initial umap group by batch
fname=paste0(outFolder,"Figure5.2_UMAP_Harmony-res0.5_group_batch",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="BATCH")+ #, cols=col0)+
        theme_bw()+
        theme(legend.position="right",
              axis.title=element_text(size=25),
              axis.text=element_text(size=25),
              legend.text=element_text(size=25), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(), 
              plot.title=element_text(size=35),           
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()


# make initial umap group by library
fname=paste0(outFolder,"Figure5.3_UMAP_Harmony-res0.5_group_Library",Sys.Date(),".png");
png(fname,width=5000,height=6000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=F, group.by="Library")+ #, cols=col0)+
        theme_bw()+
        theme(legend.position="bottom",
              axis.title=element_text(size=25),
              axis.text=element_text(size=25),
              legend.text=element_text(size=25), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(), 
              plot.title=element_text(size=35),           
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()



# make initial umap group by library
fname=paste0(outFolder,"Figure5.4_UMAP_Harmony-res0.5_group_seurat_cluster_with_names",Sys.Date(),".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", label=T, group.by="seurat_clusters", label.size=15)+ #, cols=col0)+
        theme_bw()+
        theme(legend.position="right",
              axis.title=element_text(size=25),
              axis.text=element_text(size=25),
              legend.text=element_text(size=25), 
              #legend.position=c(0.1, 0.85),
              #legend.text=element_text(size=12),
              #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
              #legend.key=element_blank(), #legend.key=element_rect(fill=NA),
              #legend.box.background=element_blank(), 
              plot.title=element_text(size=35),  
              text=element_text(size=35),         
              panel.border=element_rect(colour="black", fill=NA))
print(fig1)
dev.off()



aa <- FetchData(sc,c("UMAP_1","UMAP_2","BATCH","EXP","treats","Sample_ID", "seurat_clusters"))


#fname=paste0(outFolder,"UMAP_Harmony-res0.5.grid",".pdf");
#pdf(fname,width=14,height=5)
fname=paste0(outFolder,"Figure6.1_UMAP_Harmony-res0.5.grid-batch",Sys.Date(),".png");
png(fname,width=5000,height=4000, res=240)
    #fname=paste0(outFolder,"UMAP_LocationHarmony.Origin.png");
    #png(fname,width=1600,height=1200)
    p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
                geom_point(size=0.1) +
                ##    scale_color_manual(values=group.colors) +
##                guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
##                scale_color_manual("Origin",values=c("M"="#D1D1D1","F"="#A61BB5"))+
                facet_wrap(~BATCH) +
                theme_bw() +
                theme(legend.position="right",
                axis.title=element_text(size=25),
                axis.text=element_text(size=25),
                legend.text=element_text(size=25), 
                strip.text.x = element_text(size = 30),
                #legend.position=c(0.1, 0.85),
                legend.title=element_text(size=25),
                #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
                legend.key=element_blank(), #legend.key=element_rect(fill=NA),
                #legend.box.background=element_blank(), 
                plot.title=element_text(size=35),           
                panel.border=element_rect(colour="black", fill=NA))
                #guides(shape = guide_legend(override.aes = list(size = 35)))
    p2
    ##    theme_black()
dev.off()


#fname=paste0(outFolder,"UMAP_Harmony-res0.5.grid",".pdf");
#pdf(fname,width=14,height=5)
fname=paste0(outFolder,"Figure6.2_UMAP_Harmony-res0.5.grid-treats",Sys.Date(),".png");
png(fname,width=5000,height=4000, res=240)
    #fname=paste0(outFolder,"UMAP_LocationHarmony.Origin.png");
    #png(fname,width=1600,height=1200)
    p2 <- ggplot(aa,aes(UMAP_1,UMAP_2,color=seurat_clusters)) +
                geom_point(size=0.1) +
                ##    scale_color_manual(values=group.colors) +
##                guides(colour = guide_legend(override.aes = list(size=10),title="Cell origin")) +
##                scale_color_manual("Origin",values=c("M"="#D1D1D1","F"="#A61BB5"))+
                facet_wrap(~treats) +
                theme_bw() +
                theme(legend.position="right",
                axis.title=element_text(size=25),
                axis.text=element_text(size=25),
                legend.text=element_text(size=25), 
                strip.text.x = element_text(size = 30),
                #legend.position=c(0.1, 0.85),
                legend.title=element_text(size=25),
                #legend.background=element_blank(),#legend.background=element_rect(colour=NA, fill=NA),
                legend.key=element_blank(), #legend.key=element_rect(fill=NA),
                #legend.box.background=element_blank(), 
                plot.title=element_text(size=35),           
                panel.border=element_rect(colour="black", fill=NA))
                #guides(shape = guide_legend(override.aes = list(size = 35)))
    p2
    ##    theme_black()
dev.off()




### chromosme Y

#sc[["percent.Y"]] <- PercentageFeatureSet(sc,features=anno[anno$Chr=="chrY",]$kbid)
#sc[["percent.Y"]] %>% summary()



