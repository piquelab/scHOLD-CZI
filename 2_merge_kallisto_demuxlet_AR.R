#
library(Matrix)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(RColorBrewer)
theme_set(theme_grey())
library(furrr)
library(diem) #github install
library(diemr)
library(annotables)
library(data.table)
library("AnnotationHub")

ah <- AnnotationHub()
if(length(ah["AH98047"]) == 0) {
  edb <- ah[["AH75011"]]
} else {
  edb <- ah[["AH98047"]]
}
geneIDs <- genes(edb) %>%
  as.data.frame() %>% 
  setDT(keep.rownames = "ensembl_gene_id") %>%
  .[, c("ensembl_gene_id","entrezid","symbol","seqnames","start","end","strand","gene_biotype", "description")]
names(geneIDs)[c(1,2,4,8)] <- c("ensgene","entrez","chr","biotype")
geneIDs.mt <- subset(geneIDs, chr=="MT")


## Function to read kallisto?
readKallisto  <- function (run, prefixFile="cells_x_genes",expPrefix=NULL) 
{
    if (!dir.exists(paths = run)) {
        stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, paste0(prefixFile,".barcodes.txt"))
    gene.loc <- file.path(run, paste0(prefixFile,".genes.txt"))
##    features.loc <- file.path(run, "features.tsv.gz")
    matrix.loc <- file.path(run, paste0(prefixFile,".nascent.mtx"))
    if (!file.exists(barcode.loc)) {
        stop("Barcode file missing")
    }
    if (!file.exists(gene.loc)) {
        stop("Gene name or features file missing")
    }
    if (!file.exists(matrix.loc)) {
        stop("Expression matrix file missing")
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if(is.null(expPrefix)){
        rownames(data) <- cell.names
    } else {
        rownames(data) <- paste0(expPrefix,"_",cell.names)
    }    
    feature.names <- readLines(gene.loc)
    colnames(data) <- feature.names
    t(data)
}

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_HOLD01-HOLD14_dbgap.ID_cziexp_02_16_2024.txt","CZ1_group.txt") #for testing
base <- args[1]
outFolder=paste0(base,"2b_mergeKallistoAndDemuxlet/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

basefolder=gsub("analysis/","counts_kallisto/nascent/",base)

cov_file=args[2]
# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#read in samples file (just list of samples to run, each sample on newline)
if(!is.na(args[3])){
samples=read.table(args[3],header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
}


## read data including barcodes.txt, genes.txt and mtx
## 2_merge_kb2 for downstream analysis

################################################
### 1, generate folders containing h5ad data ###
################################################

expNames <- dir(basefolder,"^HOLD*")
folders <- paste0(basefolder, expNames, "", sep="")
ind <- dir.exists(folders) #ind <- file.info(folders)$isdir;ind[is.na(ind)]<- FALSE
folders <- folders[ind]
expNames <- expNames[ind]
names(folders) <- expNames
expNames <- names(folders)


###########################################################
### 2, read h5ad data into seurat then merge 39 objects ###
###########################################################

#library(furrr)
#future::plan(strategy = 'multicore', workers = 5)
#options(future.globals.maxSize = 10 * 1024 ^ 3)

## 2.2, read mtx data into seurat object

adata <- future_map(expNames,function(ii){
    ##
    #ii <- expNames[1]
    message(paste0("running ", ii))
    expPrefix = ii;
    cat("#Loading ",paste0(folders[ii], "/counts_unfiltered/"), " ...")    
    sFull <- readKallisto(folders[ii], prefixFile="counts_unfiltered/cells_x_genes")# , expPrefix)  #"/spliced/s"
    cat(dim(sFull),"\n")
    ##
    scs <- colSums(sFull)
    cat(dim(sFull),"\n")
    sel <- sFull
    count0 <- sFull
    cat(dim(sFull),"\n")
    rownames(count0) <- gsub("[SU]-|\\.[0-9]*","",rownames(count0))

    #anno <- merge(data.frame(ensgene=rownames(count0)),geneIDs,by="ensgene",all.x=T)                  
    #metadata <- anno[match(rownames(count0), anno$ensgene),]
    #metadata <- transform(metadata, symbol=ifelse(symbol=="" | is.na(symbol),ensgene,symbol))
    #rownames(count0) <- metadata$symbol

    sc <- CreateSeuratObject(counts = count0, project = "kallisto-CZI1",min.cells = 3, min.features=200)#
    cat(dim(sc),"\n")
    sc@meta.data$Library<-rep(ii,nrow(sc@meta.data))
    #sc
    cat(dim(sc),"\n")
    cat("#Final: ",dim(sc),"\n")
    sc
})


opfn <- paste0(outFolder,"seuratObj-merge-all.",Sys.Date(),".rds") 
write_rds(adata, opfn)

# read the kallisto seurat obj
#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-merge-all."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#adata <- read_rds(opfn)

sc2 <- merge(adata[[1]],adata[-1], add.cell.ids = expNames, project="kall-CZI1")
opfn <- paste0(outFolder,"seuratObj-merge-all-libs-unlist-with-expNames-cell-id.",Sys.Date(),".rds") 
write_rds(sc2, opfn)

##  
###       
sc2@meta.data$NEW_BARCODE <- colnames(sc2)

#anno <- tibble(rn=rownames(sc2)) %>% mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), uns=grepl("U-",rn)) %>% left_join(grch38)
#sc2[["percent.mt"]] <- PercentageFeatureSet(sc2, features = anno %>% dplyr::filter(chr=="MT") %>% dplyr::pull(rn) )

anno <- merge(data.frame(ensgene=rownames(sc2)),geneIDs,by="ensgene",all.x=T)                  
sc2[["percent.mt"]] <- PercentageFeatureSet(sc2, features = rownames(sc2) %in% geneIDs.mt$ensgene  )

opfn <- paste0(outFolder,"1_Seurat_kb.",Sys.Date(),".rds") 
write_rds(sc2, opfn)

### 1_Seurat_kb.rds, unfiltered data and removing diem ## default data, 301,637 barcodes
### 1_Seurat_kb2.rds, unfilered data and removing dime, scs>0 and ucs>0, 304,360 barcodes
### 1_Seurat_kb3.rds, filtered data ## filtered by bustools, 306,299 barcodes 


###########################################
### 2.2, show summary stats of raw data ###
###########################################
###(1)
sc <- sc2

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^1_Seurat_kb."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)

count <- sc@assays$RNA@counts
anno <- data.frame(rn=rownames(count))%>%
        mutate(ensgene=gsub("[SU]-|\\.[0-9]*","",rn), 
               rnz=rowSums(count))

##number of genes
#tmp <- anno%>%filter(uns,rnz>0)

meta <- sc@meta.data
#meta$nCount_spliced <- nCount_spliced
#meta$nFeature_spliced <- nFeature_spliced

#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell=n(),
                       reads=mean(nCount_RNA),
                       ngene=mean(nFeature_RNA),
                       percent.mt=mean(percent.mt),
                       #align="kallisto"
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       .groups="drop")
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))


#### rename the dd

#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell_KL=n(),
                       reads_KL=mean(nCount_RNA),
                       ngene_KL=mean(nFeature_RNA),
                       percent.mt_KL=mean(percent.mt),
                       #align="kallisto"
                       #S_reads=mean(nCount_spliced),
                        #S_ngene=mean(nFeature_spliced), 
                       .groups="drop")
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))

write.csv(as.data.frame(dd), paste0(outFolder,"raw-stats-kallisto-lib.csv"), row.names=F, quote=FALSE)

fig0 <- VlnPlot(sc, features = "percent.mt", ncol = 1, group.by='Library')
png(paste0(figuredir,"Figure0.1_violin_percent_mt.png"), width=4000, height=1000, res=120)
print(fig0)
dev.off()


fig0 <- VlnPlot(sc, features = "nFeature_RNA", ncol = 1, group.by='Library',pt.size = FALSE)
png(paste0(figuredir,"Figure0.2_violin_nfeatures_genes.png"), width=4000, height=1000, res=120)
print(fig0)
dev.off()


fig0 <- VlnPlot(sc, features = "nCount_RNA", ncol = 1, group.by='Library',pt.size = FALSE)+ 
    scale_y_continuous(limits = c(0,20000))
png(paste0(figuredir,"Figure0.3_violin_ncount_umi.png"), width=4000, height=1000, res=120)
print(fig0)
dev.off()


###(1), barcodes for each experiment           
fig0 <- ggplot(dd,aes(x=ident, y=ncell_KL, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("kallisto: #Barcodes per experiment")+
        geom_text(aes(label=ncell_KL),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure1.1_barcodes.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads_KL,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads_KL)
dd2 <- dd%>%dplyr::select(ident,ngene_KL,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene_KL)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="kallisto: #UMIs per cell", "2"="kallisto: #Genes per cell"))
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
fig0 <- ggplot(dd,aes(x=ident, y=reads_KL, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("kallisto: #UMI per cell")+
        geom_text(aes(label=round(reads_KL)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure1.3_UMI_numb.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

fig0 <- ggplot(dd,aes(x=ident, y=ngene_KL, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("kallisto: #Gene per cell")+
        geom_text(aes(label=round(ngene_KL)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure1.4_Gene_numb.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

mean(dd$reads_KL)

fig0 <- ggplot(dd,aes(x=ident, y=percent.mt_KL, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("kallisto: percent mitochondria (mean)")+
        geom_text(aes(label=round(percent.mt_KL)),vjust=-0.7, size=2.)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure2.1_percent_mt.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()


sc[["percent.mt"]] %>% summary()

# prorportion of cells with < xx percent mitochondria for kallisto
mean(sc[["percent.mt"]]<10)
mean(sc[["percent.mt"]]<15)
mean(sc[["percent.mt"]]<20)

sc[["nCount_RNA"]] %>% summary()
sc[["nFeature_RNA"]] %>% summary()

#######################################
### 3, filter data by demux results ###
#######################################

##########################
### 3.1, filtered data ###
##########################
cat("3.1.", "filter data by removing mismatching barcodes", "\n\n")
gc()

#identical(rownames(meta),meta$NEW_BARCODE)
infolder<- paste0(base,"1_demux_output/1_demux_New.SNG.rds")
demux <- read_rds(infolder)

demux <- demux %>% mutate(BATCH=gsub("-.*", "", EXP), treats=gsub(".*[0-9].{,2}-","",EXP)) 
head(demux)

# correct the EtOH ETOH discrepency in treat column 
#demux$treats <- gsub("EtOH", "ETOH", demux$treats)

# filter 
demux <- demux %>% dplyr::filter(NUM.READS>10,NUM.SNPS>10) %>%
  dplyr::select(NEW_BARCODE,NUM.READS,NUM.SNPS,EXP,BATCH,treats,Sample_ID=SNG.BEST.GUESS) 
dim(demux) #647853

matching_bc = intersect(colnames(sc),demux$NEW_BARCODE)

length(matching_bc) # 602137
length(colnames(sc)) # 959956
length(rownames(demux)) # 647853


sc <- subset(sc, cells = matching_bc)

md <- sc@meta.data #%>% rownames_to_column("NEW_BARCODE")
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


###########################################################
### 2.2, show summary stats of post merging with demux ###
###########################################################
###(1)

# load the merged seurat object
#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^seuratObj-merge.md.post-merge-demux."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc2 <- read_rds(opfn)
#sc <- sc2

count <- sc@assays$RNA@counts
anno <- data.frame(ensgene=rownames(count), uns=grepl("S-",rownames(count)),rnz=rowSums(count))

##number of genes
tmp <- anno%>%dplyr::filter(uns,rnz>0)

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
                       .groups="drop")
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))


###(1), barcodes for each experiment           
fig0 <- ggplot(dd,aes(x=ident, y=ncell, fill=factor(batch)))+
        geom_bar(stat="identity")+
        xlab("")+
        scale_y_continuous("", expand=expansion(mult=c(0,0.2)))+
        ggtitle("kallisto: #Barcodes per experiment post filter")+
        geom_text(aes(label=ncell),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure4.1_barcodes_post_filter.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads)
dd2 <- dd%>%dplyr::select(ident,ngene,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="kallisto: #UMIs per cell post filtr", "2"="kallisto: #Genes per cell post filtr"))
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
png(paste0(figuredir,"Figure4.2_genes_post_filter.png"), width=3000, height=2500, res=240)
print(fig0)
dev.off() 


#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell_KL=n(),
                       reads_KL=mean(nCount_RNA),
                       ngene_KL=mean(nFeature_RNA),
                       percent.mt_KL=mean(percent.mt),
                       #align="kallisto"
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       .groups="drop")
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))

write.csv(as.data.frame(dd), paste0(outFolder,"post-merge-kallisto-lib.csv"), row.names=F, quote=FALSE)


# kallisto stats post merge before mt filter
sum(dd$ncell_KL)
median(dd$ncell_KL)
median(dd$reads_KL)
median(dd$ngene_KL)
dim(dd)
dim(sc)
###################### Mitochondria filter
sc[["percent.mt"]] %>% summary()

# prorportion of cells with < xx percent mitochondria
mean(sc[["percent.mt"]]<10) # 0.7819271
mean(sc[["percent.mt"]]<15) # 0.9535661
mean(sc[["percent.mt"]]<20) # 0.9735754


# filter out high mitochondrial percentage
sc <- subset(sc, subset = nFeature_RNA > 200 & percent.mt < 20) #& nFeature_RNA < 20000

dim(sc) #414578

opfn <- paste0(outFolder,"1_Seurat_kb_demux_merged_high_percent_mt_removed.",Sys.Date(),".rds") 
write_rds(sc, opfn)

#opfn_i <- file.info(dir(outFolder, full.names=T, pattern="^1_Seurat_kb_demux_merged_high_percent_mt_removed."))
#opfn <- rownames(opfn_i)[which.max(opfn_i$mtime)]
#sc <- read_rds(opfn)


#############################################################
### 2.2, show summary stats post merge post filtering  ###
#############################################################
###(1)



count <- sc@assays$RNA@counts
anno <- data.frame(rn=rownames(count))%>%
        mutate(rnz=rowSums(count))

##number of genes
#tmp <- anno%>%filter(uns,rnz>0)

meta <- sc@meta.data
#meta$nCount_spliced <- nCount_spliced
#meta$nFeature_spliced <- nFeature_spliced

#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell=n(),
                       reads=mean(nCount_RNA),
                       ngene=mean(nFeature_RNA),
                       percent.mt=mean(percent.mt),
                       #align="kallisto"
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       .groups="drop")
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))


# kallisto stats post merge after mt filter
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
        ggtitle("kallisto: #Barcodes per experiment post filtr")+
        geom_text(aes(label=ncell),vjust=-0.7, size=2.5)+
        theme_bw()+
        theme(legend.title=element_blank(),
              axis.text.x=element_text(angle=90,hjust=1, vjust=0.5, size=7),
              plot.title=element_text(hjust=0.5))  

png(paste0(figuredir,"Figure5.1_barcodes_post_filtr.png"), width=1000, height=600, res=120)
print(fig0)
dev.off()

### (2), reads and number of genes (total, including spliced and unspliced)          
dd1 <- dd%>%dplyr::select(ident,reads,batch)%>%mutate(stats=1)%>%dplyr::rename(y=reads)
dd2 <- dd%>%dplyr::select(ident,ngene,batch)%>%mutate(stats=2)%>%dplyr::rename(y=ngene)
ddnew <- rbind(dd1,dd2)      

stats <- as_labeller(c("1"="kallisto: #UMIs per cell post filtr", "2"="kallisto: #Genes per cell post filtr"))
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
png(paste0(figuredir,"Figure5.2_genes_post_filter.png"), width=3000, height=2500, res=240)
print(fig0)
dev.off() 



#dd <- meta%>%group_by(orig.ident)%>%
dd <- meta%>%group_by(Library)%>%
             summarise(ncell_KL=n(),
                       reads_KL=mean(nCount_RNA),
                       ngene_KL=mean(nFeature_RNA),
                       percent.mt_KL=mean(percent.mt),
                       #align="kallisto"
                       #S_reads=mean(nCount_spliced),
                       #S_ngene=mean(nFeature_spliced), 
                       .groups="drop")
dd <- dd%>%dplyr::rename(ident=Library)%>%
           mutate(batch=gsub("-.*","",ident))

write.csv(as.data.frame(dd), paste0(outFolder,"post-merge-pos-mt-filter-kallisto-lib.csv"), row.names=F, quote=FALSE)

