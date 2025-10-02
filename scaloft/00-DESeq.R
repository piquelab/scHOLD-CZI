library(DESeq2)
library(qvalue)
library(annotables)
library(tidyr)
library(tidyverse)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(plyr);library(dplyr)
library(parallel)
library(ggseurat)
library(sva)

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_uniq-n265_psesl-a2_fixed_12-19-2024.txt","ALL","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt",0.1,50) #old (all treat) was res 0.2 and dim 50

base <- args[1]

cov_file=fread(args[2]) #this is the psych cov file

#if(!is.na(args[3])){
#removeind <- fread(paste0(base,args[3]),header=F)$V1
#project <- sapply(strsplit(args[3],"_"),function(y)y[1])
#}
project=args[3]
method=args[4]
sample_batch <- args[5]
resset <- args[6]
dimset=args[7]
#filter <- "noDEX"
filter <- "CTRLonly" #ALOFT used

#uncoment this section if no longer loading in eigenvec pc file (contains all cov info already)
#cov_file=args[2]
#read in samples file (just list of samples to run, each sample on newline)
#samples=read.table(args[3],header=F)
#samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
#}
#exp <- read.table(cov_file, row.names=NULL,header=T)
#exp <- transform(exp, Sample_ID=dbgap.ID)
#exp$Batch <- gsub("HOLD0","HOLD",exp$Batch)
#if(!is.na(args[3])){
#  exp <- exp %>% dplyr::filter(Batch %in% samples$Batch)
#}
#all treat folders
outFolder=paste0(base,method,"_pseudobulk_ctrl/")
#for resolution 0.2 and including V2 chem
outFolder=paste0(outFolder,"lessfilt/")
outdir=paste0(base,"5b_IdenCelltype_",method,"/")

#ctrl only folders
outFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
outdir=paste0(base,"5b_IdenCelltype_",method,"/",filter,"/")


#outFolder=paste0(outFolder,"dim13res0.3/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
#eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
#eigenvec2 <- merge(eigenvec2_o,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)

#exp=fread(paste0(base,"filtered.incV2.",gsub(".*/","",sample_batch)))
exp <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt")
#eigenvec2 <- merge(exp,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
eigenvec2 <- merge(exp,cov_file,by="dbgap.ID",all.x=T)
#Ali shared that not all waves were run for singlecell and to make sure everything is grabbed correctly this filter is needed
eigenvec2 <- eigenvec2 %>% dplyr::filter(SCAIP7_18 == 1 | SCAIP1_6 == 1)
eigenvec2$Sample_ID <- eigenvec2$dbgap.ID
eigenvec2 <- transform(eigenvec2, SCAIP1_6=ifelse(is.na(SCAIP1_6),0,SCAIP1_6)) #weird case where some are NA ... Ali not sure why
#eigenvec2[eigenvec2$dbgap.ID=="AL-103",c("dbgap.ID","Wave","SCAIP1_6","SCAIP7_18")]

#notrun_var <- c("DSES_01","DSES_03","PWaist","PHip")
#colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
#psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
#psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
#psychvarstorun <- colnames(eigenvec2[,22:57])
psychvarstorun <- colnames(eigenvec2[,24:59])
variables <- c("pedu", "pincme", 
                "psesl", 
                "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasg", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq", 
                "genPC1", "genPC2", "genPC3",  "Sex", "cage1", "ceth1", "cwght1", "chght1", "csex1",
                "cgpd5", "cgpd", "cbpd"
                )
 
variable_names <- c("Parental Education", "Parental Income",
                "Subjective SES", 
                "Neighborhood Stress", "Self-disclosure", 
                "Perceived responsiveness", "Parental Warmth", "YR Parent Child Conflict",  "Risky Family", "Youth Depression", 
                "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Stimulated IL-5 cortisol trt", "Stimulated IL-13 cortisol trt", "Stimulated IFNy cortisol trt",    
                "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", 
                "DD Asthma symptoms", "Nightly Asthma", "Asthma Frequency", "Asthma Severity", 
                "Nightly Awakenings",  "Sleep Duration", "Sleep quality",
                "New Genotype PC1", "New Genotype PC2", "New Genotype PC3","Sex_alph", "Age", "Ethnicity", "Weight", "Height", "sex",
                "Female Menarche Status", "Female Puberty Score", "Male puberty score"
                )
variables_df <- data.frame(variable=variables, description=variable_names)

opfn <- paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
sc <- read_rds(opfn)

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
geneIDs.sex <- subset(geneIDs, chr=="Y" | chr=="X")
geneIDs.male <- subset(geneIDs, chr=="Y" )
geneIDs.female <- subset(geneIDs, chr=="X")

#making cov file for eqtl mapping
#First row gives the sample ID and each additional one corresponds to a single covariate
#First column gives the covariate ID and each additional one corresponds to a sample
#this selects the batchwave seen in the experimental data
u_eigenvec2 <- merge(sc@meta.data[,c("Sample_ID","BATCH")],unique(eigenvec2[,-2]),by.x=c("Sample_ID","BATCH"),by.y=c("Sample_ID","Batch2"))
u_eigenvec2 <- u_eigenvec2 %>%
  dplyr::select(Sample_ID, everything())
fwrite(u_eigenvec2, sep='\t', quote=F, row.names=F, col.names=T, file=paste0("/rs/rs_grp_scaloft/scALOFT_2024/covariates/",filter,".eigenvec2_u.txt"))
#u_eigenvec2[u_eigenvec2$Sample_ID=="AL-103",c("Sample_ID","Wave","SCAIP1_6","SCAIP7_18")]

#from https://www.biostars.org/p/9482789/
sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
colOrd <- data.frame(NEW_BARCODE=names(sc$orig.ident))
coldata_ex <- merge(sc@meta.data,unique(eigenvec2[,-2]),by.x=c("Sample_ID","BATCH"),by.y=c("Sample_ID","Batch2")) #first 2 columns are dbgapid ie. sampleid and batch
RcolData <- merge(colOrd, coldata_ex, by="NEW_BARCODE",all.x=T)
RcolData <- RcolData[match(rownames(sc@meta.data), RcolData$NEW_BARCODE), ]
sc@meta.data <-cbind(sc@meta.data,RcolData[,which(!colnames(RcolData) %in% colnames(sc@meta.data))])
sc@meta.data$orig.ident <- "scaloft_comb"
#check wavefilt worked:
#unique(sc@meta.data[sc@meta.data$Sample_ID=="AL-030","Wave"]) #should be the data from wave B1
#unique(sc@meta.data[sc@meta.data$Sample_ID=="AL-103",c("Sample_ID","Wave","SCAIP1_6","SCAIP7_18")])
fwrite(sc@meta.data, file=paste0(outFolder,"scmetadata_allind.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

sc_genesdf <- ldply(lapply(unique(sc@meta.data$letter_clusters),function(c){
    cat("running",c,"\n")
    sc_c <- subset(sc, subset=letter_clusters==c)
    sc_genes <- data.frame(cluster=c,genes=rownames(sc[["RNA"]]))
    rm(sc_c)
    gc()
    return(sc_genes)
}),data.frame)
    fwrite(sc_genesdf, file=paste0(outFolder,"scmetadata_allgenes.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#subset for 20cell count filt (just cell/ind)
nind <- ddply(sc@meta.data, c("letter_clusters"), plyr::summarize,
  ind=length(unique(Sample_ID)))
fwrite(nind, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".indcount.txt"))
counts <- plyr::count(sc@meta.data,c("letter_clusters","Sample_ID"))
names(counts)[3] <- "ic_cellcounts"
sc@meta.data <- left_join(sc@meta.data,counts,by=c("letter_clusters","Sample_ID"))
rownames(sc@meta.data) <- sc@meta.data$NEW_BARCODE
before=dim(sc@meta.data)[1]
sc <- subset(sc, subset=ic_cellcounts>=20)
after=dim(sc@meta.data)[1]
cat("removed ", before-after, "combos","\n") #35,523
#head(rownames(sc@meta.data))
table(sc@meta.data$letter_clusters,sc@meta.data$treats)
       CTRL   LPS LPS-DEX   PHA PHA-DEX
  C0  48528 49956    6627 27166   29898
  C1  29820 30360    4342 14351   16765
  C10  2772  2881     346  2361    2149
  C11  1577  1444      93  1671    1518
  C12    50    68       5    46      53
  C13    11     2       0     7       4
  C2  25334 26909    3266 16662   16584
  C24    37    43      42    23      40
  C3  21629 22272    2150 16297   16545
  C4  17942 18712    2036 17283   17773
  C5  15454 16300    2256 15861   15321
  C6   1327  1598      88 32345   23901
  C7  14236  7781    2195  6922   11776
  C8   6550  6495     888  4402    4627
  C9   3829  3961     355  3792    3867

nind <- ddply(sc@meta.data, c("letter_clusters"), plyr::summarize,
  ind=length(unique(Sample_ID)))
nind_ctrl <- subset(nind, treats=="CTRL")
fwrite(nind, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".indcount_icfilt.txt"))

sc[["sex_male"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.male$symbol])
sc[["sex_female"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.female$symbol])

opfn <- paste0(outFolder,project,".",resset,".",dimset,".wavefilt.icfilt.seurat.RDS")
saveRDS(sc, file=opfn)
#sc <- read_rds(opfn)

#subset for 20cell count filt (combo without batch)
nind <- ddply(sc@meta.data, c("letter_clusters","treats"), plyr::summarize,
  ind=length(unique(Sample_ID)))
fwrite(nind, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".indcount.txt"))
counts <- plyr::count(sc@meta.data,c("treats","letter_clusters","Sample_ID"))
names(counts)[4] <- "tic_cellcounts"
sc@meta.data <- left_join(sc@meta.data,counts,by=c("treats","letter_clusters","Sample_ID"))
rownames(sc@meta.data) <- sc@meta.data$NEW_BARCODE
before=dim(sc@meta.data)[1]
sc <- subset(sc, subset=tic_cellcounts>=20)
after=dim(sc@meta.data)[1]
cat("removed ", before-after, "combos","\n") #35,523
#head(rownames(sc@meta.data))
table(sc@meta.data$letter_clusters,sc@meta.data$treats)
       CTRL   LPS LPS-DEX   PHA PHA-DEX
  C0  48479 49898    6627 27096   29833
  C1  29760 30268    4331 14102   16627
  C10  1476  1414      21   834     837
  C11   493   347       0   734     455
  C2  25260 26818    3247 16529   16414
  C3  21456 22100    1933 16076   16189
  C4  17613 18311    1872 16909   17394
  C5  15100 15930    2065 15499   14920
  C6    595   584       0 32002   23515
  C7  13729  7093    2050  6228   11306
  C8   5451  5528     438  3311    3481
  C9   2668  2873      60  2702    2794

nind <- ddply(sc@meta.data, c("letter_clusters","treats"), plyr::summarize,
  ind=length(unique(Sample_ID)))
nind_ctrl <- subset(nind, treats=="CTRL")
fwrite(nind, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".indcount_ticfilt.txt"))

sc[["sex_male"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.male$symbol])
sc[["sex_female"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.female$symbol])

opfn <- paste0(outFolder,project,".",resset,".",dimset,".wavefilt.ticfilt.seurat.RDS")
saveRDS(sc, file=opfn)

#subset for 20cell count filt (combo with batch)
counts <- plyr::count(sc@meta.data,c("Library","letter_clusters","Sample_ID"))
names(counts)[4] <- "btic_cellcounts"
sc@meta.data <- left_join(sc@meta.data,counts,by=c("Library","letter_clusters","Sample_ID"))
rownames(sc@meta.data) <- sc@meta.data$NEW_BARCODE
before=dim(sc@meta.data)[1]
sc <- subset(sc, subset=btic_cellcounts>=20)
after=dim(sc@meta.data)[1]
cat("removed ", before-after, "combos","\n") #35,523
#head(rownames(sc@meta.data))
table(sc@meta.data$letter_clusters,sc@meta.data$treats)
       CTRL   LPS LPS-DEX   PHA PHA-DEX
  C0  48479 49898    6627 27096   29833
  C1  29760 30268    4331 14102   16627
  C10  1476  1414      21   834     837
  C11   493   347       0   734     455
  C2  25260 26818    3247 16529   16414
  C3  21456 22100    1933 16076   16189
  C4  17613 18311    1872 16909   17394
  C5  15100 15930    2065 15499   14920
  C6    595   584       0 32002   23515
  C7  13729  7093    2050  6228   11306
  C8   5451  5528     438  3311    3481
  C9   2668  2873      60  2702    2794

nind <- ddply(sc@meta.data, c("letter_clusters","treats"), plyr::summarize,
  ind=length(unique(Sample_ID)))
nind_ctrl <- subset(nind, treats=="CTRL")
   letter_clusters treats ind
1               C0   CTRL 217
6               C1   CTRL 216
11             C10   CTRL  48
16             C11   CTRL  19
20              C2   CTRL 216
25              C3   CTRL 209
30              C4   CTRL 195
35              C5   CTRL 191
40              C6   CTRL  16
44              C7   CTRL 154
49              C8   CTRL 120
54              C9   CTRL  78
fwrite(nind, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".indcount_bticfilt.txt"))

#if(!is.na(args[3])){
#sc@meta.data$keep <- ifelse(sc@meta.data$Sample_ID %in% removeind, FALSE, TRUE)
#sc <- subset(sc, subset = keep == TRUE)
#sc$keep <- NULL
#}

sc[["sex_male"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.male$symbol])
sc[["sex_female"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.female$symbol])
#sc@meta.data <- transform(sc@meta.data, SCAIP1_6=ifelse(is.na(SCAIP1_6),0,SCAIP1_6)) #weird case where some are NA ... Ali not sure why

fwrite(sc@meta.data, file=paste0(outFolder,"scmetadata_bticfilt.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

sc_genesdf <- ldply(lapply(unique(sc@meta.data$letter_clusters),function(c){
    cat("running",c,"\n")
    sc_c <- subset(sc, subset=letter_clusters==c)
    sc_genes <- data.frame(cluster=c,genes=rownames(sc[["RNA"]]))
    rm(sc_c)
    gc()
    return(sc_genes)
}),data.frame)
    fwrite(sc_genesdf, file=paste0(outFolder,"scmetadata_bticfiltgenes.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

opfn <- paste0(outFolder,project,".",resset,".",dimset,".wavefilt.bticfilt.seurat.RDS")
saveRDS(sc, file=opfn)

#sc <-readRDS(opfn)
var="pedu"
sc@meta.data$SCAIP1_6 <- as.factor(sc@meta.data$SCAIP1_6)
for (var in psychvarstorun[c(1:10)]){
  p <- ggplot(data = sc@meta.data) +
    geom_boxplot(aes(x=SCAIP1_6, y=get(var),colour=SCAIP1_6)) 
png(paste0(figuredir,project,".box_scaip1vs2_",var,".png"), width=1000, height=1000, res=120)
print(p)
dev.off()
}

###############  WIP
sharedcols <- Reduce(intersect, list(colnames(eigenvec2),colnames(sc@meta.data)))
eigen_shared <- eigenvec2[,..sharedcols]
sc_shared <- sc@meta.data[,sharedcols]
identical(eigenvec2[,..sharedcols],as.data.table(sc@meta.data)[,..sharedcols])
identical_cols <- sapply(names(eigen_shared), function(x) identical(eigen_shared[[x]], sc_shared[[x]]))


##################

# 1 is male and 0 is female
sc@meta.data$Sex_c <- ifelse(sc@meta.data$Sex==0, "female","male")
p <- ggplot(data = sc@meta.data) +
  xlab("sex by questionaire") +
  cowplot::theme_cowplot()
p1 <- p + geom_boxplot(aes(x=Sex_c, y=sex_male,colour=Sex_c)) + ylab("sex by male chromosome")
p2 <- p + geom_boxplot(aes(x=Sex_c, y=sex_female,colour=Sex_c)) + ylab("sex by male chromosome")
fig <-plot_grid(p1, p2, nrow=1, ncol=2, align="h")
png(paste0(figuredir,project,".box_ncount_sex_averageind.png"), width=3000, height=1000, res=120)
print(fig)
dev.off()


#additional plots (Sex, cage1, Wave, ceth1,factor_age,genPC1,pincme,cgpd,cbpd)
sc@meta.data$integer_age=as.integer(sc@meta.data$cage1)
sc@meta.data$factor_age=cut.default(sc@meta.data$cage1, 3, labels = NULL, include.lowest = FALSE, right = TRUE, dig.lab = 3)

umapfiguredir="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/2.1_mergeCellRangerAnddemux/figures/"
plotgroup="genPC1"
fname=paste0(umapfiguredir,project,".Figure5.1_UMAP_Harmony-res",resset,".",dimset,"_group_",plotgroup,".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- DimPlot(sc, reduction = "umap", group.by = plotgroup, pt.size = .5)+
        ggtitle("")+
  theme(legend.key.size = unit(50,"point"),panel.background = element_rect(fill="white",colour = "black"),
    legend.title=element_blank(),axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)), legend.text=element_text(size = rel(1.8)))
print(fig1)
dev.off()

plotgroup="cgpd"
fname=paste0(umapfiguredir,project,".Figure5.1_UMAP_Harmony-res",resset,".",dimset,"_group_",plotgroup,".png");
png(fname,width=5000,height=5000, res=240)
fig1 <- FeaturePlot(sc, features = plotgroup,  reduction = "umap", 
                        cols = c("lightgrey", "darkblue")) & 
  theme(panel.background = element_rect(fill="white",colour = "black"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),
    axis.text.y = element_text(colour = "black",size = rel(1.3)), axis.title.x = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.3)))
print(fig1)
dev.off()


#findmarkers keep all clusters >10 cells
cellcount <- as.data.frame(table(sc@meta.data$letter_clusters, sc@meta.data$orig.ident))
lowcell <- subset(cellcount, Freq<10)
if(dim(lowcell)[1]==0){cat( "no low cell counts")}else{sc <-subset(x = sc, subset = letter_clusters %in% unique(lowcell$Var1), invert = TRUE)}

library(MAST)
lapply(sc@meta.data$letter_clusters, function(type){
    if(!isTRUE(file.size(paste0(outFolder,project,".",resset,".",dimset,".",type,"vsAll",".findmarkers.txt")) > 0)){
  message(paste0("running ",type,"vsAll"))
  df <- sc
  Idents(df) <- df$letter_clusters
cluster.markers <- FindMarkers(object = df, ident.1 = type, min.pct = 0.25)
  cluster.markers <-transform(cluster.markers, gene=rownames(cluster.markers))
    fwrite(cluster.markers, sep='\t', quote=F, row.names=F, file=paste0(outFolder,project,".",resset,".",dimset,".",type,"vsAll",".findmarkers.txt"))
  rm(df)
  gc(reset=TRUE)
    }
})

#get number of cells per individual
cellcount <- as.data.frame(table(sc@meta.data$Sample_ID, sc@meta.data$orig.ident))
png(width = 8, height = 8, file=paste0(figuredir,project,".",resset,".",dimset,".cellcounts_perindividual_hist.png"), pointsize=12, 
  bg = "transparent", units = "in", res = 1200)
ggplot(data=cellcount, aes(Freq)) + 
theme_bw()+
geom_histogram(bins=50)
dev.off()

length(which(cellcount$Freq<100))

#this was abandonded
#lapply(sc@meta.data$letter_clusters, function(type){
#  treat="PHA-DEX"
#  control="PHA"
#    if(!isTRUE(file.size(paste0(outFolder,project,".",resset,".",dimset,".",type,".",treat,"vs",control,".findmarkers.txt")) > 0)){
#  message(paste0("running ",treat,"vs",control," ",type))
#  df <- subset(sc, subset = letter_clusters == type)
#  Idents(df) <- df$treats
#  #if(length(subset(df, subset=treats=="CTRL", return.null = TRUE)$Barcode)>5 & length(subset(df, subset=treats=="LPS", return.null = TRUE)$Barcode)>5) {
#  condition_diff1 <- FindMarkers(df, ident.1 = treat, ident.2 = control, logfc.threshold = 0.1, min.pct = .05, test.use = "MAST", latent.vars = c('percent.mt','nFeature_RNA','Sex','cage1'))
#  fwrite(condition_diff1, sep='\t', quote=F, row.names=F, file=paste0(outFolder,project,".",resset,".",dimset,".",type,".",treat,"vs",control,".findmarkers.txt"))
#  rm(df)
#  gc(reset=TRUE)
#    }
#})


cellcount <- as.data.frame(table(sc@meta.data$letter_clusters, sc@meta.data$orig.ident))
lowcell <- subset(cellcount, Freq<3000) #had it at 3k, checking if I can use 1k -- 1k not good - 1sample/batch error running combat
#basded on this should remove C11 as well -- holding off for now
if(dim(lowcell)[1]==0){cat( "no low cell counts")}else{sc <-subset(x = sc, subset = letter_clusters %in% unique(lowcell$Var1), invert = TRUE)}
# [1] C12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24

cellcount <- transform(cellcount, removed=ifelse(Var1 %in% unique(lowcell$Var1), "Y", "N"))
fwrite(cellcount[order(cellcount$Freq),-2], sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".cellcount_bticfilt.txt"))

fig0 <- VlnPlot(sc, features = c("sex_male", "sex_female"), ncol = 2,pt.size = FALSE) 
png(paste0(figuredir,project,".violin_ncount_sex.png"), width=3000, height=1000, res=120)
print(fig0)
dev.off()

mean_sex <- ddply(sc@meta.data, "Sample_ID", plyr::summarize,
    sex_male=mean(sex_male, na.rm=T),
    sex_female=mean(sex_female, na.rm=T),
    Sex=unique(Sex))

p <- ggplot(data = mean_sex) +
  cowplot::theme_cowplot()
p1 <- p + geom_boxplot(aes(x=Sex, y=sex_male,colour=Sex))
p2 <- p + geom_boxplot(aes(x=Sex, y=sex_female,colour=Sex))
fig <-plot_grid(p1, p2, nrow=2, ncol=1, align="h")
png(paste0(figuredir,project,".box_ncount_sex_averageind.png"), width=3000, height=1000, res=120)
print(fig)
dev.off()

sce <- as.SingleCellExperiment(sc)
#seurat_clusters, treats, BATCH, Library
#opfn <- paste0(outFolder,project,".",resset,".",dimset,".wavefilt.icfilt.singlecellexp.RData")
#save(sce, file=opfn)
rm(sc)
gc()
#load(opfn)

library(biomaRt)  
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",1:22), mart = mart)
genesplusxy <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",c(1:22,"X","Y")), mart = mart)

#/ aggregate by cluster,library info and covariate of interest:
#sce <- subset(sce, ,!letter_clusters %in% unique(lowcell$Var1))

sum_by <- c("letter_clusters", "BATCH","Sample_ID","treats")
summed <- aggregateAcrossCells(sce, id=colData(sce)[,sum_by])
# Number of cells per sample and cluster
t <- table(colData(sce)$Sample_ID,
           colData(sce)$letter_clusters)
metadata <- as.data.frame(colData(sce))

#/ add rownames using the information from the colData:
colnames(summed) <- apply(colData(summed)[,sum_by], 1, function(x) paste(x, collapse="_"))
#head(assay(summed, "counts"))
raw <- assay(summed, "counts")

rm(sce)
gc()
counts_ls_XY <- lapply(unique(summed$letter_clusters), function(i){
    cat("running ",i,"\t")
    #i <- unique(summed$letter_clusters)[1]
    cell_idx <- which(tstrsplit(colnames(summed), "_")[[1]] == i)
    keep <- which(colnames(raw) %in% colnames(summed[, cell_idx]))
    data <- raw[, keep]
    filtered_data <- data > 0 
    #filtered_data[filtered_data < 0] <- NA
    cat(dim(data),"\t")
    if(is.matrix(filtered_data)){
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(data) / 4) #filter to keep genes expressed in 25% of samples
    data <- unlist(data[keep, ])
    cat(dim(data),"\t")
    data <- data[rownames(data) %in% genesplusxy$hgnc_symbol, ]
    cat(dim(data),"\n")
    summed_filt <- summed[rownames(summed) %in% rownames(data), cell_idx]
    #table(rownames(summed_filt) %in% geneIDs.sex$symbol)
    return(summed_filt)
    } else {
    NULL
    }
  #counts_ls[[i]] <- summed[, cell_idx]
  #names(counts_ls) <-
    })
names(counts_ls_XY) <-unique(summed$letter_clusters)
counts_ls_XY[sapply(counts_ls_XY, is.null)] <- NULL

metadata_ls_XY <- lapply(counts_ls_XY, function(i){
    #i <- counts_ls[[1]]
    cat("running ",unique(i$letter_clusters),"\t")
    df <- data.frame(cluster_sample_id = colnames(i)) ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$letter_clusters <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$BATCH  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    df$Sample_ID  <- tstrsplit(df$cluster_sample_id, "_")[[3]]
    df$treats  <- tstrsplit(df$cluster_sample_id, "_")[[4]]
    test <- metadata %>% dplyr::select(-c(sex_male,sex_female,orig.ident,NEW_BARCODE,percent.mt,nCount_RNA,nFeature_RNA,NUM.READS,NUM.SNPS))
    #need to remove all RNA_snn_res columns so one per each resolution
    test <- dplyr::select(test, -contains("RNA_snn_res"))
    #test <- test %>% dplyr::select(-c(RNA_snn_res.0.3,RNA_snn_res.0.4))
    test <- unique(transform(test,cluster_sample_id=paste(letter_clusters,BATCH,Sample_ID,treats,sep="_")))
    #test[test$cluster_sample_id=="C0_SCAIP10_AL-021_CTRL",c(1:20)]
    #test <- test[!duplicated(test$cluster_sample_id),]
    df_n <- merge(df, unique(test), all.x=T)
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df_n) <- df_n$cluster_sample_id
    #df_n <- df_n[complete.cases(df_n), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    return(df_n)
})

all(names(counts_ls_XY) == names(metadata_ls_XY))

opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.plusXY.RData")
save(counts_ls_XY,metadata_ls_XY, file=opfn)

counts_ls <- lapply(unique(summed$letter_clusters), function(i){
    cat("running ",i,"\t")
    #i <- unique(summed$letter_clusters)[1]
    cell_idx <- which(tstrsplit(colnames(summed), "_")[[1]] == i)
    keep <- which(colnames(raw) %in% colnames(summed[, cell_idx]))
    data <- raw[, keep]
    filtered_data <- data > 0 
    #filtered_data[filtered_data < 0] <- NA
    cat(dim(data),"\t")
    if(is.matrix(filtered_data)){
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(data) / 4) #filter to keep genes expressed in 25% of samples
    data <- unlist(data[keep, ])
    cat(dim(data),"\t")
    data <- data[rownames(data) %in% genes$hgnc_symbol, ]
    cat(dim(data),"\n")
    summed_filt <- summed[rownames(summed) %in% rownames(data), cell_idx]
    #table(rownames(summed_filt) %in% geneIDs.sex$symbol)
    return(summed_filt)
    } else {
    NULL
    }
  #counts_ls[[i]] <- summed[, cell_idx]
  #names(counts_ls) <-
    })
names(counts_ls) <-unique(summed$letter_clusters)
counts_ls[sapply(counts_ls, is.null)] <- NULL

metadata_ls <- lapply(counts_ls, function(i){
    #i <- counts_ls[[1]]
    cat("running ",unique(i$letter_clusters),"\t")
    df <- data.frame(cluster_sample_id = colnames(i)) ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$letter_clusters <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$BATCH  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    df$Sample_ID  <- tstrsplit(df$cluster_sample_id, "_")[[3]]
    df$treats  <- tstrsplit(df$cluster_sample_id, "_")[[4]]
    test <- metadata %>% dplyr::select(-c(sex_male,sex_female,orig.ident,NEW_BARCODE,percent.mt,nCount_RNA,nFeature_RNA,NUM.READS,NUM.SNPS))
    #need to remove all RNA_snn_res columns so one per each resolution
    test <- dplyr::select(test, -contains("RNA_snn_res"))
    #test <- test %>% dplyr::select(-c(RNA_snn_res.0.3,RNA_snn_res.0.4))
    test <- unique(transform(test,cluster_sample_id=paste(letter_clusters,BATCH,Sample_ID,treats,sep="_")))
    #test[test$cluster_sample_id=="C0_SCAIP10_AL-021_CTRL",c(1:20)]
    #test <- test[!duplicated(test$cluster_sample_id),]
    df_n <- merge(df, unique(test), all.x=T)
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df_n) <- df_n$cluster_sample_id
    #df_n <- df_n[complete.cases(df_n), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    return(df_n)
})

all(names(counts_ls) == names(metadata_ls))

opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
save(counts_ls,metadata_ls, file=opfn)

load(opfn)

clusters <- names(counts_ls)
fwrite(clusters, sep='\t', quote=F, row.names=F, col.names=F, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))

cluster_celltype <- ldply(lapply(metadata_ls,function(i){
    df <- unique(data.frame(cluster=i$letter_clusters,cell_type=i$customclassif))
    return(df)
}),data.frame)[,-1]
fwrite(cluster_celltype, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".cluster_celltype.txt"))

#combatrun="income_PCs_sex_age_and_treats_adjusted"
combatrun="income_PCs_sex_age_adjusted" #CTRL only doesnt need treatment
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
if (!file.exists(baseoutFolder)) dir.create(baseoutFolder, showWarnings=F)
if (!file.exists(paste0(baseoutFolder,"figures/"))) dir.create(paste0(baseoutFolder,"figures/"), showWarnings=F)

mclapply(names(counts_ls), function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    #highcell <- cluster_metadata_bf[cluster_metadata_bf$cell_count>=20,c("comb","letter_clusters","cell_count")]
    #cluster_metadata <- cluster_metadata_bf[rownames(cluster_metadata_bf) %in% rownames(highcell),]
    #cluster_metadata_var <- unique(cluster_metadata[,c("Sample_ID","BATCH","treats","Sex","cage1","genPC1","genPC2","genPC3","pincme")])
    cluster_metadata_var <- unique(cluster_metadata[,c("Sample_ID","BATCH","Sex","cage1","genPC1","genPC2","genPC3","pincme")]) #CTRL only doesnt need treat
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    x <- subset(as.data.frame(table(cluster_metadata_var$BATCH)),Freq>1)
    cluster_metadata_var <- subset(cluster_metadata_var,BATCH %in% unique(x$Var1))

    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
    #adjusted <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","Sex","cage1","genPC1","genPC2","genPC3","pincme")])
    adjusted <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("Sex","cage1","genPC1","genPC2","genPC3","pincme")])#CTRL only doesnt need treat

    #cat(length(rownames(adjusted)),"\n")    
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    save(adjusted, file=opfn)
})

genestoremove <- ldply(lapply(names(counts_ls), function(cluster){
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
max=unique(rownames(which(adjusted >= .Machine$integer.max, arr.ind = TRUE)))
if(is.null(max)){
    max=NA
}
return(data.frame(clus=cluster,genes=max))
}),data.frame)
genestoremove #35
#none so move on

#treatment results lps v ctrl
#combatrun="income_PCs_sex_age_and_treats_adjusted"
combatrun="income_PCs_sex_age_adjusted" #CTRL only doesnt need treatment
#run="income_PCs_sex_age_and_treats_COMBAT_treatment"
run="income_PCs_sex_age_COMBAT_treatment" #CTRL only doesnt need treatment
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

contrastdf_base <- data.frame(control=c("CTRL","LPS","CTRL","PHA"),treatment=c("LPS","LPS.DEX","PHA","PHA.DEX"))
contrastdf <- transform(contrastdf_base, contrast1=paste0(control),contrast2=paste0("treats_",treatment))
contrastdf <- transform(contrastdf,control=gsub("[.]","-",control),treatment=gsub("[.]","-",treatment))

for (cluster in names(counts_ls)){
  #[c(7:9)] #to run clusters if job terminates part way through
  #cluster="C6"
      cat("running ", cluster, "\n")
  fdr=0.05
  #cluster_counts_sce <- counts_ls[[cluster]]
  cluster_metadata_sce <- metadata_ls[[cluster]]
  #cluster_counts <- assay(cluster_counts_sce, "counts")
  cluster_metadata <- data.frame(cluster_metadata_sce)
  cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
  opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
  load(opfn)
  lapply(1:nrow(contrastdf),function(x) {
    con=contrastdf[x,]
    contrast=paste0(con$contrast2,"_vs_",con$contrast1)
    if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",run,".txt")) > 0)){
    cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
    cluster_metadata_var <- within(cluster_metadata_var, treats <- factor(relevel(treats, ref = con$control)))
    cluster_metadata_var <- cluster_metadata_var[,c("Sample_ID","treats")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
    cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts)),]
    all(colnames(cluster_counts) == rownames(cluster_metadata_var))

    cat("running",contrast,"deseq \n")
    # Create DESeq2 object        
    dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cluster_metadata_var, 
                                  design = ~ Sample_ID + treats)
    dds <- DESeq(dds,parallel=TRUE)
    opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",contrast,".",cluster,".",run,".RData")
    save(dds, file=opfn)

    res <- results(dds, name =contrast)
    sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
    names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
    sub.table <- sub.table[!is.na(sub.table$padj), ]
    cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
    sub.table$cluster=cluster
    sub.table$contrast=contrast
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",contrast,".",run,".txt"))
    sigDEGs <- subset(sub.table,padj<fdr)
    sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
    fwrite(sigDEGs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",contrast,".",run,".txt"))
    table <- data.frame(cluster=cluster,contrast=contrast)
    table$number_samples <- paste(nrow(cluster_metadata_var))
    table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts))
    table$tested_genes <- paste(nrow(sub.table))
    table$DEGs_FDR <- paste(nrow(sigDEGs))
    table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",run,".txt"))
    }
    })
}

png(width = 6, height = 6, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",cluster,".",contrast,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
plotMA(res)
dev.off()


for (x in c(1:nrow(contrastdf))){
    con=contrastdf[x,]
    contrast=paste0(con$contrast2,"_vs_",con$contrast1)
    cat("running ",contrast,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",contrast,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",contrast,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",contrast,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",contrast,".",run,".txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".sigDEGs.",contrast,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",contrast,".",run,".txt"))
}

library(flextable)
library(ftExtra)
library(rlist)
library(officer)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)

#[!psychvarstorun %in% puberty]
subvars <- ldply(lapply(c(1:nrow(contrastdf)), function(x){
    con=contrastdf[x,]
    contrast=paste0(con$contrast2,"_vs_",con$contrast1)
    cat("running ",contrast,"\n")
    if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",contrast,".",run,".txt"))){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",contrast,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"),select=c(1,6))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newls <- newls[sapply(newls, nrow)>0]
            newtable <- data.frame(contrast=contrast,list.cbind(newls))
            return(newtable)
        }
    }
}),data.frame)

dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=max)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("numind|ngenes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() 
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = 1, border = border, part = "all")
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(subsubvars_degonly)), align = "center", part = "body")

#dft <- align(dft, i = 1, j = NULL, align = "left", part = "body")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".ndegonly.png"))


filenames <- list.files(outFolder) #file list from directory
filenames1 <- filenames[grep("DESeq_output-", filenames)] #pick specific files from list
x <- na.omit(sapply(strsplit(filenames1,"(?<=.)(?=C[0-9])",perl=TRUE),function(y)y[2]))
clusters <- unique(sapply(strsplit(x,"[.]"),function(z)z[1]))

#for (cluster in clusters){
cluster=clusters[1]
    cat("running ", cluster, "\n")
    var="LPSvsCTRL"
    load(paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-treats",cluster,".RData"))

    vsd <- vst(dds, blind=FALSE)
    # Plot PCA
    fig0 <- DESeq2::plotPCA(vsd, intgroup = c("treats"))+ggtitle("treats")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    #fig1 <- DESeq2::plotPCA(vsd, intgroup = c("PC1"))+ggtitle("PC1")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    #fig2 <- DESeq2::plotPCA(vsd, intgroup = c("PC2"))+ggtitle("PC2")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig3 <- DESeq2::plotPCA(vsd, intgroup = c("Sex"))+ggtitle("sex")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig4 <- DESeq2::plotPCA(vsd, intgroup = c("cage1"))+ggtitle("age")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig <-plot_grid(fig0,fig3,fig4, nrow=3, ncol=1, align="h")
    png(paste0(figuredir,project,".",resset,".",dimset,".deseqPCA-treats","-",cluster,".png"), width=1000, height=1000, res=120)
    print(fig)
    dev.off()
#}

combatrun="income_PCs_sex_age_adjusted"
run="income_PCs_sex_age_adjusted_withWave"
#combatrun="income_PCs_sex_age_and_treats_adjusted"
#run="income_PCs_sex_age_and_treats_adjusted_withWave_poscount"
#run="income_PCs_sex_age_and_treats_adjusted_withWave"
#combatrun="PCs_sex_age_and_treats_adjusted_wavefilt"
#run="PCs_sex_age_and_treats_adjusted_withWave"
#combatrun="SES_PCs_sex_age_and_treats_adjusted"
#run="SES_PCs_sex_age_and_treats_adjusted_withWave"
#run="sex_age_and_treats_adjusted_wavefilt"
#run="SES_PCs_sex_age_and_treats_adjusted"
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/",filter,"/")
combatfolder=paste0(baseoutFolder,combatrun,"/")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleidwave/"))) dir.create(paste0(outFolder,"sampleidwave/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid_all/"))) dir.create(paste0(outFolder,"sampleid_all/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

runage=FALSE
withWave=TRUE
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
firstrunvars <- psychvarstorun[c(1:10)]
#firstrunvars <- "psesl"
for (cluster in names(counts_ls)){
    #cluster="C0"
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    for(i in unique(cluster_metadata_sce$treats)){
            #i <- "CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, Sex=as.factor(Sex))
        if(all(c("0","1") %in% levels(cluster_metadata_t$Sex))){
            cluster_metadata_t <- within(cluster_metadata_t, Sex <- relevel(Sex, ref = "0"))
        }
        if(runage==TRUE){
        if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-age-",i,".",run,".txt")) > 0)){
            var="age"
            #d <- data.frame(number_ind=length(unique(cluster_metadata[,c("Sample_ID")])))
            #fwrite(d, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sampleid_all.",run,".txt"))
            #d <- data.frame(number_ind=length(unique(cluster_metadata_t[,c("Sample_ID")])),cluster=cluster,treats=i)
            #fwrite(d, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sampleid_all","-",i,".",run,".txt"))

            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1")]
            cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
            #fwrite(cluster_metadata_var[,c("Sample_ID","Wave")], sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sampleidwave_",var,"-",i,".",run,".txt"))
            cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
            cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
            all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
            if(withWave==TRUE){
            design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + cage1")
            } else {
            design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + cage1")
            }
            cat("running deseq ",var,cluster,i," \n")
            dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                          colData = cluster_metadata_var, 
                                          design = as.formula(design))
            dds <- DESeq(dds,parallel=TRUE)
            opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RDS")
            saveRDS(dds, file=opfn)

            res <- results(dds)
            sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
            names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
            sub.table <- sub.table[!is.na(sub.table$padj), ]
            cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
            sub.table$var=var
            sub.table$cluster=cluster
            sub.table$treats =i
            fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))
            sigDEGs <- subset(sub.table,padj<fdr)
            sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
            fwrite(sigDEGs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"))
            table <- data.frame(symb=var, variable= var, cluster=cluster)
            table$number_samples <- paste(nrow(cluster_metadata_var))
            table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
            table$gene_number <- paste(nrow(cluster_counts_t))
            table$DEGs_FDR <- paste(nrow(sigDEGs))
            table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
        }
        }
        mclapply(firstrunvars[c(1:length(firstrunvars))][!firstrunvars[c(1:length(firstrunvars))] %in% puberty],function(var){
        #var="pnsi"
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1",var)]
                cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
                fwrite(cluster_metadata_var[,c("Sample_ID","Wave")], sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".sampleidwave_",var,"-",i,".",run,".txt"))
                cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
                cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
                all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
                if(withWave==TRUE){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + Sex + cage1 + ",var)
                } else {
                design <-  paste0("~ genPC1 + genPC2 + genPC3 + Sex + cage1 + ",var)
                }
                cat("running deseq ",var,cluster,i," \n")
                dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                              colData = cluster_metadata_var, 
                                              design = as.formula(design))
                dds <- estimateSizeFactors(dds, type = 'poscounts')
                #design = as.formula(paste0("~ LibBatch + Wave + Sex + cage1 + SCAIP1_6_genPC1 + SCAIP1_6_genPC2 + SCAIP1_6_genPC3 + as.numeric(",myvar, ")")))
                dds <- DESeq(dds,parallel=TRUE)
                opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RDS")
                saveRDS(dds, file=opfn)

                res <- results(dds)
                sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
                names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
                sub.table <- sub.table[!is.na(sub.table$padj), ]
                cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
                sub.table$var=var
                sub.table$cluster=cluster
                sub.table$treats =i
                fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))
                sigDEGs <- subset(sub.table,padj<fdr)
                sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
                fwrite(sigDEGs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"))
                table <- data.frame(symb=var, variable= var, cluster=cluster)
                table$number_samples <- paste(nrow(cluster_metadata_var))
                table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
                table$gene_number <- paste(nrow(cluster_counts_t))
                table$tested_genes <- paste(nrow(sub.table))
                table$DEGs_FDR <- paste(nrow(sigDEGs))
                table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
                fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
            }
        })
    }
}

#this was some QC work for checking individuals 
#var="pnsi"
#cluster="C6"
#i="CTRL"
#opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RDS")
#dds <- read_rds(opfn)
#res <- results(dds)
#metadata(res)$filterThreshold

#png(width = 8, height = 8, file=paste0(figuredir,project,".filterrejections-",var,"-",cluster,"-",i,".",run,".png"), pointsize=12, 
#  bg = "transparent", units = "in", res = 1200)
#plot(metadata(res)$filterNumRej, 
#     type="b", ylab="number of rejections",
#     xlab="quantiles of filter")
#lines(metadata(res)$lo.fit, col="red")
#abline(v=metadata(res)$filterTheta)
#dev.off()

#m <- melt(cluster_counts_t)
#png(width = 8, height = 8, file=paste0(figuredir,project,".countsxind-",var,"-",cluster,"-",i,".",run,".png"), pointsize=12, 
#  bg = "transparent", units = "in", res = 1200)
#ggplot(m,aes(x=Var1,y=value))+
#  theme_bw()+
#  geom_point()+ #aes(color=sig)
#    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
#    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
#    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
#    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
#dev.off()

#png(width = 8, height = 8, file=paste0(figuredir,project,".counts_hist-",var,"-",cluster,"-",i,".",run,".png"), pointsize=12, 
#  bg = "transparent", units = "in", res = 1200)
#ggplot(data=m, aes(value)) + 
#theme_bw()+
#geom_histogram(bins=50)
#dev.off()

#var="pedu"
#cluster="C6"
#i="CTRL"
#m <- melt(cluster_counts_t)
#png(width = 8, height = 8, file=paste0(figuredir,project,".countsxind-",var,"-",cluster,"-",i,".",run,".png"), pointsize=12, 
#  bg = "transparent", units = "in", res = 1200)
#ggplot(m,aes(x=Var1,y=value))+
#  theme_bw()+
#  geom_point()+ #aes(color=sig)
#    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
#    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
#    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
#    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
#dev.off()

#png(width = 8, height = 8, file=paste0(figuredir,project,".counts_hist-",var,"-",cluster,"-",i,".",run,".png"), pointsize=12, 
#  bg = "transparent", units = "in", res = 1200)
#ggplot(data=m, aes(value)) + 
#theme_bw()+
#geom_histogram(bins=50)
#dev.off()


sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sampleid_all/",project,".",resset,".",dimset,".",cluster,".sampleid_all.",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sampleid_all/",project,".",resset,".",dimset,".",cluster,".sampleid_all.",run,".txt"),header=T)
        }
        sub.table <- transform(sub.table,cluster=cluster)
        }),data.frame)
fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid_all/",project,".",resset,".",dimset,".sampleid_all_combtreat.",run,".txt"))

sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        ldply(lapply(c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX"), function(i){
        sub.table <- fread(paste0(outFolder,"sampleid_all/",project,".",resset,".",dimset,".",cluster,".sampleid_all","-",i,".",run,".txt"),header=T)
        }),data.frame)
        }),data.frame)
fwrite(reshape2::dcast(sub.table, cluster ~ treats,value.var="number_ind"), sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid_all/",project,".",resset,".",dimset,".sampleid_all.",run,".txt"))

for (var in c("age",psychvarstorun)[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]){
    #for (var in c(firstrunvars)){
    for (i in c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,".sampleidwave_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,".sampleidwave_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".sampleidwave_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".sigDEGs_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    }
}
#filter: ALL.stats_all_cell_types*.nomissvars_PCs*






#puberty variables run separately as they run per each sex
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")[c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds") %in% colnames(cov_file)]
for (cluster in names(counts_ls)){
#for (cluster in c("C8","C9")){
    #cluster="C0"
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "CTRL"))
    cluster_metadata <- transform(cluster_metadata, Sex=as.factor(ifelse(Sex==0,"female","male")))
    if(all(c("male","female") %in% levels(cluster_metadata$Sex))){
        cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "male"))
    }
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    #opfn <- paste0(combatfolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
            #run covariates separately for each treatment condition
    for(i in unique(cluster_metadata_sce$treats)){
    #for(i in c("LPS","LPS-DEX","PHA","PHA-DEX")){
    #i <- "CTRL"
        for(s in c("male","female")){
        if (s=="female"){
        puberty_sex <- c("cgpd5","cgpd")
            }else {
        puberty_sex <- c("cbpd")
            }
        cluster_metadata_t <- subset(cluster_metadata, treats==i & Sex==s)
            mclapply(puberty_sex, function(var){
        #var="pnsi"
        missing=all(is.na(cluster_metadata_t[,var]))
            if(isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,"_",s,".stats_all_ce+ _types-",var,"-",i,".",run,".txt")) >0)){
                cat(var,cluster,i,s, "has run already \n")
                next
                } else if (missing){
                    cat("var is all missing information")
                    next
                    } else{
                cat("running")
                cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3",var)]
                cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
                fwrite(cluster_metadata_var[,c("Sample_ID","Wave")], sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,"_",s,".sampleidwave_",var,"-",i,".",run,".txt"))
                cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
                cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
                all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
                if(withWave==TRUE){
                design <-  paste0("~ Wave + genPC1 + genPC2 + genPC3 + ",var)
                } else {
                design <-  paste0("~ genPC1 + genPC2 + genPC3 + ",var)
                }
                cat("running deseq ",var,cluster,i,s," \n")
                dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                              colData = cluster_metadata_var, 
                                              design = as.formula(design))
                dds <- DESeq(dds,parallel=TRUE)
                opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,"_",s,".",run,".RDS")
                saveRDS(dds, file=opfn)

                res <- results(dds)
                sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
                names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
                sub.table <- sub.table[!is.na(sub.table$padj), ]
                cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
                sub.table$var=var
                sub.table$cluster=cluster
                sub.table$treats =i
                sub.table$sex =s
                fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,"_",s,".deseqres_",var,"-",i,".",run,".txt"))
                sigDEGs <- subset(sub.table,padj<fdr)
                sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
                fwrite(sigDEGs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,"_",s,".sigDEGs_",var,"-",i,".",run,".txt"))
                table <- data.frame(symb=var, variable= var, cluster=cluster,sex=s)
                table$number_samples <- paste(nrow(cluster_metadata_var))
                table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
                table$gene_number <- paste(nrow(cluster_counts_t))
                table$tested_genes <- paste(nrow(sub.table))
                table$DEGs_FDR <- paste(nrow(sigDEGs))
                table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
                fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,"_",s,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
                }
            })
        }
    }
}

for (var in puberty){
    #for (var in c(firstrunvars[c(5:10)],"age")){
    for (i in c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")){
        for (s in c("male","female")){
        cat("running ",var," ",i,s,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,"_",s,".sampleidwave_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".",cluster,"_",s,".sampleidwave_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleidwave/",project,".",resset,".",dimset,".sampleidwave_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,"_",s,".deseqres_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,"_",s,".deseqres_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,"_",s,".sigDEGs_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,"_",s,".sigDEGs_",var,"-",i,".",run,".txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".sigDEGs_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,"_",s,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,"_",s,".stats_all_cell_types-",var,"-",i,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    }
}
}


#automated table making
#- top is cluster; sub top is num ind. | num genes | num DEG
#- side is variable
#one per each treatment
# assign variabeles
variables <- c("pedu", "pincme", 
                "psesl", 
                "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasg", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq", 
                "genPC1", "genPC2", "genPC3",  "Sex", "cage1", "ceth1", "cwght1", "chght1", "csex1",
                "cgpd5", "cgpd", "cbpd"
                )
 
variable_names <- c("Parental Education", "Parental Income",
                "Subjective SES", 
                "Neighborhood Stress", "Self-disclosure", 
                "Perceived responsiveness", "Parental Warmth", "YR Parent Child Conflict",  "Risky Family", "Youth Depression", 
                "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Stimulated IL-5 cortisol trt", "Stimulated IL-13 cortisol trt", "Stimulated IFNy cortisol trt",    
                "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", 
                "DD Asthma symptoms", "Nightly Asthma", "Asthma Frequency", "Asthma Severity", 
                "Nightly Awakenings",  "Sleep Duration", "Sleep quality",
                "New Genotype PC1", "New Genotype PC2", "New Genotype PC3","Sex_alph", "Age", "Ethnicity", "Weight", "Height", "sex",
                "Female Menarche Status", "Female Puberty Score", "Male puberty score"
                )
variables_df <- data.frame(variable=variables, description=variable_names)

library(flextable)
library(ftExtra)
library(rlist)
library(officer)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)

#[!psychvarstorun %in% puberty]
for (i in c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")){
subvars <- ldply(lapply(psychvarstorun, function(var){
#subvars <- ldply(lapply(firstrunvars, function(var){
    cat("running",i,var,"\n")
    if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"),select=c(1,7))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,list.cbind(newls))
            return(newtable)
        }
    }
}),data.frame)

dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=max)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("numind|ngenes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() 
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = 2, border = border, part = "all")
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(3:length(subsubvars_degonly)), align = "center", part = "body")

#dft <- align(dft, i = 1, j = NULL, align = "left", part = "body")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".ndegonly.png"))


}


#total DEGs
run="income_PCs_sex_age_and_treats_adjusted_withWave"
outFolder=paste0(outFolder,run,"/")

subvars <- ldply(lapply(c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX"), function(i){
  psy <- ldply(lapply(c(variables_df$variable,"age"), function(var){
#psy <- ldply(lapply(variables_df$variable, function(var){
    cat("running",i,var,"\n")
    if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    stats=transform(stats, treat=i)
    deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    DEGs_FDR_10=unique(subset(deseqres,padj<0.1)[,c("identifier","treats","var")])
    return(DEGs_FDR_10)
    }
    }),data.frame)
  return(psy)
  }), data.frame)

ctrl <- subset(subvars, treats=="CTRL")
ctrl_s <- ddply(ctrl, c("var"), plyr::summarize,
  sigDEGs=length(unique(identifier)))
fwrite(ctrl_s, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-CTRL.",run,".txt"))

stim <- subset(subvars, treats=="LPS" | treats=="PHA")
stim_s <- ddply(stim, c("var"), plyr::summarize,
  sigDEGs=length(unique(identifier)))
fwrite(stim_s, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-stim.",run,".txt"))

dex <- subset(subvars, treats=="LPS-DEX" | treats=="PHA-DEX")
dex_s <- ddply(dex, c("var"), plyr::summarize,
    sigDEGs=length(unique(identifier)))
fwrite(dex_s, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-dex.",run,".txt"))
