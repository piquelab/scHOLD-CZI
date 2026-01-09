require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
library(annotables)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(data.table)
library(plyr); library(dplyr)
library(parallel)
library("AnnotationHub")
library(ggseurat)
library(cowplot)
library(sva)
library(ggpubr)
future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 15 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.15) #for testing
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.1) #for testing
#old args
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",11,0.2) #for testing

base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
outFoldertemp=paste0(base,method,"_pseudobulk_ctrl/")
if (!file.exists(outFoldertemp)) dir.create(outFoldertemp, showWarnings=F)
# set new output dir for filtered out unmatched figures
outFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
#outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)
leadvar <- fread("/rs/rs_grp_schold/covariates/other_covariates/HOLD LEAD 5.27.25.csv")
cov_pluslead <- merge(cov_file,leadvar,by="pID",all=T)
cov_pluslead <- transform(cov_pluslead, Lead=ifelse(Lead==-99,NA,Lead))
#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)

eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o[,-c("sex","sex_alph","age")],cov_pluslead,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip","SNI_NoP","age","sex","sex_alph","isel","pID") #SNI_NoP is the only variable that should be excluded based on the observed issues in score distributions that don’t make sense (negative values and extreme outliers)
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
old_cytokines <- psychvarstorun[c(8,14:24)] #old cytokines
psychvarstorun <- psychvarstorun[!psychvarstorun %in% old_cytokines]
firstrunvars=c("SES","pr_comp","ISEL_Mean","PSS_all_mean","BPd_avg","Chol_HDL","Chol_LDL","nii_mean","SNI_NumPeople_r","BPs_avg","DED_all_mean","LogCRP","HVS_mean","LivingAlone")
reordered_psychvarstorun <- psychvarstorun[c(which(psychvarstorun %in% firstrunvars),which(!psychvarstorun %in% firstrunvars))]
secondrunvars <- c(reordered_psychvarstorun[c(1:20)],"age","Lead")
zcytokines <- c("cytocomp","z_ifny_0_log_w",   
"z_il10_0_log_w","z_il12p70_0_log_w",
"z_il13_0_log_w",   
"z_il1b_0_log_w","z_il2_0_log_w",    
"z_il4_0_log_w",    
"z_il6_0_log_w","z_il8_0_log_w",    
"z_tnfa_0_log_w")
treatments=c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")
treatmentsfirst=c("RNA-CTRL","RNA-LPS")
allPFAS <- psychvarstorun[c(34:54)]
combPFAS <- allPFAS[seq(1,21,3)]
PFASvars <- c(allPFAS,"sumPFAS","log_sumPFAS")
allvars <- c(reordered_psychvarstorun[c(1:83)],"age","Lead","sumPFAS","log_sumPFAS")

# Variables to focus on
    variables <- c("DSES_07","DSES_09","FCDEM_12","chronic_sum","smoke","CVDRISK",
                "HS_CRP","il6","PSQI_total","HVS_mean","EDS_mean","SES","LogCRP","Logifny","Logil10","Logil12",
                "Logil13","Logil1b","Logil2","Logil4","Logil6","Logil8","Logtnfa","cytocomp","NAIscr2019","NAIscr2017",
                "LQ2019","StressSev","StressCount","nii_mean","NII_fam","SLS","LivingAlone","SNI_HCG","SNI_NoP","PSS_all_mean",
                "DED_all_mean","Chol","Trig","Chol_HDL","Chol_Ratio","Chol_LDL","HbA1C","BPs_avg","BPd_avg","ISEL_Mean", "pr_comp",
                "SNI_NumPeople_r","LogDED","LogTrig","chronic_sum_categ",
                "age","factor_HS_CRP"
                )
variable_names <- c("Participant's highest level of education","Pre-tax household income","Are you currently taking any prescription medications",
                    "Total number of chronic conditions","Smoking status","Cardiovascular disease risk",
                    "High Sensitivity CRP (mg/L) raw values","IL-6 raw values","Total sleep score (PSQI)","Vigilance Scale - Mean","Everyday Discrimination Scale - Mean",
                    "SES composite [mean(ZDSES_07, ZDSES_09)]","LogCRP","Logifny","Logil10","Logil12",
                    "Logil13","Logil1b","Logil2","Logil4","Logil6",
                    "Logil8","Logtnfa","Cytokines composite [mean(ZLogtnfa, ZLogi6, ZLogifny)]","Neighborhood Adversity Score","Neighborhood Adversity Score",
                    "Neighborhood Segregation Score","Core: Total Severity of Stressors","Core: Total Count of Stressors","Negative interactions mean","Negative interactions with family","Avg on loneliness scale",
                    "Dichotomous variable measuring living alone (1 = living alone)","High contact groups: number of categories of social groups p interacted with at least once every two weeks",
                    "This is the mean score representing the number of people with whom the respondent has regular contact","Mean of daily pss measures","Daily experience of discrimination cumulative average",
                    "Cholesterol in mg/dL","Triglycerides in mg/dL","HDL Cholesterol in mg/dL","CHOL/HDL ratio","LDL Cholesterol in mg/dL",
                    "Glycated HGB Affinity HPLC HbA1C","Average systolic blood pressure (wonky BP3_s measure corrected)","Average diastolic blood pressure","social support","psychological resources composite",
                    "number of pople in contact with","Log of Daily Discrimination", "LogTrig", "Chronic conditions transformed",
                    "age","High Sensitivity CRP (mg/L) factorized 1 to 5"
                    )
variables_df <- data.frame(variable=variables, description=variable_names)
additionalvars <- fread("/rs/rs_grp_schold/covariates/other_covariates/chronic_chond_lables_attributes_varnames.txt")
variables_df <- merge(variables_df,additionalvars[!additionalvars$Variable %in% c("HS_CRP"),-1],by.x=c("variable","description"),by.y=c("Variable","Column_Label"),all=T)
allvarsdf <- data.frame(variable=c(allvars,PFASvars,"Lead"))
variables_df <- merge(variables_df,allvarsdf,by="variable",all=T)

cat_cov <- fread(paste0(base,"categories_cov.txt")) 

cat_cov_m <- merge(variables_df,cat_cov[,-2],by.y="Column_Name",by.x="variable",all.y=T)

#Lead histogram
myvar <- "Lead"
keep <- as.vector(!is.na(eigenvec2[, ..myvar]))
cov <- eigenvec2[keep,]
cov$variable <- myvar
var_name=variables_df[variables_df$variable==myvar,]
# Function to insert newlines into long titles
wrap_title <- function(title, width = 55) {
  paste(strwrap(title, width = width), collapse = "\n")
}

p1 <- ggplot(cov, aes(x=get(myvar)))+#, y=get(myvar))) +
    geom_histogram(color = "black", fill = "dodgerblue4", alpha=0.5) +
    ggtitle(wrap_title(paste0(var_name), width = 55)) +  # Apply the wrapping function here
    theme(plot.title = element_text(hjust = 0.5, size=20), 
          axis.text.y = element_text(size = 20), 
          axis.text.x = element_text(size = 20), 
          axis.title = element_text(size = 20),
          plot.subtitle = element_text(size = 18, hjust = 0.5)) +
    #theme(legend.position="none", axis.title.x = element_blank()) +
    labs(x=paste0(myvar), subtitle = paste0("n=",nrow(cov))) 
figfn <- paste(base,method,"_pseudobulk_ctrl/nodex/figures/", myvar,"_histogram_distribution",".png",sep="")
png(figfn, width=1000, height=1000, res=120)
#print(plot_grid(fig1, fig3, ncol=2, labels="AUTO", label_fontface="plain"))
print(p1)
dev.off()

options(future.globals.maxSize = 60 * 1024 ^ 3)

## RPR AR changing from annotation hub to annotables back. 
##ah <- AnnotationHub()
##if(length(ah["AH98047"]) == 0) {
##  edb <- ah[["AH75011"]]
##} else {
##  edb <- ah[["AH98047"]]
##}

##geneIDs <- genes(edb) %>%
##  as.data.frame() %>% 
##  setDT(keep.rownames = "ensembl_gene_id") %>%
##  .[, c("ensembl_gene_id","entrezid","symbol","seqnames","start","end","strand","gene_biotype", "description")]
##names(geneIDs)[c(1,2,4,8)] <- c("ensgene","entrez","chr","biotype")
geneIDs <- grch38 ## %>% select(ensgene,entrez,symbol,chr,start,end,bio)
geneIDs.sex <- subset(geneIDs, chr=="Y" | chr=="X")
geneIDs.male <- subset(geneIDs, chr=="Y" )
geneIDs.female <- subset(geneIDs, chr=="X")

outdir=paste0(base,"5b_IdenCelltype_",method,"/noDEX/")
#outdir=paste0(base,"5b_IdenCelltype_",method,"/")
opfn <- paste0(outdir,project,".seuratObj-.harmony-sctype-",resset,".",dimset,".rds")
sc <- read_rds(opfn)
#from https://www.biostars.org/p/9482789/
#note that this object has old isel not iselmean
sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
colOrd <- data.frame(NEW_BARCODE=names(sc$orig.ident))
coldata_ex <- merge(as.data.frame(sc@meta.data),eigenvec2,by="Sample_ID") #first 2 columns are dbgapid ie. sampleid and batch
RcolData <- merge(colOrd, coldata_ex, by="NEW_BARCODE",all.x=T)
RcolData <- RcolData[match(rownames(sc@meta.data), RcolData$NEW_BARCODE), ]
sc@meta.data <-cbind(sc@meta.data,RcolData[,which(!colnames(RcolData) %in% colnames(sc@meta.data))])
fwrite(sc@meta.data, file=paste0(outFolder,"scmetadata_allind.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#sc_genesdf <- ldply(lapply(unique(sc@meta.data$letter_clusters),function(c){
#    cat("running",c,"\n")
#    sc_c <- subset(sc, subset=letter_clusters==c)
#    sc_genes <- data.frame(cluster=c,genes=rownames(sc[["RNA"]]))
#    rm(sc_c)
#    gc()
#    return(sc_genes)
#}),data.frame)
#    fwrite(sc_genesdf, file=paste0(outFolder,"scmetadata_allgenes.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#subset for 20cell count filt
counts <- plyr::count(sc@meta.data,c("Library","letter_clusters","Sample_ID"))
names(counts)[4] <- "btic_cellcounts"
sc@meta.data <- left_join(sc@meta.data,counts,by=c("Library","letter_clusters","Sample_ID"))
rownames(sc@meta.data) <- sc@meta.data$NEW_BARCODE
before=dim(sc@meta.data)[1]
sc <- subset(sc, subset=btic_cellcounts>=20)
after=dim(sc@meta.data)[1]
cat("removed ", before-after, "combos","\n") #removed 3797 combos
table(sc@meta.data$letter_clusters,sc@meta.data$treats)
#  RNA-CTRL RNA-LPS
#  C0   134245  139961
#  C1    49547   55992
#  C2    31138   35991
#  C3    37146   28506
#  C4    20917   23281
#  C5      684     811
#  C6       72       0

# Nov 2025 re-run
#     RNA-CTRL RNA-LPS
#  C0   132000  145679
#  C1    52141   61731
#  C2    31857   35901
#  C3    37207   28699
#  C4    21185   23934
#  C5      668     847
#  C6      108     131


cellcount <- as.data.frame(table(sc@meta.data$letter_clusters, sc@meta.data$orig.ident))
fwrite(cellcount, file=paste0(outFolder,"cellcount_cluster.txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

lowcell <- subset(cellcount, Freq<1000) #was 3k, going to see if 1k works
if(dim(lowcell)[1]==0){cat( "no low cell counts")}else{sc <-subset(x = sc, subset = letter_clusters %in% unique(lowcell$Var1), invert = TRUE)}

#sce <- as.SingleCellExperiment(sc) ### GetAssayData doesn't work for multiple layers in v5 assay.
# a slution is convert v5 assay to a v3 assay with the following code. It seems like this is basically trying to merge the layers.
#sc <- sc2
sc[["RNA"]] <- as(sc[["RNA"]], Class="Assay")

sc[["sex_male"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.male$symbol])
sc[["sex_female"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.female$symbol])

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

opfn <- paste0(outFolder,project,".",resset,".",dimset,".bticfilt.seurat.RDS")
saveRDS(sc, file=opfn)

#sc <- read_rds(opfn)

p <- ggplot(data = sc) +
  cowplot::theme_cowplot()
p1 <- p + geom_boxplot(aes(x=sex_alph, y=sex_male,colour=sex_alph))
p2 <- p + geom_boxplot(aes(x=sex_alph, y=sex_female,colour=sex_alph))
fig <-plot_grid(p1, p2, nrow=2, ncol=1, align="h")
png(paste0(figuredir,project,".box_ncount_sex.png"), width=3000, height=1000, res=120)
print(fig)
dev.off()

fig0 <- VlnPlot(sc, features = c("sex_male", "sex_female"), ncol = 2,pt.size = FALSE) 
png(paste0(figuredir,project,".violin_ncount_sex.png"), width=3000, height=1000, res=120)
print(fig0)
dev.off()

mean_sex <- ddply(sc@meta.data, "Sample_ID", plyr::summarize,
    sex_male=mean(sex_male, na.rm=T),
    sex_female=mean(sex_female, na.rm=T),
    sex_alph=unique(sex_alph))

p <- ggplot(data = mean_sex) +
  cowplot::theme_cowplot()
p1 <- p + geom_boxplot(aes(x=sex_alph, y=sex_male,colour=sex_alph))
p2 <- p + geom_boxplot(aes(x=sex_alph, y=sex_female,colour=sex_alph))
fig <-plot_grid(p1, p2, nrow=2, ncol=1, align="h")
png(paste0(figuredir,project,".box_ncount_sex_averageind.png"), width=3000, height=1000, res=120)
print(fig)
dev.off()

# save thhe seurat obkect post filtering out clusters with less than 1k, and ag
opfn <- paste0(outFolder,project,".",resset,".",dimset,".seuratObj-all-layers-prior-to-sce-conversion",Sys.Date(),".rds") 
write_rds(sc, opfn)

sce <- as.SingleCellExperiment(sc)

#opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".SingleCellExperiment.RDS")
#write_rds(sce, opfn)
rm(sc)
rm(sc2)

gc()

#seurat_clusters, treats, BATCH, Library
#/ aggregate by cluster,library info and covariate of interest:
sum_by <- c("letter_clusters", "BATCH","Sample_ID","treats")
summed <- aggregateAcrossCells(sce, id=colData(sce)[,sum_by])

#/ add rownames using the information from the colData:
colnames(summed) <- apply(colData(summed)[,sum_by], 1, function(x) paste(x, collapse="_"))
raw <- assay(summed, "counts")

library(biomaRt)  
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",1:22), mart = mart)

counts_ls <- lapply(unique(summed$letter_clusters), function(i){
    cat("running ",i,"\n")
    #i <- unique(summed$letter_clusters)[1]
    cell_idx <- which(tstrsplit(colnames(summed), "_")[[1]] == i)
    keep <- which(colnames(raw) %in% colnames(summed[, cell_idx]))
    data <- raw[, keep]
    filtered_data <- data > 0 
    #filtered_data[filtered_data < 0] <- NA
    if(is.matrix(filtered_data)){
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(data) / 4) #filter to keep genes expressed in 25% of samples
    data <- unlist(data[keep, ])
    data <- data[rownames(data) %in% genes$hgnc_symbol, ]
    summed_filt <- summed[rownames(summed) %in% rownames(data), cell_idx]
    return(summed_filt)
    } else {
    NULL
    }
    })
names(counts_ls) <-unique(summed$letter_clusters)
counts_ls[sapply(counts_ls, is.null)] <- NULL

# Number of cells per sample and cluster
t <- table(colData(sce)$Sample_ID,
           colData(sce)$letter_clusters)

metadata_ls <- lapply(counts_ls, function(i){
    #i <- counts_ls[[1]]
        cat("running ",unique(i$letter_clusters),"\n")
    df <- data.frame(cluster_sample_id = colnames(i)) ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$letter_clusters <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$BATCH  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    df$Sample_ID  <- tstrsplit(df$cluster_sample_id, "_")[[3]]
    df$treats  <- tstrsplit(df$cluster_sample_id, "_")[[4]]
    metadata <- as.data.frame(colData(sce))
    test <- metadata %>% dplyr::select(-c(sex_male,sex_female,orig.ident,NEW_BARCODE,percent.mt,nCount_RNA,nFeature_RNA,NUM.READS,NUM.SNPS))
    test <- dplyr::select(test, -contains("RNA_snn_res"))
    df_n1 <- merge(df, unique(test), all.x=T)
    rownames(df_n1) <- df_n1$cluster_sample_id
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    #df_n <- df_n[complete.cases(df_n), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    return(df_n1)
})
all(names(counts_ls) == names(metadata_ls))

#opfn <- paste0(outFolder,project,".DESeq_countlists.RData")
opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt_proteincoding.RData")
save(counts_ls,metadata_ls, file=opfn)

## RPR AR
## WE stopped here because we do not want to use combatset or anything from combatseq. 

counts_ls <- lapply(unique(summed$letter_clusters), function(i){
#counts_ls <- lapply(unique(summed$letter_clusters)[!unique(summed$letter_clusters) %in% "C5"], function(i){
    cat("running ",i,"\n")
    #i <- unique(summed$letter_clusters)[1]
    cell_idx <- which(tstrsplit(colnames(summed), "_")[[1]] == i)
    keep <- which(colnames(raw) %in% colnames(summed[, cell_idx]))
    data <- raw[, keep]

    celltype_combo <- data.frame(cluster.withdex=c("C0","C2","C3","C4","C5","C6","C1"),celltype.withdex=c("CD4_T","CD8_T","NK","Monocyte","B","DC","CD4_T"),
        cluster.withoutdex=c("C0","C1","C2","C3","C4","C5","C6"),celltype.withoutdex=c("CD4_T","CD8_T","NK","Monocyte","B","DC","Monocyte"))
    celltype_comboc <- subset(celltype_combo,cluster.withoutdex==i )
    combatrun="SES_PCs_sex_age_and_treats_adjusted"
    baseoutFolder_dex=paste0(base,method,"_pseudobulk_ctrl/adjusted/",0.2,".",11,"/cell20filt/")
    opfn <- paste0(baseoutFolder_dex,project,".",0.2,".",11,".ComBat_seq.",celltype_comboc$cluster.withdex,".",combatrun,".RData")
    load(opfn)
    sc_genes_oldmatch <- data[rownames(data) %in% rownames(adjusted_counts), ]
    all(rownames(sc_genes_oldmatch) == rownames(adjusted_counts))
    filtered_data <- sc_genes_oldmatch > 0 
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(sc_genes_oldmatch) / 4) #filter to keep genes expressed in 25% of samples
    sc_genes_oldmatch_25 <- unlist(sc_genes_oldmatch[keep, ]) #lose AL596223.1
    #rownames(sc_genes_oldmatch)[!rownames(sc_genes_oldmatch) %in% rownames(sc_genes_oldmatch_25)]
    sc_genes_oldmatch_25 <- sc_genes_oldmatch_25[!rownames(sc_genes_oldmatch_25) %in% geneIDs.sex$symbol, ]
    summed_filt_oldmatch <- summed[rownames(summed) %in% rownames(sc_genes_oldmatch_25), cell_idx]
    fwrite(data.frame(cluster=i,genes=rownames(summed_filt_oldmatch)), file=paste0(outFolder,"pseudobulk_genes_oldmatch",i,".txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

    filtered_data <- data > 0 
    #filtered_data[filtered_data < 0] <- NA
    if(is.matrix(filtered_data)){
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(data) / 4) #filter to keep genes expressed in 25% of samples
    data <- unlist(data[keep, ])
    data <- data[!rownames(data) %in% geneIDs.sex$symbol, ]
    summed_filt <- summed[rownames(summed) %in% rownames(data), cell_idx]
    diffgenes <- rownames(summed_filt)[!rownames(summed_filt) %in% rownames(sc_genes_oldmatch_25)] #1161
    fwrite(data.frame(cluster=i,genes=rownames(summed_filt)), file=paste0(outFolder,"pseudobulk_genes",i,".txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    fwrite(data.frame(cluster=i,genes=diffgenes), file=paste0(outFolder,"pseudobulk_genes_diffto_olddex",i,".txt"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    return(summed_filt)
#    return(summed_filt_oldmatch)
    } else {
    NULL
    }
    })
names(counts_ls) <-unique(summed$letter_clusters)
#names(counts_ls) <-unique(summed$letter_clusters)[!unique(summed$letter_clusters) %in% "C5"]
counts_ls[sapply(counts_ls, is.null)] <- NULL

# Number of cells per sample and cluster
t <- table(colData(sce)$Sample_ID,
           colData(sce)$letter_clusters)

metadata_ls <- lapply(counts_ls, function(i){
    #i <- counts_ls[[1]]
        cat("running ",unique(i$letter_clusters),"\n")
    df <- data.frame(cluster_sample_id = colnames(i)) ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$letter_clusters <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$BATCH  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    df$Sample_ID  <- tstrsplit(df$cluster_sample_id, "_")[[3]]
    df$treats  <- tstrsplit(df$cluster_sample_id, "_")[[4]]
    metadata <- as.data.frame(colData(sce))
    test <- metadata %>% dplyr::select(-c(sex_male,sex_female,orig.ident,NEW_BARCODE,percent.mt,nCount_RNA,nFeature_RNA,NUM.READS,NUM.SNPS))
    test <- dplyr::select(test, -contains("RNA_snn_res"))
    df_n1 <- merge(df, unique(test), all.x=T)
    rownames(df_n1) <- df_n1$cluster_sample_id
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    #df_n <- df_n[complete.cases(df_n), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    return(df_n1)
})
all(names(counts_ls) == names(metadata_ls))

#opfn <- paste0(outFolder,project,".DESeq_countlists.RData")
opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt_oldmatch.RData")

save(counts_ls,metadata_ls, file=opfn)

#this checks what combos are affected by the filter
lcf <- ldply(lapply(names(counts_ls), function(cluster){
  cluster_metadata_sce <- metadata_ls[[cluster]]
  cluster_metadata <- data.frame(cluster_metadata_sce)
  lowcell <- cluster_metadata[cluster_metadata$cell_count<20,c("comb","letter_clusters","cell_count")]
  return(lowcell)
  }), data.frame)
table(lcf$letter_clusters)

#metadata_ls <- lapply(metadata_ls, function(i){
#    ii <- merge(i, unique(eigenvec2[,c("Sample_ID","ISEL_Mean")]), by=c("Sample_ID"),all=T)
#    return(ii)
#    })
#opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.RData")
#save(counts_ls,metadata_ls, file=opfn)

################################################################################################3
##################################################################################################
opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
load(opfn)

#only using SES var, also accounting for PC1 and PC2
combatrun="SES_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
if (!file.exists(baseoutFolder)) dir.create(baseoutFolder, showWarnings=F)

baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/test/")
if (!file.exists(baseoutFolder)) dir.create(baseoutFolder, showWarnings=F)


remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 3 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
}


mclapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
    cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata_bf <- data.frame(cluster_metadata_sce)
    #highcell <- cluster_metadata_bf[cluster_metadata_bf$cell_count>=20,c("comb","letter_clusters","cell_count")]
    #cluster_metadata <- cluster_metadata_bf[rownames(cluster_metadata_bf) %in% rownames(highcell),]
    cluster_metadata_var <- cluster_metadata_bf[,c("Sample_ID","BATCH","PC1","PC2","treats","sex_alph","age","SES")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

#checkrm <- remove_outliers(cluster_counts_t)
##checkrm2 <- rm.outlier(cluster_counts_t, fill = FALSE, median = FALSE, opposite = FALSE)
#check <- outlier(cluster_counts_t) #"HBB"    "MALAT1"
#f <- apply(check, 1, any)
#outlier_genes <- rownames(cluster_counts_t)[f]
#    cluster_counts_t <- cluster_counts_t[!rownames(cluster_counts_t) %in% outlier_genes,]

#            genesdexvnodex <- fread(paste0("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/SES_PCs_sex_age_and_treats_generem/genes_dexvsnodex.",var,".txt"))
#            bothgenesdexvnodex <- subset(genesdexvnodex, gene_present=="both")
#            cluster_counts_t <- cluster_counts_t[rownames(cluster_counts_t) %in% bothgenesdexvnodex$identifier, ]


    adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age","PC1","PC2","SES")])
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    save(adjusted_counts, file=opfn)
})

genestoremove <- ldply(lapply(names(counts_ls), function(cluster){
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
max=unique(rownames(which(adjusted_counts >= .Machine$integer.max, arr.ind = TRUE)))
return(data.frame(genes=max))
}),data.frame)
genestoremove #35

combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
mclapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
    cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    #highcell <- cluster_metadata_bf[cluster_metadata_bf$cell_count>=20,c("comb","letter_clusters","cell_count")] #duplicate filter
    #cluster_metadata <- cluster_metadata_bf[rownames(cluster_metadata_bf) %in% rownames(highcell),]#duplicate filter
    cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","PC1","PC2","treats","sex_alph","age","SES")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    cluster_counts_t <- cluster_counts_t[!rownames(cluster_counts_t) %in% genestoremove$genes,]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

    adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age","PC1","PC2","SES")])
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    save(adjusted_counts, file=opfn)
})

combatrun="SES_PCs_sex_age_adjusted"
mclapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
    cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata_bf <- data.frame(cluster_metadata_sce)
    cluster_metadata_var <- cluster_metadata_bf[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age","SES")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))
    adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("sex_alph","age","PC1","PC2","SES")])
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    save(adjusted_counts, file=opfn)
})

genestoremove <- ldply(lapply(names(counts_ls), function(cluster){
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
max=unique(rownames(which(adjusted_counts >= .Machine$integer.max, arr.ind = TRUE)))
if(is.null(max)){
    max=NA
}
return(data.frame(clus=cluster,genes=max))
}),data.frame)
genestoremove #35

combatrun="SES_PCs_sex_age_adjusted_generem"
mclapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
    cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    genestoremove_c <- subset(genestoremove, clus==cluster)
    cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","PC1","PC2","sex_alph","age","SES")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    cluster_counts_t <- cluster_counts_t[!rownames(cluster_counts_t) %in% genestoremove_c$genes,]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

    adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("sex_alph","age","PC1","PC2","SES")])
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    save(adjusted_counts, file=opfn)
})

#############################################################
############################################################
opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")

load(opfn)

combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
run="SES_PCs_sex_age_and_treats_generem"
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

#to test
outFolder=paste0(baseoutFolder,run,"/test/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

combatrun="SES_PCs_sex_age_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/test/")
run="SES_PCs_sex_age_generem"
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

runPFAS=TRUE
zingeRrun=FALSE
iselcov=FALSE
if(iselcov){
    run=paste0(run,"_iselcov")
}
SEScov=FALSE
if(SEScov){
    run=paste0(run,"_SEScov")
}
WHRcov=FALSE
if(WHRcov){
    run=paste0(run,"_WHRcov")
}
PSScov=FALSE
if(PSScov){
    run=paste0(run,"_PSScov")
}
if(zingeRrun){
    run=paste0(run,"_zingeR")
}

lapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    #cluster_metadataL <- left_join(cluster_metadata,unique(eigenvec2[,c("Sample_ID","Lead","WHR")]),by="Sample_ID")
    #rownames(cluster_metadataL) <- cluster_metadataL$rowid
    #cluster_metadata <- cluster_metadataL

    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)

    if(runPFAS){
        cluster_metadata_pfas <- cluster_metadata %>% 
        rowwise() %>% 
        mutate(
        sumPFAS = sum(!!!syms(combPFAS),na.rm=T)
        )
    cluster_metadata <- as.data.frame(cluster_metadata_pfas)
    cluster_metadata <- transform(cluster_metadata, log_sumPFAS=log2(sumPFAS+1))
    rownames(cluster_metadata) <- cluster_metadata$rowid
    #adjusted_countsP <- adjusted_counts[rownames(adjusted_counts) %in% genes$hgnc_symbol, ]
    #keep <- rowSums(adjusted_countsP,na.rm=T) >= (ncol(adjusted_countsP) / 2) #filter to keep genes expressed in 50% of samples
    #adjusted_counts <- unlist(adjusted_countsP[keep, ])
    }

    for(i in treatmentsfirst){
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        lapply(c("ISEL_Mean","PSS_all_mean","cytocomp"),function(var){
            #var="pr_comp"
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
                cat("running deseq ",cluster, var,i," \n")
        #if(var=="Lead"){
        #    adjusted_counts <- adjusted_counts[rownames(adjusted_counts) %in% genes$hgnc_symbol, ]
        #}

        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)
        if(var=="factor_HS_CRP"){
        cluster_metadata_t <- subset(cluster_metadata_t, HS_CRP<10) #advised to remove as likely an infection
        cluster_metadata_t <- transform(cluster_metadata_t, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
        cluster_metadata_t <- within(cluster_metadata_t, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
        }
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        if(SEScov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + SES + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","SES",var)]
        }
        if(iselcov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ISEL_Mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","ISEL_Mean",var)]
        }
        if(WHRcov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + WHR + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","WHR",var)]
        }
        if(PSScov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + PSS_all_mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","PSS_all_mean",var)]
        }
        if(var=="CVDRISK"){
        design <-  paste0("~ PC1 + PC2 + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2",var)]
        if(SEScov){
        design <-  paste0("~ PC1 + PC2 + SES + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","SES",var)]
        }
        if(WHRcov){
        design <-  paste0("~ PC1 + PC2 + WHR + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","WHR",var)]
        }
        if(iselcov){
        design <-  paste0("~ PC1 + PC2 + ISEL_Mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","ISEL_Mean",var)]
        }
        } #end of cvdrisk
        if(var=="age"){
        design <-  paste0("~ PC1 + PC2 + sex_alph + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph",var)]
        if(SEScov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + SES + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","SES",var)]
        }
        if(WHRcov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + WHR + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","WHR",var)]
        }
        if(iselcov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + ISEL_Mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","ISEL_Mean",var)]
        }
        if(PSScov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + PSS_all_mean + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","PSS_all_mean",var)]
        }
        } # end of age
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- try(DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        )
        if(inherits(dds, "try-error")){
          print("An error occurred")
        } else {
        if(zingeRrun){
                    designb = model.matrix( as.formula(design) , cluster_metadata_var)
        weights <- zeroWeightsLS(counts=cluster_counts_t, design=designb, maxit=200, normalization="DESeq2_poscounts", colData=cluster_metadata_var, designFormula=as.formula(design))
        assays(dds)[["weights"]]=weights
        dds = DESeq2::estimateSizeFactors(dds, type="poscounts")
        dds = estimateDispersions(dds)
        dds = nbinomWaldTest(dds, betaPrior=FALSE, useT=TRUE, df=rowSums(weights)-2)
        } else {
        dds <- DESeq(dds,parallel=TRUE)
        }

        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
        save(dds, file=opfn)

        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'baseMean', res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier','baseMean','padj','pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$var=var
        sub.table$cluster=cluster
        sub.table$treats =i
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$tested_genes <- paste(nrow(sub.table))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
        rm(table,sigDEGs,sub.table,res,dds,cluster_metadata_var)
        gc()
            }
            }
        })
    }
})


for (var in c(zcytokines)){
    #for (var in c("PSS_all_mean","ISEL_Mean","cytocomp")){
    for (i in c(treatmentsfirst)){
    if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
        cat("running ",var," ",i,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    }
    }
}

library(flextable)
library(ftExtra)
library(rlist)
library(officer)
#[!psychvarstorun %in% puberty]
 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
zoom="allvars"

for (i in c(treatmentsfirst)){
#subvars <- ldply(lapply(c(variables_df$variable), function(var){
subvars <- ldply(lapply(c(allvars), function(var){
    cat("running",i,var,"\n")
    if(isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"),select=c(1,8))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newls <- newls[sapply(newls, nrow)>0]
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,list.cbind(newls))
            newtable$description <- ifelse(is.na(newtable$description) | newtable$description=="", var,newtable$description)
            return(newtable)
        }
    }
}),data.frame)
subvars <- unique(subvars)

if( ncol(subvars) > 0){
dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
    #path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,"PFAS.png"))
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",zoom,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=my.max)<50],drop=F]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",zoom,".",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("numind|ngenes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() 
#dft <- align(dft, i = 1, j = NULL, align = "left", part = "body")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",zoom,".",run,".ndegonly.png"))
}
}

zoom="allvars"
for (i in c(treatmentsfirst)){
subvars <- ldply(lapply(c(allvars), function(var){
    cat("running",i,var,"\n")
    if(isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
        deseqres <- transform(deseqres, regulation=ifelse(logFC>0, "up", "down"))
        dessig <- subset(deseqres, padj<0.1)
        return(dessig)
    }
    }
}),data.frame)

# Count the number of genes in each category
data_summary <- ddply(subvars, c("cluster","var","regulation"),plyr::summarize,
    nDEG=length(regulation))
data_summary <- transform(data_summary, nDEG=ifelse(regulation=="up",nDEG,-(nDEG)))
group_colors <- c("up" = "orange","down" = "purple")

p <- ggplot(data_summary, aes(x = cluster, y = nDEG, fill = regulation)) +
  geom_col() + #position = "dodge"
  coord_flip() +
  geom_hline(yintercept=0)+
  #facet_wrap(.~var)+
  scale_fill_manual(values = group_colors) +
  labs(title = paste0("Up- and Down-Regulated Genes in ",i), x = "Cluster", y = "Number of DEGs", fill = "Regulation") +
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 8, height = 8, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",zoom,".",run,".DEGupvdownbar.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
}



for (i in c(treatmentsfirst)){
#subvars <- ldply(lapply(c(variables_df$variable), function(var){
subvars <- ldply(lapply(c(PFASvars), function(var){
    cat("running",i,var,"\n")
    if(isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"),select=c(1,7))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newls <- newls[sapply(newls, nrow)>0]
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,list.cbind(newls))
            newtable$description <- ifelse(is.na(newtable$description) | newtable$description=="", var,newtable$description)
            return(newtable)
    }
}),data.frame)
subvars <- unique(subvars)
if( ncol(subvars) > 0){
dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
    #path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,"PFAS.png"))
  path = paste0(outFolder,"figures/",i,".",zoom,".",run,".png"))
}
}

#############NOT UPDATED , maybe updated ?
library(Vennerable)

#want to plot isel vs no isel in DESeq model
run="SES_PCs_sex_age_and_treats_generem"
outFolder=paste0(baseoutFolder,run,"/")
sesrun=paste0(run,"_iselcov")
var="PSS_all_mean"
threshold=0.10

for (i in treatmentsfirst){
#lapply(c("PSS_all_mean"),function(var){
    allvardeseq <- ldply(lapply(names(counts_ls),function(cluster) {
        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
        load(opfn)
        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'pvalue', 'logFC','SE')
        deseqres <- transform(sub.table, padj=p.adjust(pvalue,method="BH"),var=var,cluster=cluster,treats=i)
        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",sesrun,".RData")
        load(opfn)
        res <- results(dds)
        sub.table <- data.frame(res@rownames, res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'pvalue', 'logFC','SE')
        ses_deseqres <- transform(sub.table, padj=p.adjust(pvalue,method="BH"),var=var,cluster=cluster,treats=i)
        deseqres <- transform(deseqres, model="noSocialSupportCov")
        ses_deseqres <- transform(ses_deseqres, model="SocialSupportCov")
        merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var","treats"))
        subvars <- transform(merged_treat, noisel_zscore=logFC.x/SE.x, isel_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3noisel_sig", ifelse(padj.y<=threshold, "2isel_sig", "1Not_Sig")))))
        return(subvars)
        }),data.frame)
group_colors <- c("1Not_Sig" = "grey","3noisel_sig" = "#E34234", "2isel_sig" = "turquoise4", "4both" = "#7851A9")
p <- ggplot(allvardeseq, aes(x=isel_zscore, y=noisel_zscore)) +
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".iselvsnoisel.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
#per cluster
p <- ggplot(allvardeseq, aes(x=isel_zscore, y=noisel_zscore)) +
  theme_bw()+
  facet_wrap(.~cluster,ncol=3)+
  geom_point(aes(color=sig))+ #aes(color=sig)
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".iselvsnoisel_bycluster.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

noSESressig <- subset(allvardeseq,padj.x<0.1)
noSESressig <- transform(noSESressig, model=model.x)
SESressig <- subset(allvardeseq,padj.y<0.1)
SESressig <- transform(SESressig, model=model.y)

if(dim(noSESressig)[1]>0 & dim(SESressig)[1]>0){
df <- unique(rbind(noSESressig[,c("identifier","model")],SESressig[,c("identifier","model")]))
d_list <- split(df$identifier,df$model)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0(outFolder,"figures/",project,".",var,"-",i,".noiselvsisel_venn.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
print(p)
dev.off()
}

#bar plot version
subvars_sig <- subset(allvardeseq, !sig=="1Not_Sig")
count_df <- subvars_sig %>%
  group_by(cluster, sig) %>%
  dplyr::summarise(count = n(), .groups = "drop")

group_colors <- c("3noisel_sig" = "#E34234", "2isel_sig" = "turquoise4", "4both" = "#7851A9")

fig0 <- ggplot(count_df, aes(x = count, y = cluster, fill = sig)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
    #facet_grid(.~term)+
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("DEGs Overlap Between SocialSupport and no SocialSupport covariate model", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
png(width = 12, height = 7, file=paste0(outFolder,"figures/noiselvsisel",".",var,"-",i,".stackedbar.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(fig0)
dev.off()
#}) #var loop
}

#SES vs no SES
run="SES_PCs_sex_age_and_treats_generem"
outFolder=paste0(baseoutFolder,run,"/")
sesrun=paste0(run,"_SEScov")
threshold=0.10

for (i in treatmentsfirst){
    allvardeseq <- ldply(lapply(c("PSS_all_mean","ISEL_Mean"),function(v){
    cat("running",i,"\n")
    deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",run,".txt"))
    ses_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",sesrun,".txt"))
    merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var","treats"))
    subvars <- transform(merged_treat, noses_zscore=logFC.x/SE.x, ses_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3noses_sig", ifelse(padj.y<=threshold, "2ses_sig", "1Not_Sig")))))
}), data.frame)
group_colors <- c("1Not_Sig" = "grey","3noses_sig" = "#E34234", "2ses_sig" = "turquoise4", "4both" = "#7851A9")

lapply(c("PSS_all_mean","ISEL_Mean"),function(v){
    subvardeseq <- subset(allvardeseq, var==v)
p <- ggplot(subvardeseq, aes(x=ses_zscore, y=noses_zscore)) +
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",v,".",run,".sesvsnoses.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
#per cluster
p <- ggplot(subvardeseq, aes(x=ses_zscore, y=noses_zscore)) +
  theme_bw()+
  facet_wrap(.~cluster,ncol=3)+
  geom_point(aes(color=sig))+ #aes(color=sig)
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",v,".",run,".sesvsnoses_bycluster.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
#venn
deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",run,".txt"))
ses_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",sesrun,".txt"))
deseqres <- transform(deseqres, model="noSESCov")
ses_deseqres <- transform(ses_deseqres, model="SESCov")
noSESressig <- subset(deseqres,padj<0.1)
SESressig <- subset(ses_deseqres,padj<0.1)
if(dim(noSESressig)[1]>0 & dim(SESressig)[1]>0){
df <- unique(rbind(noSESressig[,c("identifier","model")],SESressig[,c("identifier","model")]))
d_list <- split(df$identifier,df$model)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0(outFolder,"figures/",project,".",v,"-",i,".nosesvsses_venn.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
print(p)
dev.off()
}
#barplot
subvars_sig <- subset(subvardeseq, !sig=="1Not_Sig")
count_df <- subvars_sig %>%
  group_by(cluster, sig) %>%
  dplyr::summarise(count = n(), .groups = "drop")

group_colors <- c("3noses_sig" = "#E34234", "2ses_sig" = "turquoise4", "4both" = "#7851A9")

fig0 <- ggplot(count_df, aes(x = count, y = cluster, fill = sig)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
    #facet_grid(.~term)+
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("DEGs Overlap Between SES and no SES covariate model", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
png(width = 12, height = 7, file=paste0(outFolder,"figures/nosesvsses",".",v,"-",i,".stackedbar.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(fig0)
dev.off()
})
}

#want to plot WHR vs no WHR in DESeq model
run="SES_PCs_sex_age_and_treats_generem"
outFolder=paste0(baseoutFolder,run,"/")
sesrun=paste0(run,"_WHRcov")
var="PSS_all_mean"
threshold=0.10

for (i in treatmentsfirst){
    allvardeseq <- ldply(lapply(c("cytocomp","PSS_all_mean","ISEL_Mean"),function(var){
    cat("running",i,var,"\n")
    deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    ses_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",sesrun,".txt"))
    merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var","treats"))
    subvars <- transform(merged_treat, noWHR_zscore=logFC.x/SE.x, WHR_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3noWHR_sig", ifelse(padj.y<=threshold, "2WHR_sig", "1Not_Sig")))))
}), data.frame)
group_colors <- c("1Not_Sig" = "grey","3noWHR_sig" = "#E34234", "2WHR_sig" = "turquoise4", "4both" = "#7851A9")

lapply(c("cytocomp","PSS_all_mean","ISEL_Mean"),function(v){
        cat("plotting",i,v,"\n")
    subvardeseq <- subset(allvardeseq, var==v)
p <- ggplot(subvardeseq, aes(x=WHR_zscore, y=noWHR_zscore)) +
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",v,".",run,".WHRvsnoWHR.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
#per cluster
p <- ggplot(subvardeseq, aes(x=WHR_zscore, y=noWHR_zscore)) +
  theme_bw()+
  facet_wrap(.~cluster,ncol=3)+
  geom_point(aes(color=sig))+ #aes(color=sig)
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",v,".",run,".WHRvsnoWHR_bycluster.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
})
}

for (i in treatmentsfirst){
    lapply(c(zcytokines,firstrunvars),function(var){
    cat("running",i,"\n")
    deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    ses_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",sesrun,".txt"))
    merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var","treats"))
    subvars <- transform(merged_treat, noWHR_zscore=logFC.x/SE.x, WHR_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3noWHR_sig", ifelse(padj.y<=threshold, "2WHR_sig", "1Not_Sig")))))

deseqres <- transform(deseqres, model="noWHRCov")
ses_deseqres <- transform(ses_deseqres, model="WHRCov")

noSESressig <- subset(deseqres,padj<0.1)
SESressig <- subset(ses_deseqres,padj<0.1)

if(dim(noSESressig)[1]>0 & dim(SESressig)[1]>0){
df <- unique(rbind(noSESressig[,c("identifier","model")],SESressig[,c("identifier","model")]))

d_list <- split(df$identifier,df$model)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0(outFolder,"figures/",project,".",var,"-",i,".noWHRvsWHR_venn.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 600)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
print(p)
dev.off()
}

#bar plot version
group_colors <- c("3noWHR_sig" = "#E34234", "2WHR_sig" = "turquoise4", "4both" = "#7851A9")
subvars_sig <- subset(subvars, !sig=="1Not_Sig")

count_df <- subvars_sig %>%
  group_by(cluster, sig) %>%
  dplyr::summarise(count = n(), .groups = "drop")

fig0 <- ggplot(count_df, aes(x = count, y = cluster, fill = sig)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
    #facet_grid(.~term)+
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("DEGs Overlap Between WHR and no WHR covariate model", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
png(width = 12, height = 7, file=paste0(outFolder,"figures/WHRvsnoWHR",".",var,"-",i,".stackedbar.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(fig0)
dev.off()

})
}

#for (i in c(treatmentsfirst)){
#    lapply(names(counts_ls), function(cluster){
#        lapply(c(PFASvars), function(var){
        cluster="C4"
        var="PFOA"
        i="RNA-CTRL"
        cat("running",i,var,cluster,"\n")
        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
        load(file=opfn)
        res <- results(dds)
        png(width = 6, height = 6, file=paste0(outFolder,"figures/",var,".",cluster,".",i,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
        plotMA(res)
        dev.off()
#        })
#    })
#}

#also MA plots -- I didn't output basemean values ...

        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))


png(width = 6, height = 6, file=paste0(outFolder,"figures/",var,"-",i,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggplot(dat, aes(baseMean, log2FoldChange, col=sig, label=SYMBOL)) + 
        geom_point() + scale_x_log10() + 
        xlab("mean of normalized counts") + 
        ylab("log fold change") + 
        geom_hline(yintercept=0, col="grey40") + 
        scale_color_manual(values=c("grey60", "blue")) +
        geom_point(data=tab, shape=1, size=5, show.legend=FALSE) + 
        geom_label_repel(data=tab, nudge_x = 1, nudge_y = 2*sign(tab$log2FoldChange), show.legend=FALSE)
dev.off()

#volcano
library(EnhancedVolcano)

        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
#for (c in names(counts_ls)){
    c="C4"
    cat("running ", c, "\n")
    sub.table <- subset(deseqres, cluster==c)
    sub.table <- sub.table[order(sub.table$padj,-abs(sub.table$logFC)), ]
    xmax=max(sub.table$logFC)+0.5
    xmin=min(sub.table$logFC)-0.5
    topp1 <- min(subset(sub.table, padj > quantile(padj, prob = 1 -99/100,na.rm=T))$pvalue)
    tolab <- c(unique(head(sub.table,n=20)$identifier))
    png(width = 8, height = 8, file=paste0(outFolder,"figures/dge_volcano-",var,"-",i,".",run,"-",c,".png"), pointsize=12, 
          bg = "transparent", units = "in", res = 1200)
    p <- EnhancedVolcano(sub.table,
      lab = sub.table$identifier,
      x = 'logFC',
      y = 'pvalue',
      title=paste0(c),
      subtitle = NULL,
      xlim = c(xmin, xmax),
      pCutoffCol= 'padj',
      pCutoff = 0.1,
      FCcutoff = 0.5,
      labSize = 3.0,
      #pointSize = c(ifelse(sub.table$gene_symbol %in% tolab, 5,1.5)),
      col=c('black', 'black', 'blue', 'red3'),
      hline = c(topp1),
      hlineCol = c('green'),
      selectLab = tolab,
      legendPosition = 'none',
      drawConnectors = TRUE)
    print(p)
    dev.off()
#}

#plotting LPS vs CONTROL
run="SES_PCs_sex_age_and_treats_generem"
outFolder=paste0(baseoutFolder,run,"/")
threshold=0.10
treat="RNA-LPS"
control="RNA-CTRL"
subvars <- ldply(lapply(c("PSS_all_mean","ISEL_Mean"), function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",control,".",run,".txt"))
        treat_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",treat,".",run,".txt"))
        merged_treat <- rbind(ctrl_deseqres,treat_deseqres)
        return(merged_treat)
}),data.frame)
#intrun="SES_PCs_sex_age_and_treats_generem_treatint"
##intrun="SES_PCs_sex_age_and_treats_generem_zingeR_treatint"
#intoutFolder=paste0(baseoutFolder,intrun,"/")
#myDir <- paste0(intoutFolder,"deseqres/")
#filenames <- list.files(myDir) #file list from directory
#filenames <- filenames[grep(".treatinteraction.txt",filenames)]
#data_names <- gsub(paste0(".",run,".treatinteraction.txt"), "", filenames) #remove file ending
#data_names <- gsub(paste0(project,".",resset,".",dimset,".deseqres_"), "", data_names)
#for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
#ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
#names(ddf) <- c(data_names)
#deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#names(deseqres)[8] <- c("treats")
#deseqres <- transform(deseqres, treats="interaction")
#deseqresL <- subset(deseqres, var %in% variablesL)

#all <- rbind(subvars,deseqresL)
library(gridExtra)
all <- subvars
all_sig <- subset(all, padj<0.1)
lapply(split(all_sig,all_sig$var),function(i){
    myvar=unique(i$var)
    cat("running ",myvar)
    df <- unique(i[,c("identifier","treats")])
    d_list <- split(df$identifier,df$treats)
    V_d_list <- Venn(d_list)
    Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
    #Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
    Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
    gp <- VennThemes(Vennlist)
    png(width =6, height = 6, file=paste0(outFolder,"figures/",resset,".",dimset,".",myvar,".ctrlvslps_venn.png"), pointsize=12, 
    #png(width =6, height = 6, file=paste0(outFolder,"figures/",resset,".",dimset,".",myvar,".ctrlvsplsvsint_venn.png"), pointsize=12, 
          bg = "transparent", units = "in", res = 1200)
    par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
    p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE)); grid.text(paste0(myvar), y=0.92, gp=gpar(col="black", cex=2))
    print(p)
    dev.off()
})

#plotting LPS vs CONTROL scatter
subvars <- ldply(lapply(c("PSS_all_mean","ISEL_Mean"), function(v){
    cat("running",treat,"vs",control,v,"\n")
    ctrl_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",control,".",run,".txt"))
    treat_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",treat,".",run,".txt"))
    merged_treat <- merge(ctrl_deseqres,treat_deseqres,by=c("identifier","cluster","var"),all=T)
    merged_treat_all <- transform(merged_treat, CTRL_zscore=logFC.x/SE.x, LPS_zscore=logFC.y/SE.y)
    merged_treat_all <- merged_treat_all %>% mutate(sig = case_when(
    padj.x<threshold & padj.y<threshold ~ "4CTRLandLPS_sig",
    padj.x < threshold ~ "3CTRL_sig",
    padj.y < threshold ~ "2LPS_sig",    
    ))
    merged_treat_all$sig[is.na(merged_treat_all$sig)] <- "1Not_Sig"
    return(merged_treat_all)
}),data.frame)

cols <- c("4CTRLandLPS_sig" = "green","1Not_Sig" = "grey","3CTRL_sig" = "red", "2LPS_sig" = "blue")
p <- ggplot(subvars, aes(x=CTRL_zscore, y=LPS_zscore)) +
  facet_wrap(.~var)+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl_z.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
#ploting by var and cluster
p <- ggplot(subvars, aes(x=CTRL_zscore, y=LPS_zscore)) +
  facet_grid(cluster~var)+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 14, height = 14, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl_z_mat.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

#for logfc
subvars_sig1 <- subset(subvars, padj.x<threshold | padj.y<threshold)
p <- ggplot(subvars_sig1, aes(x=logFC.x, y=logFC.y)) +
  facet_wrap(.~var, scales="free")+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  xlab("CTRL logFC") + ylab("LPS logFC")+
  geom_abline()+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
#ploting by var and cluster
p <- ggplot(subvars, aes(x=logFC.x, y=logFC.y)) +
  facet_grid(cluster~var, scales="free")+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  xlab("CTRL logFC") + ylab("LPS logFC")+
  geom_abline()+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 14, height = 14, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl_mat_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
