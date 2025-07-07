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
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",11,0.2) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
# set new output dir for filtered out unmatched figures
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
figfn <- paste(base,method,"_pseudobulk_ctrl/adjusted/figures/", myvar,"_histogram_distribution",".png",sep="")
png(figfn, width=1000, height=1000, res=120)
#print(plot_grid(fig1, fig3, ncol=2, labels="AUTO", label_fontface="plain"))
print(p1)
dev.off()

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

outdir=paste0(base,"5b_IdenCelltype_",method,"/")
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

#subset for 20cell count filt
counts <- plyr::count(sc@meta.data,c("Library","letter_clusters","Sample_ID"))
names(counts)[4] <- "btic_cellcounts"
sc@meta.data <- left_join(sc@meta.data,counts,by=c("Library","letter_clusters","Sample_ID"))
rownames(sc@meta.data) <- sc@meta.data$NEW_BARCODE
before=dim(sc@meta.data)[1]
sc <- subset(sc, subset=btic_cellcounts>=20)
after=dim(sc@meta.data)[1]
cat("removed ", before-after, "combos","\n") #3971
table(sc@meta.data$letter_clusters,sc@meta.data$treats)
    RNA-CTRL RNA-LPS RNA-LPS-DEX
  C0    77284   83358       19300
  C1    57994   58139       13412
  C2    44634   49984       14593
  C3    34712   40317       10741
  C4    37170   28522       11316
  C5    20908   23289        5282
  C6      603     830         204

sc[["sex_male"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.male$symbol])
sc[["sex_female"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.female$symbol])

opfn <- paste0(outFolder,project,".",resset,".",dimset,".wavefilt.bticfilt.seurat.RDS")
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

cellcount <- as.data.frame(table(sc@meta.data$letter_clusters, sc@meta.data$orig.ident))
lowcell <- subset(cellcount, Freq<3000) #normally 3000 which would remove C6 too
if(dim(lowcell)[1]==0){cat( "no low cell counts")}else{sc <-subset(x = sc, subset = letter_clusters %in% unique(lowcell$Var1), invert = TRUE)}

sce <- as.SingleCellExperiment(sc)
#opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".SingleCellExperiment.RDS")
#write_rds(sce, opfn)
rm(sc)
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
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt_proteincoding.RData")
save(counts_ls,metadata_ls, file=opfn)

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
    data <- data[!rownames(data) %in% geneIDs.sex$symbol, ]
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
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
save(counts_ls,metadata_ls, file=opfn)
load(opfn)

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

#only using SES var, also accounting for PC1 and PC2
combatrun="SES_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
if (!file.exists(baseoutFolder)) dir.create(baseoutFolder, showWarnings=F)

mclapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
    cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","PC1","PC2","treats","sex_alph","age","SES")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

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
genestoremove #237 genes

combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
mclapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
    cat("running ", cluster, "\n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata_var <- cluster_metadata[,c("Sample_ID","BATCH","PC1","PC2","treats","sex_alph","age","SES")]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_var))]
    cluster_counts_t <- cluster_counts_t[!rownames(cluster_counts_t) %in% genestoremove$genes,]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

    adjusted_counts <- ComBat_seq(cluster_counts_t, batch=cluster_metadata_var$BATCH, group=NULL, covar_mod=cluster_metadata_var[,c("treats","sex_alph","age","PC1","PC2","SES")])
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    save(adjusted_counts, file=opfn)
})

combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
allPFAS <- psychvarstorun[c(45:65)]
combPFAS <- allPFAS[seq(1,21,3)]
PFASvars <- c(allPFAS,"sumPFAS","log_sumPFAS")
run="SES_PCs_sex_age_and_treats_generem"
outFolder=paste0(baseoutFolder,run,"/")
runPFAS=FALSE

run=paste0(run,"_proteincoding")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

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
    cluster_metadataL <- left_join(cluster_metadata,unique(eigenvec2[,c("Sample_ID","Lead")]),by="Sample_ID")
    rownames(cluster_metadataL) <- cluster_metadataL$rowid
    cluster_metadata <- cluster_metadataL

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
    adjusted_countsP <- adjusted_counts[rownames(adjusted_counts) %in% genes$hgnc_symbol, ]
    keep <- rowSums(adjusted_countsP,na.rm=T) >= (ncol(adjusted_countsP) / 2) #filter to keep genes expressed in 25% of samples
    adjusted_counts <- unlist(adjusted_countsP[keep, ])
    }

    #d <- data.frame(number_ind=length(unique(cluster_metadata[,c("Sample_ID")])))
    #fwrite(d, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid_all/",project,".",cluster,".sampleid_all.",run,".txt"))

    for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        lapply("Lead",function(var){
            #c(psychvarstorun,"factor_HS_CRP",updated_vars,FCDEM_11_8) variables_df$variable
            #var="pr_comp"
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
                cat("running deseq ",var,i," \n")
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
        if(var=="CVDRISK"){
        design <-  paste0("~ PC1 + PC2 + ",var)
        if(SEScov){
        design <-  paste0("~ PC1 + PC2 + SES + ",var)
        }
        }
        if(iselcov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + isel + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","isel",var)]
        }
        if(var=="age"){
        design <-  paste0("~ PC1 + PC2 + sex_alph + ",var)
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph",var)]
        if(SEScov){
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","SES",var)]
        }
        }
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
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

        #fwrite(data.frame(Sample_ID=cluster_metadata_var[,c("Sample_ID")]), sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid/",project,".",resset,".",dimset,".",cluster,".sampleid_",var,"-",i,".",run,".txt"))

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
})


for (var in c(variables_df$variable)){
    #for (var in "Lead"){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
    if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
        cat("running ",var," ",i,"\n")
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
}

library(flextable)
library(ftExtra)
library(rlist)
library(officer)
#[!psychvarstorun %in% puberty]
 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
subvars <- ldply(lapply(c(variables_df$variable), function(var){
#subvars <- ldply(lapply(PFAS, function(var){
    cat("running",i,var,"\n")
    if(isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
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
            newls <- newls[sapply(newls, nrow)>0]
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,topvar=ifelse(var %in% variablesL,"Y",""),list.cbind(newls))
            return(newtable)
        }
    }
}),data.frame)

if( ncol(subvars) > 0){
dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
    #path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,"PFAS.png"))
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=my.max)<50],drop=F]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("numind|ngenes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() 
#dft <- align(dft, i = 1, j = NULL, align = "left", part = "body")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".ndegonly.png"))
}
}

#want to plot isel vs no isel in DESeq model
sesrun="SES_PCs_isel_sex_age_and_treats_generem"
run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
sesoutFolder=paste0(baseoutFolder,sesrun,"/")
threshold=0.10

for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
    cat("running",i,"\n")
    deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    ses_deseqres <- fread(paste0(sesoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",sesrun,".txt"))
    merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var"))
    subvars <- transform(merged_treat, noisel_zscore=logFC.x/SE.x, isel_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3noisel_sig", ifelse(padj.y<=threshold, "2isel_sig", "1Not_Sig")))))
#didnt update this one for the colored by sig
p <- ggplot(subvars, aes(x=isel_zscore, y=noisel_zscore)) +
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".iselvsnoisel.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
}
#want to plot SES vs no SES in DESeq model
sesrun="SES_PCs_SES_sex_age_and_treats_generem"
run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
sesoutFolder=paste0(baseoutFolder,sesrun,"/")
threshold=0.10
variablesL <- variablesL[-c(1,3)]

for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
subvars <- ldply(lapply(variablesL, function(var){
    cat("running",i,var,"\n")
    if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
        ses_deseqres <- fread(paste0(sesoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",sesrun,".txt"))
        merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var"))
        merged_treat <- transform(merged_treat, noSES_zscore=logFC.x/SE.x, SES_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3noSES_sig", ifelse(padj.y<=threshold, "2SES_sig", "1Not_Sig")))))
        return(merged_treat)
    }
}),data.frame)
#didnt update this one for the colored by sig
p <- ggplot(subvars, aes(x=SES_zscore, y=noSES_zscore)) +
  facet_wrap(.~var)+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".sesvsnoses.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
}
#plotting DSES09 sep as it was less correlated - due to specific cluster?
subvars <- ldply(lapply(c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX"), function(i){
var="DSES_09"
deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
ses_deseqres <- fread(paste0(sesoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",sesrun,".txt"))
merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var","treats"))
merged_treat <- transform(merged_treat, noSES_zscore=logFC.x/SE.x, SES_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3noSES_sig", ifelse(padj.y<=threshold, "2SES_sig", "1Not_Sig")))))
return(merged_treat)
}),data.frame)
p <- ggplot(subvars, aes(x=noSES_zscore, y=SES_zscore)) +
  facet_wrap(treats~cluster, scales="free")+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",var,".",run,".sesvsnoses.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

#identifying SES only sig DEGs
sesonlydegs <- ldply(lapply(c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX"), function(i){
    subvars <- ldply(lapply(variablesL, function(var){
        cat("running",i,var,"\n")
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
            deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
            ses_deseqres <- fread(paste0(sesoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",sesrun,".txt"))
            merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var","treats"))
            merged_treat <- transform(merged_treat, noSES_zscore=logFC.x/SE.x, SES_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3noSES_sig", ifelse(padj.y<=threshold, "2SES_sig", "1Not_Sig")))))
            return(merged_treat)
        }
    }),data.frame)
    sesonly_subvars <- subset(subvars, sig=="2SES_sig")
    return(sesonly_subvars)
}),data.frame)
length(unique(sesonlydegs$identifier))
table(sesonlydegs$var,sesonlydegs$treats)
sesonlydegs_ctrl <- subset(sesonlydegs,treats=="RNA-CTRL")
table(sesonlydegs$var,sesonlydegs$cluster)
head(unique(sesonlydegs_ctrl[order(sesonlydegs_ctrl$padj.y),]$identifier),n=20)

#fishers test for health vs psych vars in SES test
myDir <- paste0(sesoutFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0(project,".",resset,".",dimset,".deseqres_"),filenames)]
data_names <- gsub(".SES_PCs_SES_sex_age_and_treats_generem.txt", "", filenames) #remove file ending
data_names <- gsub("ALL.0.2.11.deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')[,analysis:=data_names[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
df <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) 

ses_deseqresd <- merge(df,cat_cov_m,by.x="var",by.y="variable")
ses_deseqresd <- transform(ses_deseqresd,identifier=gsub("-","_",identifier))
ses_deseqresdm <- melt(ses_deseqresd[,c("var","identifier","cluster","treats","health_vs_psych","padj")])
test <- reshape2::dcast(identifier ~ health_vs_psych,ses_deseqresdm, value.var="value")


ses_deseqresd_h <- subset(ses_deseqresd,health_vs_psych=="health")
ses_deseqresd_p <- subset(ses_deseqresd,health_vs_psych=="psych")
ses_deseqresd_sig <- subset(ses_deseqresd, padj<0.1)

t <- table(health=ses_deseqresd_sig$health_vs_psych=="health",psych=ses_deseqresd_sig$health_vs_psych=="psych", useNA = "always")[-3,-3]

t <- table(health=ses_deseqresd[ses_deseqresd$padj<0.1 & ses_deseqresd$health_vs_psych=="health",],psych=ses_deseqresd[ses_deseqresd$padj<0.1 & ses_deseqresd$health_vs_psych=="psych",], useNA = "always")[-3,-3]
df_fisher <- tryCatch(fisher.test(t), error=function(e) data.frame(estimate=NA,p.value=NA,conf.int=c(NA,NA)))
t1 <- unique(data.frame(contrast=paste0(i,"_",c),truetrue=tryCatch(t[2,2],error=function(x)NA),oddsratio=df_fisher$estimate,pval=df_fisher$p.value,
    CI_low=df_fisher$conf.int[1], CI_high=df_fisher$conf.int[2]))


# running treatment DEGs
withCOMBAT=FALSE
run="treatment_withCOMBAT"
run="treatment_noCOMBAT"
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid/"))) dir.create(paste0(outFolder,"sampleid/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid_all/"))) dir.create(paste0(outFolder,"sampleid_all/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

contrastdf <- data.frame(control=c("RNA-CTRL","RNA-LPS"),treatment=c("RNA-LPS","RNA-LPS-DEX"))
contrastdf <- transform(contrastdf, contrast=paste0(treatment,"_vs_",control))

lapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
#    if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt")) > 0)){
    cat("running deseq ",cluster," \n")
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    if(withCOMBAT){
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    } else {
        adjusted_counts <- cluster_counts
    }

    d <- data.frame(number_ind=length(unique(cluster_metadata[,c("Sample_ID")])))
    fwrite(d, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid_all/",project,".",cluster,".sampleid_all.",run,".txt"))

    mclapply(1:nrow(contrastdf),function(x) {
        con=contrastdf[x,]
        contrast=con$contrast
        cat("running",contrast,"\n")
        cluster_metadata_var <- cluster_metadata[,c("Sample_ID","treats")]
        cluster_metadata_var <- subset(cluster_metadata_var, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        indcount <- plyr::count(cluster_metadata_var,"Sample_ID")
        cluster_metadata_var <- subset(cluster_metadata_var, Sample_ID %in% subset(indcount, freq>1)$Sample_ID)
        design <-  paste0("~ Sample_ID + treats")
        cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

        dds <- DESeqDataSetFromMatrix(cluster_counts_t, 
                                      colData = cluster_metadata_var, 
                                      design = as.formula(design))
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",cluster,".",contrast,".",run,".RDS")
        saveRDS(dds, file=opfn)

        res <- results(dds, contrast=c("treats",con$treatment,con$control))
        sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
        names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
        sub.table <- sub.table[!is.na(sub.table$padj), ]
        cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        sub.table$cluster=cluster
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres.",contrast,".",run,".txt"))
        sigDEGs <- subset(sub.table,padj<fdr)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sigDEGs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs.",contrast,".",run,".txt"))
        table <- data.frame(cluster=cluster, contrast=contrast)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$tested_genes <- paste(nrow(sub.table))
        table$DEGs_FDR <- paste(nrow(sigDEGs))
        table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt"))
        })    
    #}
})

for (x in c(1:nrow(contrastdf))){
    con=contrastdf[x,]
    contrast=con$contrast
    cat("running ",contrast,"\n")
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres.",contrast,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres.",contrast,".",run,".txt"),header=T)
        }
        }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs.",contrast,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs.",contrast,".",run,".txt"),header=T)
        }    
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".sigDEGs.",contrast,".",run,".txt"))
    sub.table <- ldply(lapply(names(counts_ls), function(cluster){
        if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt"))){
        sub.table <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types.",contrast,".",run,".txt"),header=T)
        }
    }),data.frame)
    fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types.",contrast,".",run,".txt"))
}

subvars <- ldply(lapply(1:nrow(contrastdf), function(x){
    con=contrastdf[x,]
    contrast=con$contrast
    cat("running",contrast,"\n")
    if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types.",contrast,".",run,".txt"))){
    stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types.",contrast,".",run,".txt"))
    if(dim(subset(stats, DEGs_FDR_10>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
        df <- ddply(deseqres, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$DEGs_FDR_10)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newtable <- data.frame(contrast=contrast,list.cbind(newls))
            return(newtable)
        }
    }
}),data.frame)

dft <- subvars %>% flextable() %>% span_header()  

dft <- align(dft, i = 1:2, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".png"))

degcols <- grep("sigDEGs",colnames(subvars))
subsubvars <- subvars[,!grepl(paste(c(sapply(strsplit(colnames(subvars[,degcols[ apply(subvars[,degcols],MARGIN=2,FUN=max)<50]]),"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(subvars))]
dft <- subsubvars %>% flextable() %>% span_header() %>% set_caption(caption = run) 
cols <- seq(1,length(subsubvars),by=3)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1:2, j = c(2:length(subsubvars)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(subsubvars)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(cols), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("numind|ngenes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() %>% split_header()
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = c(1), border = border, part = "all")
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(subsubvars_degonly)), align = "center", part = "body")

#dft <- align(dft, i = 1, j = NULL, align = "left", part = "body")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".ndegonly.png"))

baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
run="treatment_withCOMBAT"
deseqres <- fread(paste0(baseoutFolder,run,"/deseqres/ALL.0.2.11.deseqres.RNA-LPS_vs_RNA-CTRL.treatment_withCOMBAT.txt"))
figuredir <- paste0(baseoutFolder,run,"/figures/")
for (c in names(counts_ls)){
    cat("running ", c, "\n")
    sub.table <- subset(deseqres, cluster==c)
    sub.table <- sub.table[order(sub.table$padj,-abs(sub.table$logFC)), ]
    xmax=max(sub.table$logFC)+0.5
    xmin=min(sub.table$logFC)-0.5
    topp1 <- min(subset(sub.table, padj > quantile(padj, prob = 1 -99/100,na.rm=T))$pvalue)
    tolab <- c(unique(head(sub.table,n=20)$identifier))
    png(width = 8, height = 8, file=paste0(figuredir,"dge_volcano-",run,"-",c,".png"), pointsize=12, 
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
}





#New table: var|#CTRL|#LPS|#interactiontest for 
run="SES_PCs_sex_age_and_treats_generem"
run="SES_PCs_sex_age_and_treats_generem_zingeR"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
variablesL <- c("chronic_sum","cytocomp", "SES","PSS_all_mean", "Chol_HDL", "BPd_avg", "isel", "pr_comp")  #for regular runs
secondrunvars=c(psychvarstorun[c(1:13,24:45)],"age")
treat="RNA-LPS"
control="RNA-CTRL"
subvars <- ldply(lapply(secondrunvars, function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",control,".",run,".txt"))
        treat_stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",treat,".",run,".txt"))
        merged_treat <- merge(ctrl_stats,treat_stats,by=c("symb","variable","cluster"))
        return(merged_treat)
}),data.frame)

intrun="SES_PCs_sex_age_and_treats_generem_treatint"
intrun="SES_PCs_sex_age_and_treats_generem_zingeR_treatint"
intoutFolder=paste0(baseoutFolder,intrun,"/")
myDir <- paste0(intoutFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".treatinteraction.txt",filenames)]
data_names <- gsub(paste0(".",run,".treatinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".stats_all_cell_types-"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
colnames(stats) <- c("symb","description","cluster","contrast","number_samples","number_individuals","gene_number","tested_genes","DEGs_FDR","sigDEGs")
stats_int <- transform(stats, symb_contrast=paste(symb,contrast,sep=":"))

merged_treat <- merge(subvars,stats_int,by=c("symb","cluster"))
df <- merged_treat[,c("symb","description","cluster","DEGs_FDR_10.x","DEGs_FDR_10.y","sigDEGs")]
colnames(df) <- c("symb","description","cluster","CTRL","LPS","interaction")
df50 <- subset(df, CTRL>49 | LPS>49 | interaction>49)
df50 <- df50[order(df50$symb,df50$cluster),]

dft <- df50 %>% flextable() %>% span_header() %>% set_caption(caption = run) 
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1, j = c(1:length(df50)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(df50)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(2,3,4,5), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".ctrlvslpsvsintdegs_table.png"))

#New table: var|#ind|#tested genes
treat="RNA-CTRL"
subvars <- ldply(lapply(secondrunvars, function(var){
    cat("running",treat,"vs",control,var,"\n")
    ctrl_stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",control,".",run,".txt"))
    df <- merge(ctrl_stats,variables_df,by=c("variable"))
    dft <- df[,c("symb","description","cluster","number_individuals","tested_genes")]
}),data.frame)
df50_ctrl <- subset(df, CTRL>49)
subsubvars <- subset(subvars, paste0(symb,"_",cluster) %in% paste0(df50_ctrl$symb,"_",df50_ctrl$cluster))
subsubvars <- subsubvars[order(subsubvars$symb,subsubvars$cluster),]

dft <- subsubvars %>% flextable() %>% span_header() %>% set_caption(caption = run) 
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1, j = c(1:length(subsubvars)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(subsubvars)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(2,3,4,5), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".",treat,".numindandgenes_table.png"))

treat="RNA-LPS"
df50_ctrl <- subset(df, LPS>49)
subsubvars <- subset(subvars, paste0(symb,"_",cluster) %in% paste0(df50_ctrl$symb,"_",df50_ctrl$cluster))
subsubvars <- subsubvars[order(subsubvars$symb,subsubvars$cluster),]

dft <- subsubvars %>% flextable() %>% span_header() %>% set_caption(caption = run) 
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1, j = c(1:length(subsubvars)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(subsubvars)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(2,3,4,5), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".",treat,".numindandgenes_table.png"))

treat="treatint"
df50_ctrl <- subset(stats_int, sigDEGs>49)
subsubvars <- subset(subvars, paste0(symb,"_",cluster) %in% paste0(df50_ctrl$symb,"_",df50_ctrl$cluster))
subsubvars <- subsubvars[order(subsubvars$symb,subsubvars$cluster),]

dft <- subsubvars %>% flextable() %>% span_header() %>% set_caption(caption = run) 
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1, j = c(1:length(subsubvars)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(subsubvars)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(2,3,4,5), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".",treat,".numindandgenes_table.png"))

#want to plot zingeR vs no regular in DESeq model
sesrun="SES_PCs_sex_age_and_treats_generem_zingeR"
run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
sesoutFolder=paste0(baseoutFolder,sesrun,"/")
threshold=0.10

for (i in c("RNA-CTRL","RNA-LPS")){
    varsrun<-na.omit(ldply(lapply(variables_df$variable, function(var){
    if(isTRUE(file.size(paste0(sesoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",sesrun,".txt")) > 0)){
    df <- data.frame(variable=var)
    } else {
    df <- data.frame(variable=NA)
    }
    return(df)
        }),data.frame))
    subvars <- ldply(lapply(varsrun$variable, function(var){
            cat("running",i,var,"\n")
    deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    ses_deseqres <- fread(paste0(sesoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",sesrun,".txt"))
    merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var"))
    merged_treat <- transform(merged_treat, regular_zscore=logFC.x/SE.x, zingeR_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3regular_sig", ifelse(padj.y<=threshold, "2zingeR_sig", "1Not_Sig")))))
    return(merged_treat)
    }),data.frame)

#didnt update this one for the colored by sig
p <- ggplot(subvars, aes(x=zingeR_zscore, y=regular_zscore)) +
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".zingeRvsregular.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
}
#plot #DEGs scaip1 vs 2 color and shape for cluster and variable
distinct_colors <- c("#2be360",
"#7508c6",
"#eac314",
"#33046e",
"#84a300",
"#c978ff",
"#f08e00",
"#0265bd",
"#c70008",
"#4ab8ff",
"#7d0e00",
"#60d9d0",
"#ff3a89",
"#007b47",
"#680043",
"#1f4500",
"#f4aff1",
"#795900",
"#ff945d",
"#3d1e0a")
distinct_colors <-c("#db6e00",
"#4329cb",
"#5fab00",
"#c331df",
"#019b32",
"#b95aff",
"#c2a700",
"#0036b7",
"#f9bc28",
"#6d0096",
"#79db87",
"#ff7dfa",
"#bcd05c",
"#014ab7",
"#9a9200",
"#688bff",
"#ff3123",
"#01b281",
"#a80044",
"#01aea8",
"#c14600",
"#57aeff",
"#6d0004",
"#0086ce",
"#ff9b6c",
"#0155a6",
"#6b6600",
"#1f2865",
"#d7b282",
"#81004c",
"#006a34",
"#d5bbf9",
"#273a00",
"#ffabc9",
"#724500",
"#714771",
"#ff928d",
"#55133b")
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,n), col=sample(col_vector, n))
    n <- 20
    distinct=sample(col_vector, n)

for (i in c("RNA-CTRL","RNA-LPS")){
    varsrun<-na.omit(ldply(lapply(variables_df$variable, function(var){
    if(isTRUE(file.size(paste0(sesoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",sesrun,".txt")) > 0)){
    df <- data.frame(variable=var)
    } else {
    df <- data.frame(variable=NA)
    }
    return(df)
        }),data.frame))
    subvars <- ldply(lapply(varsrun$variable, function(var){
    cat("running",i,var,"\n")
    if(file.exists(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))){
        stats <- fread(paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",run,".txt"))
        normal_stats <- fread(paste0(sesoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",sesrun,".txt"))
        merged_treat <- merge(normal_stats,stats,by=c("cluster","variable","symb"))
        names(merged_treat)[c(9,15)] <- c("zingeR_DEGs","Regular_DEGs")
        return(merged_treat)
    }
}),data.frame)

p <- ggplot(subvars, aes(x=zingeR_DEGs, y=Regular_DEGs)) +
  theme_bw()+
  geom_point(size=7,aes(shape=cluster,color=variable,fill=variable))+ #aes(color=sig)
  scale_shape_manual(values = c(7,8,15:18,25))+
  scale_color_manual(values = distinct_colors) +
  #scale_fill_discrete(guide="none")+
  geom_vline(xintercept = 50,linetype="dotted",colour="red")+
  geom_hline(yintercept = 50,linetype="dotted",colour="red")+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 10, height = 10, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".zingeRvsregular_numDEGs.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
}

#PFAS plotting only
run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
PFAS <- c("Log_PFOA_1", "Log_PFHxS_1", "Log_PFNA_1", "Log_LPFHpS_1","Log_PFDA_1", "Log_PFOS_1", "Log_PFUdA_1") 
allPFAS <- psychvarstorun[c(46:66)]
combPFAS <- allPFAS[seq(1,21,3)]
zoom="allPFAS"
#zoom="zcytokines"
#zoom="old_cytokines"
#zoom="Lead"
#sumPFAS <- allPFAS[seq(1,21,3)]
for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
#subvars <- ldply(lapply(old_cytokines, function(var){
subvars <- ldply(lapply("Lead", function(var){
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
            newtable <- data.frame(variable=var, list.cbind(newls))
            return(newtable)
    }
}),data.frame)
  fwrite(subvars, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",zoom,"-",i,".",run,".txt"))


if( ncol(subvars) > 0){
dft <- subvars %>% flextable() %>% span_header()
dft <- align(dft, i = 1, j = NULL, align = "center", part = "header")
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",run,".",zoom,".png"))
}
}

##
library(ggrastr)
plotDF <- ldply(lapply(c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX"),function(i){
    resall <- ldply(lapply(allPFAS,function(var){
    res <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    ntest <- nrow(res)
       res <- res%>%
           arrange(pvalue)%>%
           mutate(observed=-log10(pvalue), expected=-log10(ppoints(ntest)))
           return(res)
        }),data.frame)
    return(resall)
    }), data.frame)
col2 <- hue_pal()(9) ### default ggplot color
var2 <- sort(unique(plotDF$var))
varSel2 <- c(combPFAS,"sumPFAS","log_sumPFAS")
#varSel2 <- allPFAS[grep("Log",allPFAS)]
#varSel2 <- "Lead"
plotDF2 <- plotDF%>%filter(var%in%varSel2)

p0 <- ggplot(plotDF2, aes(x=expected, y=observed, color=cluster))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    facet_wrap(var~treats, scales="free_y", ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          axis.title=element_text(size=12),
          axis.text=element_text(size=10),
          strip.text=element_text(size=12))

figfn <- paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".allPFAS_qqplot.png")
ggsave(figfn, p0, width=1200, height=500, units="px", dpi=120)







ci=0.95
deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
pc_results_bp <- deseqres %>% filter(!is.na(pvalue)) %>%
            arrange(pvalue) %>%
            mutate(r=rank(pvalue, ties.method = "random"),
                   pexp=r/length(pvalue),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(pvalue)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(pvalue)-r)))

    png(paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,"_qqplot.png"), width=1000, height=1000, res=120)
        p1 <- ggplot(pc_results_bp, aes(x=-log10(pexp),y=-log10(pvalue))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point(aes(color=cluster)) +
            #facet_grid(.~var)+
            geom_abline(slope=1,intercept=0) +
        ##    facet_grid(Origin ~ Location) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0(var," ", i," eGene QQ Plot")) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()


###############################################

i="RNA-CTRL"
    run="SES_PCs_sex_age_and_treats_generem"
    varrun=c(psychvarstorun,"age","factor_HS_CRP")
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
    subvars <- ldply(lapply(varrun, function(var){
        if(file.exists(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))){

        cat("running",i,var,"\n")
        deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
    }
    }),data.frame)
subvars_sig <- subset(subvars, padj<0.1)
length(unique(subvars_sig$identifier))

filenames <- c(paste0(project,".deseqres_SES-RNA-CTRL.txt"),paste0(project,".deseqres_SES-RNA-LPS.txt"))
data_names <- gsub(".txt", "", filenames) #remove file ending
shortnames <- gsub(".deseqres_", "", data_names)
shortnames <- gsub(project, "", shortnames)
for(i in 1:length(filenames)) assign(shortnames[i], fread(file.path(outFolder, filenames[i]),header = TRUE, sep='\t')[,analysis:=shortnames[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(shortnames, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(shortnames)
df <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) 




