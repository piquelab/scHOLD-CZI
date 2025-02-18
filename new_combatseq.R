require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)
require(BiocParallel)
library(Seurat)
library(scuttle)
library(tidyverse)
library(data.table)
library(plyr)
library(parallel)
library("AnnotationHub")
library(ggseurat)
library(cowplot)
library(sva)
library(ggpubr)
future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_chronic_cond_added_02_07_2024.txt","ALL","fastdemux",11,0.2) #for testing
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
#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o[,-c("sex","sex_alph","age")],cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip","SNI_NoP","age","sex","sex_alph") #SNI_NoP is the only variable that should be excluded based on the observed issues in score distributions that don’t make sense (negative values and extreme outliers)
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
focus_vars <-c("factor_HS_CRP","HS_CRP","PSQI_total","HVS_mean","EDS_mean","SES","LogCRP","cytocomp","SLS","PSS_all_mean","DED_all_mean","Chol_HDL","Chol_LDL","BPs_avg","BPd_avg")

# Variables to focus on
variablesL <- c("DSES_09","chronic_sum","SES","cytocomp","PSS_all_mean", "Chol_HDL", "BPd_avg", "isel", "pr_comp") 

variables <- c("DSES_07","DSES_09","FCDEM_12","chronic_sum","smoke","CVDRISK",
                "HS_CRP","il6","PSQI_total","HVS_mean","EDS_mean","SES","LogCRP","Logifny","Logil10","Logil12",
                "Logil13","Logil1b","Logil2","Logil4","Logil6","Logil8","Logtnfa","cytocomp","NAIscr2019","NAIscr2017",
                "LQ2019","StressSev","StressCount","nii_mean","NII_fam","SLS","LivingAlone","SNI_HCG","SNI_NoP","PSS_all_mean",
                "DED_all_mean","Chol","Trig","Chol_HDL","Chol_Ratio","Chol_LDL","HbA1C","BPs_avg","BPd_avg","isel", "pr_comp",
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
variables_df2 <- subset(variables_df, !variable %in% c("SNI_NoP"))

additionalvars <- fread("/rs/rs_grp_schold/covariates/other_covariates/chronic_chond_lables_attributes_varnames.txt")

variables_df <- merge(variables_df2,additionalvars[!additionalvars$Variable %in% c("HS_CRP"),-1],by.x=c("variable","description"),by.y=c("Variable","Column_Label"),all=T)

cat_cov <- fread(paste0(base,"categories_cov.txt"))
cat_cov_m <- merge(variables_df,cat_cov[,-2],by.y="Column_Name",by.x="variable",all.y=T)

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

sc@meta.data$letter_clusters <- paste0("C",sc@meta.data$seurat_clusters)
colOrd <- data.frame(NEW_BARCODE=names(sc$orig.ident))
coldata_ex <- merge(as.data.frame(sc@meta.data),eigenvec2,by="Sample_ID") #first 2 columns are dbgapid ie. sampleid and batch
RcolData <- merge(colOrd, coldata_ex, by="NEW_BARCODE",all.x=T)
RcolData <- RcolData[match(rownames(sc@meta.data), RcolData$NEW_BARCODE), ]
sc@meta.data <-cbind(sc@meta.data,RcolData[,which(!colnames(RcolData) %in% colnames(sc@meta.data))])

sc[["sex_male"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.male$symbol])
sc[["sex_female"]] <- PercentageFeatureSet(sc, features = rownames(sc)[rownames(sc) %in% geneIDs.female$symbol])

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
lowcell <- subset(cellcount, Freq<3000)
if(dim(lowcell)[1]==0){cat( "no low cell counts")}else{sc <-subset(x = sc, subset = letter_clusters %in% unique(lowcell$Var1), invert = TRUE)}

sce <- as.SingleCellExperiment(sc)
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".SingleCellExperiment.RDS")
write_rds(sce, opfn)

#seurat_clusters, treats, BATCH, Library
#/ aggregate by cluster,library info and covariate of interest:
sum_by <- c("letter_clusters", "BATCH","Sample_ID","treats")
summed <- aggregateAcrossCells(sce, id=colData(sce)[,sum_by])

#/ add rownames using the information from the colData:
colnames(summed) <- apply(colData(summed)[,sum_by], 1, function(x) paste(x, collapse="_"))
raw <- assay(summed, "counts")

counts_ls <- lapply(unique(summed$letter_clusters), function(i){
    cat("running ",i)
    #i <- unique(summed$letter_clusters)[1]
    cell_idx <- which(tstrsplit(colnames(summed), "_")[[1]] == i)
    keep <- which(colnames(raw) %in% colnames(summed[, cell_idx]))
    data <- raw[, keep]
    #filtered_data <- data > 0 
    filtered_data <- data
    filtered_data[filtered_data < 0] <- NA
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
    df <- data.frame(cluster_sample_id = colnames(i)) ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$letter_clusters <- tstrsplit(df$cluster_sample_id, "_")[[1]]
    df$BATCH  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
    df$Sample_ID  <- tstrsplit(df$cluster_sample_id, "_")[[3]]
    df$treats  <- tstrsplit(df$cluster_sample_id, "_")[[4]]
    metadata <- as.data.frame(colData(sce))
    test <- metadata %>% dplyr::select(-c(sex_male,sex_female,orig.ident,NEW_BARCODE,percent.mt,nCount_RNA,nFeature_RNA,NUM.READS,NUM.SNPS))
    test <- dplyr::select(test, -contains("RNA_snn_res"))
    df_n <- merge(df, unique(test), all.x=T)
    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df_n$letter_clusters))
    cell_counts <- t[, idx]
    ## Remove samples with zero cell contributing to the cluster
    cell_counts <- cell_counts[cell_counts > 0]
    ## Match order of cell_counts and sample_ids
    sample_order <- match(df_n$Sample_ID, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    ## Append cell_counts to data frame
    df_n$cell_count <- cell_counts
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df_n) <- df_n$cluster_sample_id
    #df_n <- df_n[complete.cases(df_n), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    return(df_n)
})
all(names(counts_ls) == names(metadata_ls))

#opfn <- paste0(outFolder,project,".DESeq_countlists.RData")
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.RData")
save(counts_ls,metadata_ls, file=opfn)
#load(opfn)

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

SEScov=TRUE
if(SEScov){
    run="SES_PCs_SES_sex_age_and_treats_generem"
    varrun=c(psychvarstorun[!psychvarstorun=="SES"],"age","factor_HS_CRP")
} else{
    run="SES_PCs_sex_age_and_treats_generem"
    varrun=c(psychvarstorun,"age","factor_HS_CRP")
}
combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sigDEGs/"))) dir.create(paste0(outFolder,"sigDEGs/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"deseqres/"))) dir.create(paste0(outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid/"))) dir.create(paste0(outFolder,"sampleid/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"sampleid_all/"))) dir.create(paste0(outFolder,"sampleid_all/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

firstrunvars=c("SES","pr_comp","isel","PSS_all_mean","BPd_avg")
updated_vars <- c("SNI_NumPeople_r", "chronic_sum_categ", "nii_mean", "LogDED", "LogTrig")

lapply(names(counts_ls),function(cluster){
    #cluster=names(counts_ls)[1]
        cat("running ", cluster, "\n")
    fdr=0.05
    cluster_counts_sce <- counts_ls[[cluster]]
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_counts <- assay(cluster_counts_sce, "counts")
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
	opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)

    d <- data.frame(number_ind=length(unique(cluster_metadata[,c("Sample_ID")])))
    fwrite(d, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid_all/",project,".",cluster,".sampleid_all.",run,".txt"))

    for(i in unique(cluster_metadata_sce$treats)){
            #lapply(list.df, subset, B!=2)
            #i <- "RNA-CTRL"
        cluster_metadata_t <- subset(cluster_metadata, treats==i)
        cluster_metadata_t <- transform(cluster_metadata_t, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata_t <- within(cluster_metadata_t, sex_alph <- relevel(sex_alph, ref = "Male"))

        lapply(c(variables_df$variable),function(var){
            #c(psychvarstorun,"factor_HS_CRP",updated_vars)
            #var="FCDEM_11_8"
            if(!isTRUE(file.size(paste0(outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",var,"-",i,".",run,".txt")) > 0)){
                cat("running deseq ",var,i," \n")
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + ",var)
        if(SEScov){
        design <-  paste0("~ PC1 + PC2 + sex_alph + age + SES + ",var)
        }
        if(var=="CVDRISK"){
        design <-  paste0("~ PC1 + PC2 + ",var)
        if(SEScov){
        design <-  paste0("~ PC1 + PC2 + SES + ",var)
        }
        }
        if(var=="factor_HS_CRP"){
        cluster_metadata_t <- subset(cluster_metadata_t, HS_CRP<10) #advised to remove as likely an infection
        cluster_metadata_t <- transform(cluster_metadata_t, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
        cluster_metadata_t <- within(cluster_metadata_t, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
        }
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
        if(SEScov){
        cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","SES",var)]
        }
        if(var=="age"){
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
        dds <- DESeq(dds,parallel=TRUE)
        opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",i,"-",var,cluster,".",run,".RDS")
        saveRDS(dds, file=opfn)

        fwrite(data.frame(Sample_ID=cluster_metadata_var[,c("Sample_ID")]), sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sampleid/",project,".",resset,".",dimset,".",cluster,".sampleid_",var,"-",i,".",run,".txt"))

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
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs_",var,"-",i,".",run,".txt"))
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
    #for (var in c("factor_HS_CRP","CVDRISK")){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
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

library(flextable)
library(ftExtra)
library(rlist)
library(officer)
#[!psychvarstorun %in% puberty]
 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
subvars <- ldply(lapply(c(variables_df$variable), function(var){
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
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,topvar=ifelse(var %in% variablesL,"Y",""),list.cbind(newls))
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
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"sigDEGs/",project,".",resset,".",dimset,".",cluster,".sigDEGs.",contrast,".",run,".txt"))
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

#plotting LPS vs CONTROL
run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
threshold=0.10
variablesL <- variablesL[-c(1,3)]


treat="RNA-LPS"
control="RNA-CTRL"
subvars <- ldply(lapply(variablesL, function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",control,".",run,".txt"))
        treat_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",treat,".",run,".txt"))
        merged_treat <- merge(ctrl_deseqres,treat_deseqres,by=c("identifier","cluster","var"))
        merged_treat <- transform(merged_treat, CTRL_zscore=logFC.x/SE.x, LPS_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3CTRL_sig", ifelse(padj.y<=threshold, "2LPS_sig", "1Not_Sig")))))
        return(merged_treat)
}),data.frame)
#didnt update this one for the colored by sig
p <- ggplot(subvars, aes(x=CTRL_zscore, y=LPS_zscore)) +
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
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()


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
