# Specify variables to consider
# Age is continuous so we model it as a fixed effect
# Individual and Tissue are both categorical, so we model them as random effects
form <- ~ Age + (1 | Individual) + (1 | Tissue)
# Fit model
varPart <- fitExtractVarPartModel(geneExpr[1:5, ], form, info)
# convert to matrix
as.data.frame(varPart)


# Linear mixed model
fit <- lmer(geneExpr[1, ] ~ (1 | Tissue) + Age, info)
calcVarPart(fit)
# Linear model
# Note that the two models produce slightly different results
# This is expected: they are different statistical estimates
# of the same underlying value
fit <- lm(geneExpr[1, ] ~ Tissue + Age, info)
calcVarPart(fit)

library(variancePartition)
library(data.table)
library(edgeR)
library(plyr)
future::plan(strategy = 'multicore', workers = 5)
options(future.globals.maxSize = 30 * 1024 ^ 3)

args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_converted_dbgapID_key_n183_AR_18newvar_updated_08-05.txt","ALL","fastdemux") #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
run="SES_PCs_SES_sex_age_and_treats_adjusted"

outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

# set new output dir for filtered out unmatched figures
figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)

#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip","SNI_NoP") #SNI_NoP is the only variable that should be excluded based on the observed issues in score distributions that don’t make sense (negative values and extreme outliers)
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]
focus_vars <-c("factor_HS_CRP","HS_CRP","PSQI_total","HVS_mean","EDS_mean","SES","LogCRP","cytocomp","SLS","PSS_all_mean","DED_all_mean","Chol_HDL","Chol_LDL","BPs_avg","BPd_avg")


cluster="C0"
i="RNA-CTRL"
var="SES"
combat_run="SES_PCs_sex_age_and_treats_adjusted"
opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".DESeq_countlists.RData")
load(opfn)
all_metadata <- ldply(lapply(metadata_ls,function(i){
	i$rows <- rownames(i)
	return(i)
	}), data.frame)[,-1]
rownames(all_metadata) <- all_metadata$rows

all_counts_ls <- lapply(unique(all_metadata$letter_clusters),function(cluster){
	#cluster=unique(all_metadata$letter_clusters)[1]
	opfn <- paste0(outFolder,"combat_sessexagetreatpcs/",project,".ComBat_seq.",cluster,".",combat_run,".RData")
	load(opfn)
	genes <- rownames(adjusted_counts)
	adjusted_counts <- cbind(genes, as.data.frame(adjusted_counts))

	return(adjusted_counts)
})
all_counts <- ldply(all_counts_ls,data.frame)
rownames(all_counts) <- make.names(all_counts$genes, unique = TRUE)
all_counts[is.na(all_counts)] <- 0
all_counts <- as.matrix(all_counts[,-1])
colnames(all_counts) <- gsub("[.]","-",colnames(all_counts))

    #cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    #cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = "RNA-CTRL"))
    #cluster_metadata_t <- subset(cluster_metadata, treats==i)
    cluster_metadata_var <- all_metadata[,c("letter_clusters","Sample_ID","PC1","PC2","sex_alph","age","treats","SES",var)]
 
    cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(all_counts)),]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- all_counts[,which(colnames(all_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

#so anthony used plate id here (supposed to pick top source of variation), 
#I was going to use batch and then realized that is corrected for in combat. 
#so updatting to not subset treatment. could alternatively use sample id?
design_formula <-  ~ treats + letter_clusters
form <-  ~ (1|treats) + (1|letter_clusters) + PC1 + PC2 + (1|sex_alph) + age + SES 
genes = DGEList( cluster_counts_t ) 
genes = calcNormFactors( genes)
design = model.matrix( design_formula , cluster_metadata_var)
vobj = voom( genes, design, plot=FALSE)
#form = ~ (1|Region:bbscore) + (1|individualIdentifier) + (1|Run) + (1|LIMS) + (1|Region) + (1|`Sex (final)`) + RIN + rRNA.rate + PMI + Age + mappingRate + Apo1 + Apo2 
vp = fitExtractVarPartModel( vobj, form, cluster_metadata_var)
form <- formula(paste('~ (1|treats)','(1|letter_clusters)','PC1','PC2','(1|sex_alph)','age',var, sep=" + "))
vp = fitExtractVarPartModel( vobj, form, cluster_metadata_var)


#how results look
           letter_clusters     sex_alph       treats         PC1         PC2
AL627309.1       0.5298252 3.769260e-17 0.0005242804 0.003691427 0.004324701
AL627309.5       0.5256038 3.739229e-17 0.0083216504 0.003801510 0.004342489
LINC01409        0.6763228 1.605144e-18 0.0041471979 0.002352100 0.002657541
LINC01128        0.6757485 2.549730e-17 0.0067617490 0.002485410 0.002743661
LINC00115        0.5264650 3.745356e-17 0.0083352856 0.003080519 0.003443736
FAM41C           0.5305042 3.774092e-17 0.0005249523 0.003076869 0.003720078
                    age          SES         isel Residuals
AL627309.1 1.174857e-05 6.918563e-07 6.889157e-05 0.4615531
AL627309.5 1.485839e-05 1.759004e-07 3.991157e-05 0.4578757
LINC01409  3.634262e-06 1.526560e-06 3.057309e-05 0.3144847
LINC01128  1.346118e-06 2.659357e-06 3.740307e-05 0.3122193
LINC00115  8.980227e-08 3.033465e-06 4.648228e-05 0.4586259
FAM41C     1.491856e-08 4.934233e-07 2.868609e-05 0.4621447

opfn <- paste0(outFolder,project,".",var,".variancepartition.",combat_run,".RData")
save(vp, file=opfn)

png(width = 8, height = 8, file=paste0(figuredir,project,".variancepartition-",var,".",run,".png"), pointsize=12, 
    bg = "transparent", units = "in", res = 1200)
plotVarPart(sortCols(vp))
dev.off()

#testing a run with the interaction
form <- formula(paste('~ (1|treats)','(1|letter_clusters)','(1|letter_clusters:treats)','PC1','PC2','(1|sex_alph)','age',var, sep=" + "))
vp = fitExtractVarPartModel( vobj, form, cluster_metadata_var)
opfn <- paste0(outFolder,project,".",var,".variancepartition_withInteraction.",combat_run,".RData")
save(vp, file=opfn)

png(width = 8, height = 8, file=paste0(figuredir,project,".variancepartition_withInteraction-",var,".",run,".png"), pointsize=12, 
    bg = "transparent", units = "in", res = 1200)
plotVarPart(sortCols(vp))
dev.off()

#I couldn't find anything on github doing variance partitioning so I looked up methods and found the package variancePartition 
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1323-z http://bioconductor.org/packages/variancePartition . 
#I believe this does exactly what we want but I wanted to check a few things with you about the input data and the model since this is my first go. 
#I currently set it up to run similarly to how we have DESeq, with subsetting the input data to run per cluster/treatment/variable (off of COMBAT adjusted counts). Then the model is the same as what we used for DESeq (ex/ ~ PC1 + PC2 + sex_alph + age + SES + isel). 
#Any thoughts for change or improvement? 

#trial running per cluster

lapply(unique(all_metadata$letter_clusters),function(cluster){
	#cluster=unique(all_metadata$letter_clusters)[1]
	cat("running ", cluster, "\n")
	opfn <- paste0(outFolder,"combat_sessexagetreatpcs/",project,".ComBat_seq.",cluster,".",combat_run,".RData")
	load(opfn)
	metadata <- metadata_ls[[cluster]]
	cluster_metadata_var <- metadata[,c("Sample_ID","PC1","PC2","sex_alph","age","treats","SES",var)]
    cluster_metadata_var <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(adjusted_counts)),]
    cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
    cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
    all(colnames(cluster_counts_t) == rownames(cluster_metadata_var))

	design_formula <- ~ treats
	form <- formula(paste('~ (1|treats)','PC1','PC2','(1|sex_alph)','age',var, sep=" + "))

	genes = DGEList( cluster_counts_t ) 
	genes = calcNormFactors( genes)
	design = model.matrix( design_formula , cluster_metadata_var)
	vobj = voom( genes, design, plot=FALSE)
	#form = ~ (1|Region:bbscore) + (1|individualIdentifier) + (1|Run) + (1|LIMS) + (1|Region) + (1|`Sex (final)`) + RIN + rRNA.rate + PMI + Age + mappingRate + Apo1 + Apo2 
	vp = fitExtractVarPartModel( vobj, form, cluster_metadata_var)
	opfn <- paste0(outFolder,project,".",var,".",cluster,".variancepartition.",combat_run,".RData")
	save(vp, file=opfn)

	p <- plotVarPart(sortCols(vp))
	png(width = 8, height = 8, file=paste0(figuredir,project,".variancepartition-",var,".",cluster,".",run,".png"), pointsize=12, 
    bg = "transparent", units = "in", res = 1200)
	print(p)
	dev.off()
})
