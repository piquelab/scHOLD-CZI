library(data.table)
library(plyr);library(dplyr)
library("limma")
library(edgeR)
library(ggplot2)
library(flextable)
library(ftExtra)
library(rlist)
library(officer)
 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

##########################

#job="ALOFT"
job="CZI"

if(job=="ALOFT"){
args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_uniq-n265_psesl-a2_fixed_12-19-2024.txt","ALL","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt",0.2) 
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
sample_batch <- args[5]
resset <- args[6]
dimset=50
outFolder=paste0(base,method,"_pseudobulk_ctrl/")
#for resolution 0.2 and including V2 chem
outFolder=paste0(outFolder,"lessfilt/")
figuredir=paste0(outFolder,"figures/")
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
exp <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt")
eigenvec2 <- merge(exp,cov_file,by="dbgap.ID",all.x=T)
eigenvec2 <- eigenvec2 %>% dplyr::filter(SCAIP7_18 == 1 | SCAIP1_6 == 1)
eigenvec2$Sample_ID <- eigenvec2$dbgap.ID
eigenvec2 <- transform(eigenvec2, SCAIP1_6=ifelse(is.na(SCAIP1_6),0,SCAIP1_6)) #weird case where some are NA ... Ali not sure why
psychvarstorun <- colnames(eigenvec2[,24:59])
treatments=c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")
treatmentsfirst=c("CTRL","PHA")
contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA-DEX")) #no lps vs lps-dex as too few ind
opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)
combatrun="income_PCs_sex_age_and_treats_adjusted"

puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
allvars <- psychvarstorun[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]
firstrunvars <- allvars[c(1:10)]
secondrunvars <- c(allvars[c(1:20)],"cage1")
} else if(job=="CZI"){
proteincoding=FALSE
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.15) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
contrastdf <- data.frame(control=c("RNA-CTRL"),treatment=c("RNA-LPS"))
#read in genotype PC (run only on current samples. if adding data, since I made the file 03/20/24 remake using plink_to_PC.R)
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
if(proteincoding){
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt_proteincoding/")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt_proteincoding.RData")
#load(opfn)
} else{
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt.RData")
load(opfn)
}
#counts_ls$C6 <- NULL
#metadata_ls$C6 <- NULL
combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
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
allPFAS <- psychvarstorun[c(33:53)]
combPFAS <- allPFAS[seq(1,21,3)]
PFASvars <- c(allPFAS,"sumPFAS","log_sumPFAS")
allvars <- c(reordered_psychvarstorun[c(1:82)],"age","Lead","sumPFAS","log_sumPFAS",zcytokines)
}
contrastdf <- transform(contrastdf, contrast=paste0(treatment,"_vs_",control))

baserun="treatvarint_withCOMBAT_limma"
log2cpm=FALSE
voom=TRUE
cpmqqnorm=FALSE
iselcov=FALSE
WHRcov=FALSE
runPFAS=TRUE
if(voom ){
    run=paste0(baserun,"_voom")
} else if (cpmqqnorm){
    run=paste0(baserun,"_cpmqqnorm")
} else{
    run=baserun
}
if(iselcov){
    run=paste0(run,"_iselcov")
}
if(WHRcov){
    run=paste0(run,"_WHRcov")
}
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"limmares/"))) dir.create(paste0(outFolder,"limmares/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

#create column headers for summary table to append to during loop -- do just once
if (!file.exists(paste0(outFolder,"summary.",run,".txt"))) {
    cat("making file \n")
table <- data.frame(cluster=NA, contrast=NA, variable=NA,number_individuals=NA,gene_number=NA,intDEGs=NA,treatDEGs=NA,varDEGs=NA)
fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"summary.",run,".txt"))
}

lapply(names(counts_ls),function(c){
    #c="C6"
    cluster_metadata_sce <- metadata_ls[[c]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",c,".",combatrun,".RData")
    load(opfn)
    if(job=="ALOFT"){
        adjusted_counts <- adjusted
            cluster_metadata <- transform(cluster_metadata, Sex=as.factor(Sex))
        if(all(c("0","1") %in% levels(cluster_metadata$Sex))){
            cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "0"))
        }
    } else if(job=="CZI"){
        cluster_metadata <- transform(cluster_metadata, sex_alph=as.factor(sex_alph),SNI_NumPeople_r=as.numeric(SNI_NumPeople_r))
        cluster_metadata <- within(cluster_metadata, sex_alph <- relevel(sex_alph, ref = "Male"))
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
        #keep <- rowSums(adjusted_countsP,na.rm=T) >= (ncol(adjusted_countsP) / 2) #filter to keep genes expressed in 25% of samples
        #adjusted_counts <- unlist(adjusted_countsP[keep, ])
        }
    }
    lapply(allvars,function(var){
        #var="pedu"
        if(var=="factor_HS_CRP"){
        cluster_metadata <- subset(cluster_metadata, HS_CRP<10) #advised to remove as likely an infection
        cluster_metadata <- transform(cluster_metadata, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
        cluster_metadata <- within(cluster_metadata, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
        }
        cat("running ", c, var, "\n")
        lapply(1:length(contrastdf$control),function(x){
            con=contrastdf[x,]
            contrast=paste0(con$treatment,"_vs_",con$control)
            cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
            #i=treatments[1]
            table <- fread(paste0(outFolder,"summary.",run,".txt"))
            currenttable <- subset(table, variable==var & contrast==paste0(con$treatment,"_vs_",con$control) & cluster==c)
            #this runs only if the summary table doesnt have it already
            if(dim(currenttable)[1]<1){
                cat("running",con$treatment,"_vs_",con$control,"\n")
            cluster_metadata_t <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
            if(job=="ALOFT"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]
            } else if (job=="CZI"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","treats",var)]
            }
            if(var=="age" | var=="cage1"){
            if(job=="ALOFT"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","treats",var)]
            } else if (job=="CZI"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","treats",var)]
            }}
            if(iselcov){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","ISEL_Mean","treats",var)]
            }
            if(WHRcov){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","WHR","treats",var)]
            }
            if(var=="CVDRISK"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","treats",var)]
            }
            cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
            cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
            cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
            cv_d <- transform(cv_d, treats=factor(treats))
            ## Normalization of data
            # make edgeR object
            dge <- DGEList(counts=cluster_counts_t)
            #Transform counts to counts per million
            dge <- calcNormFactors(dge)
            varcol=cv_d[,var]
            if(job=="ALOFT"){
            design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
            if(var=="cage1"){
            design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.factor(cv_d$Sex) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
            }} 
            if (job=="CZI" & !var=="CVDRISK"){
            design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            if(var=="age"){
            design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.factor(cv_d$sex_alph) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            }}
            if(var=="CVDRISK"){
            design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            }
            if(iselcov){
            design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.numeric(cv_d$ISEL_Mean) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            }
            if(WHRcov){
            design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.numeric(cv_d$WHR) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            }
            if(log2cpm){
            cpm <- log2(cpm(dge)+1)
            rownames(cpm) <- rownames(cluster_counts_t)
            fit <- lmFit(cpm, design)
            }
            if(voom){
            v_e <- voom(dge, design, plot=FALSE, normalize.method="quantile")
            #genes_normed_baseline <- data.frame(v_e$E)
            fit <- lmFit(v_e, design)
            }
            if(cpmqqnorm){
            cpm <- cpm(dge)
            rownames(cpm) <- rownames(cluster_counts_t)
            mat_qnorm <- t(apply(cpm,1,function(x){qqnorm(rank(x, ties.method = "random"), plot = F)$x}))
            colnames(mat_qnorm) <- colnames(cpm)
            fit <- lmFit(mat_qnorm, design)
            }
            fit2 <- eBayes(fit)
            tt <- topTable(fit2, n = Inf, adjust.method = "BH", coef = colnames(fit2)[grepl("varcol", colnames(fit2)) & grepl("treats", colnames(fit2))])
            tt <- transform(tt, cluster=c, contrast=contrast, variable=var, term="interaction", identifier=rownames(tt))
            intdf <- tt[,c("cluster","contrast","variable","term","identifier","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
            tt <- topTable(fit2, n = Inf, adjust.method = "BH", coef = colnames(fit2)[!grepl("varcol", colnames(fit2)) & grepl("treats", colnames(fit2))])
            tt <- transform(tt, cluster=c, contrast=contrast, variable=var, term="treat", identifier=rownames(tt))
            treatdf <- tt[,c("cluster","contrast","variable","term","identifier","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
            tt <- topTable(fit2, n = Inf, adjust.method = "BH", coef = colnames(fit2)[grepl("varcol", colnames(fit2)) & !grepl("treats", colnames(fit2))])
            tt <- transform(tt, cluster=c, contrast=contrast, variable=var, term="variable", identifier=rownames(tt))
            vardf <- tt[,c("cluster","contrast","variable","term","identifier","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
            fwrite(intdf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"limmares/",c,".",contrast,".",var,".lm_interaction.",run,".txt"))
            fwrite(treatdf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"limmares/",c,".",contrast,".",var,".lm_treat.",run,".txt"))
            fwrite(vardf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"limmares/",c,".",contrast,".",var,".lm_variable.",run,".txt"))
            table <- data.frame(cluster=c, contrast=contrast, variable=var)
            table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
            table$gene_number <- paste(nrow(cluster_counts_t))
            table$intDEGs=nrow(subset(intdf, adj.P.Val<0.1))
            table$treatDEGs=nrow(subset(treatdf, adj.P.Val<0.1))
            table$varDEGs=nrow(subset(vardf, adj.P.Val<0.1))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=F, append=T, paste0(outFolder,"summary.",run,".txt"))
            rm(currenttable,cluster_metadata_t,cluster_metadata_var,cv_d,dge,cpm,fit,fit2,tt,vardf,treatdf,intdf,table)
            gc(reset=T)
            }
        })
    })
})

     
myDir <- paste0(outFolder,"limmares/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grepl("lm_.*.treatvarint.txt",filenames)]
data_names <- gsub(".treatvarint.txt", "", filenames)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
res <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(res, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"all_limmares_treatinteraction.txt"))


lmcoefspha <- subset(res, contrast=="PHA_vs_CTRL" & variable=="pnsi")
lmcoefsphaccl20 <- subset(lmcoefspha, identifier=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, identifier=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, identifier=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, identifier=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, identifier=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, identifier=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)
fwrite(alltested, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"check.treatinteraction.txt"))


#New table: var|#CTRL|#LPS|#interactiontest for 
varrun="limma_var_voom"
varoutFolder=paste0(baseoutFolder,run,"/")


run="SES_PCs_sex_age_and_treats_generem"
#run="SES_PCs_sex_age_and_treats_generem_zingeR"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
variablesL <- c("chronic_sum","cytocomp", "SES","PSS_all_mean", "Chol_HDL", "BPd_avg", "isel", "pr_comp")  #for regular runs
#secondrunvars=c(psychvarstorun[c(1:13,24:45)],"age")
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
#intrun="SES_PCs_sex_age_and_treats_generem_zingeR_treatint"
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
treat="RNA-LPS"
control="RNA-CTRL"
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































