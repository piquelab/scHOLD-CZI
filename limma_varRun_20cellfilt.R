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

future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 30 * 1024 ^ 3)

library(biomaRt)  
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name","transcript_biotype"), filters = c("transcript_biotype","chromosome_name"),values = list("protein_coding",1:22), mart = mart)


#job="ALOFT"
job="CZI"
proteincoding=FALSE


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
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.15) #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]
contrastdf <- data.frame(control=c("RNA-CTRL","RNA-LPS"),treatment=c("RNA-LPS","RNA-LPS-DEX"))
outFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
figuredir=paste0(outFolder,"figures/")
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
opfn <- paste0(base,method,"_pseudobulk_ctrl/nodex/",project,".",resset,".",dimset,".DESeq_countlists.bticfilt_proteincoding.RData")
load(opfn)
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

baserun="limma_var"
log2cpm=FALSE
voom=TRUE
cpmqqnorm=FALSE
iselcov=FALSE
WHRcov=FALSE
runPFAS=TRUE
#proteincoding33=FALSE
#proteincoding50=FALSE
thirtythree=FALSE
fifty=FALSE
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
if(thirtythree){
    run=paste0(run,"_33")
}
if(fifty){
    run=paste0(run,"_50")
}
outFolder=paste0(baseoutFolder,run,"/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
if (!file.exists(paste0(outFolder,"stats/"))) dir.create(paste0(outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"limmares/"))) dir.create(paste0(outFolder,"limmares/"), showWarnings=F)
if (!file.exists(paste0(outFolder,"figures/"))) dir.create(paste0(outFolder,"figures/"), showWarnings=F)

#create column headers for summary table to append to during loop -- do just once
if (!file.exists(paste0(outFolder,"summary.",run,".txt"))) {
    cat("making file \n")
table <- data.frame(cluster=NA, treat=NA, variable=NA,number_individuals=NA,gene_number=NA,tested_genes=NA,varDEGs=NA)
fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"summary.",run,".txt"))
}

#-1 is just for CZI until I figure out what happened to C7
lapply(names(counts_ls)[c(1:length(names(counts_ls)))],function(c){
    #c="C0"
    cluster_metadata_sce <- metadata_ls[[c]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
    #if(job=="CZI"){
    #cluster_metadataL <- left_join(cluster_metadata,unique(eigenvec2[,c("Sample_ID","Lead","WHR")]),by="Sample_ID")
    #rownames(cluster_metadataL) <- cluster_metadataL$rowid
    #cluster_metadata <- cluster_metadataL
    #}
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
    }
    }
    #adjusted_countsP <- adjusted_counts[rownames(adjusted_counts) %in% genes$hgnc_symbol, ]
    filtered_data <- adjusted_counts > 0 
    if( thirtythree){
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(adjusted_counts) / 3) #filter to keep genes expressed in 33% of samples
    adjusted_counts <- unlist(adjusted_counts[keep, ])
    }
    if( fifty){
    keep <- rowSums(filtered_data,na.rm=T) >= (ncol(adjusted_counts) / 2) #filter to keep genes expressed in 33% of samples
    adjusted_counts <- unlist(adjusted_counts[keep, ])
    }
    lapply(allvars,function(var){
        #var="FCDEM_11_13"
        if(var=="factor_HS_CRP"){
        cluster_metadata <- subset(cluster_metadata, HS_CRP<10) #advised to remove as likely an infection
        cluster_metadata <- transform(cluster_metadata, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
        cluster_metadata <- within(cluster_metadata, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
        }
        cat("running ", c, var, "\n")
        lapply(treatmentsfirst,function(i){
            #i=treatments[1]
            table <- fread(paste0(outFolder,"summary.",run,".txt"))
            currenttable <- subset(table, variable==var & treat==i & cluster==c)
            #this runs only if the summary table doesnt have it already
            if(dim(currenttable)[1]<1){
                cat("running",i,"\n")
            cluster_metadata_t <- subset(cluster_metadata, treats==i)
            if(job=="ALOFT"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1",var)]
            } else if (job=="CZI"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",var)]
            }
            if(var=="age" | var=="cage1"){
            if(job=="ALOFT"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex",var)]
            } else if (job=="CZI"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph",var)]
            }}
            if(iselcov){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","ISEL_Mean",var)]
            }
            if(WHRcov){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","WHR",var)]
            }
            if(var=="CVDRISK"){
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2",var)]
            }
            cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
            cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
            cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
            ## Normalization of data
            # make edgeR object
            dge <- DGEList(counts=cluster_counts_t)
            #Transform counts to counts per million
            dge <- calcNormFactors(dge)
            varcol=cv_d[,var]
            if(job=="ALOFT"){
            design <- model.matrix(~ 0 + as.numeric(varcol) + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
            if(var=="cage1"){
            design <- model.matrix(~ 0 + as.numeric(varcol) + as.factor(cv_d$Sex) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
            }} 
            if (job=="CZI" & !var=="CVDRISK"){
            design <- model.matrix(~ 0 + as.numeric(varcol) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            if(var=="age"){
            design <- model.matrix(~ 0 + as.numeric(varcol) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            }}
            if(var=="CVDRISK"){
            design <- model.matrix(~ 0 + as.numeric(varcol) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            }
            if(iselcov){
            design <- model.matrix(~ 0 + as.numeric(cv_d$ISEL_Mean) + as.numeric(varcol) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
            }
            if(WHRcov){
            design <- model.matrix(~ 0 + as.numeric(cv_d$WHR) + as.numeric(varcol) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
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
            opfn <- paste0(outFolder,"limma_fit-",i,"-",var,c,".",run,".RData")
            save(fit, file=opfn)

            fit2 <- eBayes(fit)
            tt <- topTable(fit2, n = Inf, adjust.method = "BH", coef = colnames(fit2)[grepl("varcol", colnames(fit2)) & !grepl("treats", colnames(fit2))])
            tt <- transform(tt, cluster=c, treat=i, variable=var, identifier=rownames(tt))
            vardf <- tt[,c("cluster","treat","variable","identifier","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
            fwrite(vardf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"limmares/",c,".",i,".",var,".",run,".txt"))
            table <- data.frame(cluster=c, treat=i, variable=var)
            table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
            table$gene_number <- paste(nrow(cluster_counts_t))
            table$tested_genes <- nrow(subset(vardf, !is.na(adj.P.Val)))
            table$varDEGs=nrow(subset(vardf, adj.P.Val<0.1))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=F, append=T, paste0(outFolder,"summary.",run,".txt"))
            rm(currenttable,cluster_metadata_t,cluster_metadata_var,cv_d,dge,cpm,fit,fit2,tt,vardf,table)
            gc(reset=T)
            }
        })
    })
})

#for sex separate variables (ie/puberty vars in aloft)
if(job=="ALOFT"){
lapply(names(counts_ls),function(c){
    #c="C0"
    cluster_metadata_sce <- metadata_ls[[c]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",c,".",combatrun,".RData")
    load(opfn)
    adjusted_counts <- adjusted
    cluster_metadata <- transform(cluster_metadata, Sex=as.factor(ifelse(Sex==0,"female","male")))
    if(all(c("male","female") %in% levels(cluster_metadata$Sex))){
        cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "male"))
    }
    for(s in c("male","female")){
        if (s=="female"){
        puberty_sex <- c("cgpd5","cgpd")
            }else {
        puberty_sex <- c("cbpd")
            }
    lapply(puberty_sex, function(var){
        cat("running ", c, s, var, "\n")
        lapply(treatmentsfirst,function(i){
            #i=treatmentsfirst[1]
            table <- fread(paste0(outFolder,"summary.",run,".txt"))
            currenttable <- subset(table, variable==var & treat==i & cluster==c)
            #this runs only if the summary table doesnt have it already
            if(dim(currenttable)[1]<1){
                cat("running",i,"\n")
            cluster_metadata_t <- subset(cluster_metadata, treats==i & Sex==s)
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3",var)]
            cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
            cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
            cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
            ## Normalization of data
            # make edgeR object
            dge <- DGEList(counts=cluster_counts_t)
            #Transform counts to counts per million
            dge <- calcNormFactors(dge)
            varcol=cv_d[,var]
            design <- model.matrix(~ 0 + as.numeric(varcol) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
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
            opfn <- paste0(outFolder,project,".",resset,".",dimset,".limma_fit-",i,"-",var,c,".",run,".RData")
            save(fit, file=opfn)

            fit2 <- eBayes(fit)
            tt <- topTable(fit2, n = Inf, adjust.method = "BH", coef = colnames(fit2)[grepl("varcol", colnames(fit2)) & !grepl("treats", colnames(fit2))])
            tt <- transform(tt, cluster=c, treat=i, variable=var, identifier=rownames(tt))
            vardf <- tt[,c("cluster","treat","variable","identifier","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
            fwrite(vardf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"limmares/",c,".",i,".",var,".",run,".txt"))
            table <- data.frame(cluster=c, treat=i, variable=var)
            table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
            table$gene_number <- paste(nrow(cluster_counts_t))
            table$tested_genes <- nrow(subset(vardf, !is.na(adj.P.Val)))
            table$varDEGs=nrow(subset(vardf, adj.P.Val<0.1))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=F, append=T, paste0(outFolder,"summary.",run,".txt"))
            rm(currenttable,cluster_metadata_t,cluster_metadata_var,cv_d,dge,cpm,fit,fit2,tt,vardf,table)
            gc(reset=T)
            } #currentable if
        }) #treatment loop
    }) #var loop
    } #sex loop
}) #cluster loop
#want to check age effects separately in male/female
#create column headers for summary table to append to during loop -- do just once
for(s in c("male","female")){
if (!file.exists(paste0(outFolder,"summary.",s, "." ,run,".txt"))) {
    cat("making file \n")
table <- data.frame(cluster=NA, treat=NA, variable=NA,sex=NA, number_individuals=NA,gene_number=NA,tested_genes=NA,varDEGs=NA)
fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"summary.",s, "." ,run,".txt"))
}
}
lapply(names(counts_ls),function(c){
    #c="C0"
    cluster_metadata_sce <- metadata_ls[[c]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",c,".",combatrun,".RData")
    load(opfn)
    adjusted_counts <- adjusted
    cluster_metadata <- transform(cluster_metadata, Sex=as.factor(ifelse(Sex==0,"female","male")))
    if(all(c("male","female") %in% levels(cluster_metadata$Sex))){
        cluster_metadata <- within(cluster_metadata, Sex <- relevel(Sex, ref = "male"))
    }
    for(s in c("male","female")){
    var="cage1"
        cat("running ", c, s, var, "\n")
        lapply(treatmentsfirst,function(i){
            #i=treatmentsfirst[1]
            table <- fread(paste0(outFolder,"summary.",s, "." ,run,".txt"))
            currenttable <- subset(table, variable==var & treat==i & cluster==c)
            #this runs only if the summary table doesnt have it already
            if(dim(currenttable)[1]<1){
                cat("running",i,"\n")
            cluster_metadata_t <- subset(cluster_metadata, treats==i & Sex==s)
            cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3",var)]
            cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
            cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
            cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
            ## Normalization of data
            # make edgeR object
            dge <- DGEList(counts=cluster_counts_t)
            #Transform counts to counts per million
            dge <- calcNormFactors(dge)
            varcol=cv_d[,var]
            design <- model.matrix(~ 0 + as.numeric(varcol) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
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
            tt <- topTable(fit2, n = Inf, adjust.method = "BH", coef = colnames(fit2)[grepl("varcol", colnames(fit2)) & !grepl("treats", colnames(fit2))])
            tt <- transform(tt, cluster=c, treat=i, variable=var, sex=s, identifier=rownames(tt))
            vardf <- tt[,c("cluster","treat","variable","identifier","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
            fwrite(vardf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"limmares/",c,".",i,".",var,".",s,".",run,".txt"))
            table <- data.frame(cluster=c, treat=i, variable=var, sex=s)
            table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
            table$gene_number <- paste(nrow(cluster_counts_t))
            table$tested_genes <- nrow(subset(vardf, !is.na(adj.P.Val)))
            table$varDEGs=nrow(subset(vardf, adj.P.Val<0.1))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=F, append=T, paste0(outFolder,"summary.",s, "." ,run,".txt"))
            rm(currenttable,cluster_metadata_t,cluster_metadata_var,cv_d,dge,cpm,fit,fit2,tt,vardf,table)
            gc(reset=T)
            } #currentable if
        }) #treatment loop
    } #sex loop
}) #cluster loop
} #if job


var="FCDEM_11_11"
png(width = 6, height = 6, file=paste0(outFolder,"figures/MDplot_",var,"_C0.png"), pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 500)
plotMD(fit,coef=1,status=decideTests(fit))
dev.off()

var="PFOA"
c="C4"
i="RNA-CTRL"
opfn <- paste0(outFolder,"limma_fit-",i,"-",var,c,".",run,".RData")
load(opfn)
png(width = 6, height = 6, file=paste0(outFolder,"figures/MDplot_",var,"_",c,".png"), pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 500)
plotMD(fit,coef=1,status=decideTests(fit),main=paste0(i," ",var," ",c))
dev.off()

#fwrite(lmcoefs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"summary.",run,".txt"))
#make sure summary table is unique
            table <- unique(fread(paste0(outFolder,"summary.",run,".txt")))
            fwrite(table, sep='\t', quote=F, col.names=T, row.names=F, paste0(outFolder,"summary.",run,".txt"))

#for CZI, this is too memory intensive
#myDir <- paste0(outFolder,"limmares/")
#filenames <- list.files(myDir) #file list from directory
#filenames <- filenames[grepl(paste0(run,".txt"),filenames)]
#data_names <- gsub(paste0(".",run,".txt"), "", filenames)
#for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
#ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
#names(ddf) <- c(data_names)
#res <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#fwrite(res, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"all_limmares.",run,".txt"))


if(job=="ALOFT"){
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
allvarsdf <- data.frame(variable=allvars)
variables_df <- unique(merge(variables_df,allvarsdf,by="variable",all=T))

} else if(job=="CZI"){
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
variables_df <- unique(merge(variables_df,allvarsdf,by="variable",all=T))
}

varsrun = subset(variables_df, variable %in% allvars)
statsfull <- unique(fread(paste0(outFolder,"summary.",run,".txt")))

lapply(treatmentsfirst, function(i) {
    lapply(varsrun$variable, function(var){
    cat("running",i,var,"\n")
    if(!file.exists(paste0(outFolder,"limmares.",i,".",var,".",run,".txt"))){
    stats <- subset(statsfull, treat==i & variable==var)
    #ressub <- subset(res, treat==i & variable==var)
    myDir <- paste0(outFolder,"limmares/")
    filenames <- list.files(myDir) #file list from directory
    filenames <- filenames[grepl(paste0(".",i,".",var,".",run,".txt"),filenames)]
    data_names <- gsub(paste0(".",run,".txt"), "", filenames)
    for(f in 1:length(filenames)) assign(data_names[f], fread(file.path(myDir, filenames[f]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
    ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
    names(ddf) <- c(data_names)
    ressub <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
    fwrite(ressub, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"limmares.",i,".",var,".",run,".txt"))
    }
    })
})

#to plot
zoom="allvars"
lapply(treatmentsfirst, function(i) {
subvars <- ldply(lapply(varsrun$variable, function(var){
    cat("running",i,var,"\n")
    stats <- subset(statsfull, treat==i & variable==var)
    #ressub <- subset(res, treat==i & variable==var)
    ressub <- fread(paste0(outFolder,"limmares.",i,".",var,".",run,".txt"))
    #first subset variables that have at least one cluster with 50DEGs
    if(dim(subset(stats, varDEGs>50))[1]<1){ 
        cat("not enough DEGs","\n")
        } else{
        resshort <- ressub[,c("cluster","identifier")]
        df <- ddply(resshort, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$varDEGs)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newls <- newls[sapply(newls, nrow)>0]
            newtable <- data.frame(variable=var, description=variables_df[variables_df$variable==var,]$description,list.cbind(newls))
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
})


#volcano
library(EnhancedVolcano)
var="PFOA"
c="C4"
i="RNA-CTRL"

limmares <- fread(paste0(outFolder,"limmares.",i,".",var,".",run,".txt"))
#for (c in names(counts_ls)){
    c="C4"
    cat("running ", c, "\n")
    sub.table <- subset(limmares, cluster==c)
    sub.table <- sub.table[order(sub.table$adj.P.Val,-abs(sub.table$logFC)), ]
    xmax=max(sub.table$logFC)+0.5
    xmin=min(sub.table$logFC)-0.5
    topp1 <- min(subset(sub.table, adj.P.Val > quantile(adj.P.Val, prob = 1 -99/100,na.rm=T))$P.Value)
    tolab <- c(unique(head(sub.table,n=20)$identifier))
    png(width = 8, height = 8, file=paste0(outFolder,"figures/dge_volcano-",var,"-",i,".",run,"-",c,".png"), pointsize=12, 
          bg = "transparent", units = "in", res = 1200)
    p <- EnhancedVolcano(sub.table,
      lab = sub.table$identifier,
      x = 'logFC',
      y = 'P.Value',
      title=paste0(c),
      subtitle = NULL,
      xlim = c(xmin, xmax),
      pCutoffCol= 'adj.P.Val',
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
#compare deseq and limma specific genes
#C10 PHA cdres condition many sig DEGs for deseq in aloft but none in limma, also small cluster
#C0 PHA mono_av condition many sig DEGs for limma and none for deseq, also large cluster
#C11 PHA mono_av condition segDEGs for both deseq and limma, also small cluster
library(ggpubr)
library(scales)

threshold=0.10
if(job=="ALOFT"){
deseqrun="income_PCs_sex_age_and_treats_adjusted_withWave"
deseqoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/",deseqrun,"/")
i="PHA"
} else if (job=="CZI"){
i=treatmentsfirst[1]
deseqrun="SES_PCs_sex_age_and_treats_generem"
deseqoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/",deseqrun,"/")
}
#res <- fread(paste0(outFolder,"all_limmares.",run,".txt"))
#names(res)[c(2,3)] <- c("treats","var")

#nothing merged for PHA and secondrunvars
  psy <- ldply(lapply(firstrunvars[!firstrunvars %in% c("cage1")], function(v){
    deseqres <- fread(paste0(deseqoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",deseqrun,".txt"))
    ressub <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
    names(ressub)[c(2,3)] <- c("treats","var")
    m <- merge(ressub[ressub$adj.P.Val<threshold,],deseqres[deseqres$padj<threshold,],by=c("identifier","cluster","treats","var"))
    return(m)
    }),data.frame)
psy <- psy[order(psy$padj,psy$adj.P.Val),]
deseqonly <- ldply(lapply(firstrunvars[!firstrunvars %in% c("cage1")], function(v){
    deseqres <- fread(paste0(deseqoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",deseqrun,".txt"))
    ressub <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
    names(ressub)[c(2,3)] <- c("treats","var")
    m <- merge(ressub[ressub$adj.P.Val>threshold,],deseqres[deseqres$padj<threshold,],by=c("identifier","cluster","treats","var"))
    return(m)
    }),data.frame)
deseqonly <- deseqonly[order(deseqonly$padj,-abs(deseqonly$logFC.y)),]
limmaonly <- ldply(lapply(firstrunvars[!firstrunvars %in% c("cage1")], function(v){
    deseqres <- fread(paste0(deseqoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",deseqrun,".txt"))
    ressub <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
    names(ressub)[c(2,3)] <- c("treats","var")
    m <- merge(ressub[ressub$adj.P.Val<threshold,],deseqres[deseqres$padj>threshold,],by=c("identifier","cluster","treats","var"))
    return(m)
    }),data.frame)
limmaonly <- limmaonly[order(limmaonly$adj.P.Val,-abs(limmaonly$logFC.x)),]
bothexamplenum=1
conditions <- data.frame(c=c(deseqonly[bothexamplenum,]$cluster,limmaonly[bothexamplenum,]$cluster,psy[bothexamplenum,]$cluster),
    i=c(deseqonly[bothexamplenum,]$treats,limmaonly[bothexamplenum,]$treats,psy[bothexamplenum,]$treats),
    v=c(deseqonly[bothexamplenum,]$var,limmaonly[bothexamplenum,]$var,psy[bothexamplenum,]$var),
    gene=c(deseqonly[bothexamplenum,]$identifier,limmaonly[bothexamplenum,]$identifier,psy[bothexamplenum,]$identifier))

lapply(1:nrow(conditions),function(n){
    nc=conditions[n,]
c=nc$c
i=nc$i
v=nc$v
gene=nc$gene
    cat("running",c,i,v,gene,"\n")
ressub <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
deseqres <- fread(paste0(deseqoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",deseqrun,".txt"))
ressub <- subset(ressub, cluster==c)
deseqressub <- subset(deseqres, cluster==c)
d_info <- deseqressub[deseqressub$identifier==gene,]
l_info <- ressub[ressub$identifier==gene,]

cluster_metadata_sce <- metadata_ls[[c]]
cluster_metadata <- data.frame(cluster_metadata_sce)
cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
if(job=="CZI"){
cluster_metadataL <- left_join(cluster_metadata,unique(eigenvec2[,c("Sample_ID","Lead","WHR")]),by="Sample_ID")
rownames(cluster_metadataL) <- cluster_metadataL$rowid
cluster_metadata <- cluster_metadataL
}
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
if(v=="factor_HS_CRP"){
cluster_metadata <- subset(cluster_metadata, HS_CRP<10) #advised to remove as likely an infection
cluster_metadata <- transform(cluster_metadata, factor_HS_CRP=as.factor(ifelse(HS_CRP<1, "low",ifelse(HS_CRP<3,"average","high"))))
cluster_metadata <- within(cluster_metadata, factor_HS_CRP <- relevel(factor_HS_CRP, ref = "average"))
}
cluster_metadata_t <- subset(cluster_metadata, treats==i)
if(job=="ALOFT"){
cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1",v)]
} else if (job=="CZI"){
cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age",v)]
}
if(v=="age" | v=="cage1"){
if(job=="ALOFT"){
cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex",v)]
} else if (job=="CZI"){
cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph",v)]
}}
if(iselcov){
cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","ISEL_Mean",v)]
}
if(WHRcov){
cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2","sex_alph","age","WHR",v)]
}
if(v=="CVDRISK"){
cluster_metadata_var <- cluster_metadata_t[,c("Sample_ID","PC1","PC2",v)]
}
cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]

dge <- DGEList(counts=cluster_counts_t)
dge <- calcNormFactors(dge)
varcol=cv_d[,v]
if(log2cpm){
cpm <- log2(cpm(dge)+1)
rownames(cpm) <- rownames(cluster_counts_t)
method="log2cpm"
}
if(voom){
if(job=="ALOFT"){
design <- model.matrix(~ 0 + as.numeric(varcol) + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
if(v=="cage1"){
design <- model.matrix(~ 0 + as.numeric(varcol) + as.factor(cv_d$Sex) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
}} 
if (job=="CZI" & !v=="CVDRISK"){
design <- model.matrix(~ 0 + as.numeric(varcol) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
if(v=="age"){
design <- model.matrix(~ 0 + as.numeric(varcol) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
}}
if(v=="CVDRISK"){
design <- model.matrix(~ 0 + as.numeric(varcol) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
}
if(iselcov){
design <- model.matrix(~ 0 + as.numeric(cv_d$ISEL_Mean) + as.numeric(varcol) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
}
if(WHRcov){
design <- model.matrix(~ 0 + as.numeric(cv_d$WHR) + as.numeric(varcol) + as.factor(cv_d$sex_alph) + as.numeric(cv_d$age) + as.numeric(cv_d$PC1) + as.numeric(cv_d$PC2))
}
v_e <- voom(dge, design, plot=FALSE, normalize.method="quantile")
genes_normed_baseline <- data.frame(v_e$E)
gene_method <- genes_normed_baseline[rownames(genes_normed_baseline)==gene,]
method="voom"
}
if(cpmqqnorm){
cpm <- cpm(dge)
rownames(cpm) <- rownames(cluster_counts_t)
mat_qnorm <- t(apply(cpm,1,function(x){qqnorm(rank(x, ties.method = "random"), plot = F)$x}))
colnames(mat_qnorm) <- colnames(cpm)
gene_method <- mat_qnorm[rownames(mat_qnorm)==gene,]
method="cpmqqnorm"
}
gene_ind <- gsub("[.]","-", names(gene_method))
cv_d2 <- cv_d[which(rownames(cv_d) %in% gene_ind),]

#plot x=variable, y=cpm
plotdf <- data.frame(x=as.vector(cv_d2[,v]),y=as.vector(unlist(gene_method)))

p <- ggplot(plotdf, aes(x=x, y=y)) +
  theme_bw()+
  geom_point()+ #aes(color=sig)
    #geom_abline(slope=1,intercept=0) +
  geom_smooth(method = "lm", fill = NA)+
  ylab(paste0(gene," ",method)) +
  xlab(paste0(v)) +
  ggtitle(paste0(c," ",i," ",v,"\n","DESeq logFC=",round(d_info$logFC,2)," DESeq padj=",scientific(d_info$padj, digits = 3)," DESeq pvalue=",scientific(d_info$pvalue, digits = 3),
    "\n","limma logFC=",round(l_info$logFC,2)," limma padj=",scientific(l_info$adj.P.Val, digits = 3)," limma pvalue=",scientific(l_info$P.Value, digits = 3)))+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),plot.title = element_text(size = rel(1.1), face = "bold"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 6, height = 6, file=paste0(outFolder,"figures/",c,".",i,".",v,".",gene,".",run,".png"), pointsize=12, bg = "transparent", canvas = "white", units = "in", res = 500)
print(p)
dev.off()

})

#plotting number of degs
#voom vs deseq
distinct_colors <-c("#db6e00","#4329cb","#5fab00","#c331df","#019b32","#b95aff","#c2a700","#0036b7","#f9bc28","#6d0096",
"#79db87","#ff7dfa","#bcd05c","#014ab7","#9a9200","#688bff","#ff3123","#01b281","#a80044","#01aea8","#c14600","#57aeff",
"#6d0004","#0086ce","#ff9b6c","#0155a6","#6b6600","#1f2865","#d7b282","#81004c","#006a34","#d5bbf9","#273a00","#ffabc9",
"#724500","#714771","#ff928d","#55133b")

varsrun = subset(variables_df, variable %in% c(allvars))
statsfull <- unique(fread(paste0(outFolder,"summary.",run,".txt")))

for (i in treatmentsfirst){
  subvars <- ldply(lapply(varsrun$variable, function(v){
    if(isTRUE(file.size(paste0(deseqoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",v,"-",i,".",deseqrun,".txt")) > 0)){
    deseqstat <- fread(paste0(deseqoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",v,"-",i,".",deseqrun,".txt"))
    limmastat <- subset(statsfull, treat==i & variable==v)
    m <- merge(limmastat,deseqstat,by=c("cluster","variable"))
    names(m)[c(7,14)] <- c("limma_DEGs","deseq_DEGs")
    m_s <- subset(m, limma_DEGs>49 | deseq_DEGs>49)
    return(m_s)
    }
    }),data.frame)            
  #subvars <- subvars[sapply(subvars, nrow)>0]

p <- ggplot(subvars, aes(x=limma_DEGs, y=deseq_DEGs)) +
  theme_bw()+
  geom_point(size=7,aes(shape=cluster,color=variable,fill=variable))+ #aes(color=sig)
  scale_shape_manual(values = c(7,8,11,15:18,25))+
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
png(width = 10, height = 10, file=paste0(outFolder,"figures/",i,".",run,".limmavsdeseq_numDEGs.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
}

#why so many DEGs in C6, plot logFC in C6 vs C0,C1,C2,C9,C10
tcellfc <- ldply(lapply(treatmentsfirst,function(i){
  subvars <- ldply(lapply(varsrun$variable, function(v){
    cat("running",i,v,"\n")
    ressub <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
    nDEGsC6=length(unique(subset(ressub,cluster=="C6" & adj.P.Val<0.1)$identifier))
    if(nDEGsC6>49){
        dc <- dcast(ressub, identifier+treat+variable ~ cluster, value.var="logFC")
    return(dc)
    }
    }),data.frame)
  }), data.frame)


#cols <- c("4both" = "orange", "1Not_Sig" = "grey","3limma_sig" = "red", "2deseq_sig" = "blue")
for (i in treatmentsfirst){
for(c in c("C0","C1","C2","C9","C10")){
    cat("running",i,c,"\n")
    tcellfct <- subset(tcellfc, treat==i)
p <- ggplot(tcellfct, aes_string(x = "C6", y = c)) +
  theme_bw()+
    #facet_grid(cluster~var,scales="free")+
  geom_point()+ #aes(color=sig)
    #scale_colour_manual(values = cols)+
    geom_abline(slope=1,intercept=0,linetype="dotted") +
      geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
  ylab(paste0(c," logFC")) +
  xlab("C6 logFC") +
  ggtitle(paste0("C6vs",c,"_logfc in ",i))+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),plot.title = element_text(size = rel(1.3), face = "bold"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 8, height = 8, file=paste0(outFolder,"figures/",i,".C6vs",c,"_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(p)
dev.off()
}
}

#investigate single variable in control
v="cssdh"
i="CTRL"
for(c in c("C0","C1","C2","C9","C10")){
    cat("running",i,c,"\n")
    tcellfct <- subset(tcellfc, treat==i & variable==v)
p <- ggplot(tcellfct, aes_string(x = "C6", y = c)) +
  theme_bw()+
    #facet_grid(cluster~var,scales="free")+
  geom_point()+ #aes(color=sig)
    #scale_colour_manual(values = cols)+
    geom_abline(slope=1,intercept=0,linetype="dotted") +
      geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
  ylab(paste0(c," logFC")) +
  xlab("C6 logFC") +
  ggtitle(paste0("C6vs",c,"_logfc in ",v," ",i))+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),plot.title = element_text(size = rel(1.3), face = "bold"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 8, height = 8, file=paste0(outFolder,"figures/",i,".",v,".C6vs",c,"_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(p)
dev.off()
}
###################################################################################################################33
######################################################################################################################









#################
#compare limma var vs limma int marginal vs deseq

threshold=0.10
res <- fread(paste0(outFolder,"all_limmares.",run,".txt"))

limmarun="treatment_withCOMBAT_limma"
limmaoutFolder=paste0(baseoutFolder,limmarun,"/")
limmaintres <- fread(paste0(limmaoutFolder,"all_limmares_treatinteraction.txt"))
limmaintres_m <- subset(limmaintres, term=="variable")

statsfull <- fread(paste0(outFolder,"summary.",run,".txt"))
statsfull <- transform(statsfull, analysis="limma_var")
statsfull_m <- statsfull[,c(1:3,7:8)]


limmaintstatsfull <- fread(paste0(limmaoutFolder,"summary.treatinteraction.txt"))
limmaintstatsfull <- transform(limmaintstatsfull, treat=sapply(strsplit(contrast,"_"),function(y) y[1]),analysis="limmaint_varmarginal")
limmaintstatsfull_m <- limmaintstatsfull[,c("cluster","treat","variable","varDEGs","analysis")]


deseqrun="income_PCs_sex_age_and_treats_adjusted_withWave"
deseqoutFolder=paste0(baseoutFolder,deseqrun,"/")

deseq_degs <- ldply(lapply(c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX"), function(i){
  psy <- ldply(lapply(c(varsrun$variable), function(var){
#psy <- ldply(lapply(variables_df$variable, function(var){
    cat("running",i,var,"\n")
    if(file.exists(paste0(deseqoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",deseqrun,".txt"))){
    stats <- fread(paste0(deseqoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",deseqrun,".txt"))
    stats=transform(stats, treat=i, analysis="DESeq")
    return(stats)
    }
    }),data.frame)
  return(psy)
  }), data.frame)
deseq_degs_m <- deseq_degs[,c("cluster","treat","variable","DEGs_FDR_10","analysis")]
names(deseq_degs_m)[4] <- "varDEGs"

merged_all <- rbind(statsfull_m,limmaintstatsfull_m,deseq_degs_m)
merged_all_w <- reshape2::dcast(merged_all, cluster + treat + variable ~ analysis, value.var="varDEGs")

lapply(c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX"), function(i){
merged_sub <- subset(merged_all_w, treat==i)
df50 <- merged_sub %>% filter_at(vars(-c(1:3)), any_vars(. > 49))

dft <- df50 %>% flextable() %>% set_caption(caption = run) 
#dft <- df50 %>% flextable() %>% split_header(sep="_") %>% set_caption(caption = run) 
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1, j = c(1:length(df50)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(df50)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(seq(3,length(df50),by=3)), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/comparing_limmaint_limmavar_deseq.",i,".table.png"))
})




#plotting logFC between deseq and limma
cols <- c("4both" = "orange", "1Not_Sig" = "grey","3limma_sig" = "red", "2deseq_sig" = "blue")
if(job=="ALOFT"){
    focusvars=c("pnsi","cdres","criskf") #mono_av
    outliervar="mono_av"
} else if (job=="CZI"){
    focusvars=firstrunvars #can change this to reflect whats seen in above plots
}

lapply(treatments,function(i){
    cat("running",i,"\n")
  psyall <- ldply(lapply(firstrunvars, function(v){
    deseqres <- fread(paste0(deseqoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",i,".",deseqrun,".txt"))
    ressub <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
    names(ressub)[c(2,3)] <- c("treats","var")
    m <- merge(ressub,deseqres,by=c("identifier","cluster","treats","var"))
    return(m)
    }),data.frame)
psyall <- psyall[order(psyall$padj,psyall$adj.P.Val),]

    psyall <- transform(psyall, sig=as.factor(ifelse(adj.P.Val<threshold & padj<threshold, "4both", ifelse(adj.P.Val<=threshold, "3limma_sig", ifelse(padj<=threshold, "2deseq_sig", "1Not_Sig")))))

psyall_psych <- subset(psyall, var %in% focusvars)
p <- ggplot(psyall_psych, aes(x=logFC.x, y=logFC.y)) +
  theme_bw()+
    facet_grid(cluster~var,scales="free")+
  geom_point(aes(color=sig))+ #
    scale_colour_manual(values = cols)+
    geom_abline(slope=1,intercept=0,linetype="dotted") +
      geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
  ylab("DESeq logFC") +
  xlab("limma logFC") +
  ggtitle(paste0(i))+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),plot.title = element_text(size = rel(1.3), face = "bold"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 12, height = 14, file=paste0(outFolder,"figures/",i,"deseq_limma_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(p)
dev.off()
#just mono_av
mono <- subset(psyall, var==outliervar)
p <- ggplot(mono, aes(x=logFC.x, y=logFC.y)) +
  theme_bw()+
    facet_grid(.~cluster,scales="free")+
  geom_point(aes(color=sig))+ #
    scale_colour_manual(values = cols)+
    geom_abline(slope=1,intercept=0,linetype="dotted") +
      geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
  ylab("DESeq logFC") +
  xlab("limma logFC") +
  ggtitle(paste0(i))+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),plot.title = element_text(size = rel(1.3), face = "bold"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 16, height = 4, file=paste0(outFolder,"figures/mono_av.",i,"deseq_limma_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(p)
dev.off()
})
########
#PFAS qqplots
#PFAS plotting only
#run="limma_var"
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
#outFolder=paste0(baseoutFolder,run,"/")

library(ggrastr)
plotDF <- ldply(lapply(treatmentsfirst,function(i){
    resall <- ldply(lapply("CVDRISK",function(v){
    res <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
    ntest <- nrow(res)
       res <- res%>%
           arrange(P.Value)%>%
           mutate(observed=-log10(P.Value), expected=-log10(ppoints(ntest)))
           return(res)
        }),data.frame)
    return(resall)
    }), data.frame)

var2 <- sort(unique(plotDF$variable))
varSel2 <- c("sumPFAS","log_sumPFAS")
#varSel2 <- allPFAS[grep("Log",allPFAS)]
#varSel2 <- "CVDRISK"
plotDF2 <- plotDF%>%filter(variable%in%varSel2)
zoom="CVDRISK"

p0 <- ggplot(plotDF2, aes(x=expected, y=observed, color=cluster))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    facet_grid(.~treat, scales="free_y")+
    #facet_grid(variable~treat, scales="free_y")+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          axis.title=element_text(size=12),
          axis.text=element_text(size=10),
          strip.text=element_text(size=12))

figfn <- paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".",zoom,".qqplot.png")
ggsave(figfn, p0, width=1000, height=800, units="px", dpi=120)

#PFAS only
library(flextable)
library(ftExtra)
library(rlist)
library(officer)
#[!psychvarstorun %in% puberty]
 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
zoom="Lead"
statsfull <- unique(fread(paste0(outFolder,"summary.",run,".txt")))

for (i in c(treatmentsfirst)){
#subvars <- ldply(lapply(c(variables_df$variable), function(var){
subvars <- ldply(lapply(c("Lead"), function(v){
    cat("running",i,v,"\n")
    stats <- subset(statsfull, treat==i & variable==v)
        limmares <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
        df <- ddply(limmares, "cluster", plyr::summarize, ngenes=length(unique(identifier)))
        newls <- lapply(names(counts_ls),function(c){
            stats_c <- subset(stats,cluster==c)
            dc <- subset(df, cluster==c)
            d <- data.frame(numind=stats_c$number_individuals, ngenes=dc$ngenes,sigDEGs=stats_c$varDEGs)
            colnames(d) <- paste(c,colnames(d),sep=".")
            return(d)
            })
            newls <- newls[sapply(newls, nrow)>0]
            newtable <- data.frame(variable=v, description=variables_df[variables_df$variable==v,]$description,list.cbind(newls))
            newtable$description <- ifelse(is.na(newtable$description) | newtable$description=="", v, newtable$description)
            return(newtable)
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

zoom="Lead"
for (i in c(treatmentsfirst)){
subvars <- ldply(lapply(c("Lead"), function(v){
    cat("running",i,v,"\n")
    stats <- subset(statsfull, treat==i & variable==v)
    if(dim(subset(stats, varDEGs>50))[1]<1){ #first subset variables that have at least one cluster with 50DEGs
        cat("not enough DEGs","\n")
        rm(stats)
        } else{
        limmares <- fread(paste0(outFolder,"limmares.",i,".",v,".",run,".txt"))
        limmares <- transform(limmares, regulation=ifelse(logFC>0, "up", "down"))
        dessig <- subset(limmares, adj.P.Val<0.1)
        return(dessig)
    }
}),data.frame)

if(dim(subvars)[1]<1){
    cat("no DEGs>50 \n")
    } else {
# Count the number of genes in each category
data_summary <- ddply(subvars, c("cluster","variable","regulation"),plyr::summarize,
    nDEG=length(regulation))
data_summary <- transform(data_summary, nDEG=ifelse(regulation=="up",nDEG,-(nDEG)))
group_colors <- c("up" = "orange","down" = "purple")

p <- ggplot(data_summary, aes(x = cluster, y = nDEG, fill = regulation)) +
  geom_col() + #position = "dodge"
  coord_flip() +
  geom_hline(yintercept=0)+
  #facet_wrap(.~variable)+
  scale_fill_manual(values = group_colors) +
  labs(title = paste0("Up- and Down-Regulated Genes in ",i), x = "Cluster", y = "Number of DEGs", fill = "Regulation") +
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 8, height = 8, file=paste0(outFolder,"figures/",i,".",zoom,".",run,".DEGupvdownbar.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
}
}