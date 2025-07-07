require(broom) 
library(plyr)   
library(tidyverse)
library(edgeR)
library(limma)
library(annotables)
library(data.table)
library(ggpubr)

#load in original counts info for clusters and treats used
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.RData")
load(opfn)
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/treatment_withCOMBAT_lm/"

# load annotation
u_eigenvec2 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt")
cv <- unique(u_eigenvec2)
##check for individuals with more than one sample
cv %>% count(Sample_ID) %>% dplyr::filter(n>1)
eigenvec2 <- u_eigenvec2 %>% dplyr::filter(SCAIP7_18 == 1 | SCAIP1_6 == 1)
eigenvec2$Sample_ID <- eigenvec2$dbgap.ID
psychvarstorun <- colnames(eigenvec2[,24:59])

combatrun="income_PCs_sex_age_and_treats_adjusted"
combatfolder=paste0(baseoutFolder,combatrun,"/")
contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA-DEX"))

lmcoefs <- ldply(lapply(names(counts_ls),function(cluster){
    #cluster="C6"
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    if(job=="ALOFT"){
        adjusted_counts <- adjusted
    }
    vars <- ldply(lapply(allvars[c(1:20)],function(var){
        #var="pedu"
        cat("running ", cluster, var, "\n")
    df <- ldply(lapply(1:length(contrastdf$control),function(x){
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_metadata_var <- cluster_metadata_var[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        #cdf <- subset(plyr::count(cluster_metadata_var$Sample_ID),freq<2)
        #cluster_metadata_var <- subset(cluster_metadata_var, !Sample_ID %in% cdf$x)
        cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_var))]
        cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        cv_d <- transform(cv_d, treats=factor(treats))
        ## Normalization of data
        # make edgeR object
        dge <- DGEList(counts=cluster_counts_t)
        #Transform counts to counts per million
        dge <- calcNormFactors(dge)
        cpm <- log2(cpm(dge)+1)
        rownames(cpm) <- rownames(cluster_counts_t)
        samples <- dim(cluster_counts_t)[2]
        varcol=cv_d[,var]
        design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
        tcpm <- t(cpm)
        model <- lm(t(cpm) ~ design, data = data.frame(tcpm, cv_d))
        tmodel <- as.data.frame(tidy(model))
        intmodel <- tmodel[grepl("varcol", tmodel$term) & grepl("treats", tmodel$term),]
        intmodel <- transform(intmodel, adj.P.Val=p.adjust(p.value,method="BH"),identifier=response, cluster=cluster, contrast=contrast, variable=var, term="interaction")
        intdf <- intmodel[,c("cluster","contrast","variable","term","identifier","estimate","std.error","statistic","p.value","adj.P.Val")]

        treatmodel <- tmodel[!grepl("varcol", tmodel$term) & grepl("treats", tmodel$term),]
        treatmodel <- transform(treatmodel, adj.P.Val=p.adjust(p.value,method="BH"),identifier=response, cluster=cluster, contrast=contrast, variable=var, term="treat")
        treatdf <- treatmodel[,c("cluster","contrast","variable","term","identifier","estimate","std.error","statistic","p.value","adj.P.Val")]

        varmodel <- tmodel[grepl("varcol", tmodel$term) & !grepl("treats", tmodel$term),]
        varmodel <- transform(varmodel, adj.P.Val=p.adjust(p.value,method="BH"),identifier=response, cluster=cluster, contrast=contrast, variable=var, term="variable")
        vardf <- varmodel[,c("cluster","contrast","variable","term","identifier","estimate","std.error","statistic","p.value","adj.P.Val")]

        fwrite(intdf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"lmres/",cluster,".",contrast,".",var,".lm_interaction.treatvarint.txt"))
        fwrite(treatdf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"lmres/",cluster,".",contrast,".",var,".lm_treat.treatvarint.txt"))
        fwrite(vardf, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"lmres/",cluster,".",contrast,".",var,".lm_variable.treatvarint.txt"))
        table <- data.frame(cluster=cluster, contrast=contrast, variable=var)
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
        table$gene_number <- paste(nrow(cluster_counts_t))
        table$intDEGs=nrow(subset(intdf, adj.P.Val<0.1))
        table$treatDEGs=nrow(subset(treatdf, adj.P.Val<0.1))
        table$varDEGs=nrow(subset(vardf, adj.P.Val<0.1))
        return(table)
        }), data.frame)
    return(df)
    }), data.frame)
return(vars)
}), data.frame)
fwrite(lmcoefs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"summary.treatinteraction.txt"))
     
myDir <- paste0(outFolder,"lmres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grepl("lm_.*.treatvarint.txt",filenames)]
data_names <- gsub(".treatvarint.txt", "", filenames)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
res <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(res, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"all_lmres_treatinteraction.txt"))


#want to plot limma vs lm in DESeq model
limmarun="treatment_withCOMBAT_limma"
limmaoutFolder=paste0(baseoutFolder,limmarun,"/")
threshold=0.10

    lmres <- fread(paste0(outFolder,"all_lmres_treatinteraction.txt"))
    limmares <- fread(paste0(limmaoutFolder,"all_limmares_treatinteraction.txt"))
    merged_treat <- merge(lmres,limmares,by=c("identifier","cluster","variable","contrast","term"))
    subvars <- transform(merged_treat, sig=as.factor(ifelse(adj.P.Val.x<threshold & adj.P.Val.y<threshold, "4both", ifelse(adj.P.Val.x<=threshold, "3lm_sig", ifelse(adj.P.Val.y<=threshold, "2limma_sig", "1Not_Sig")))))
    subvarssig <- subset(subvars, !sig=="1Not_Sig")
lapply(split(subvars,subvars$contrast),function(c){
    lapply(split(c,c$variable),function(v){
cat("running",unique(v$contrast),":",unique(v$variable))
p <- ggplot(v, aes(x=estimate, y=logFC)) +
  theme_bw()+
  ggtitle(paste0(unique(v$contrast),":", unique(v$variable)))+
  facet_grid(cluster~term)+
  xlab("lm_estimate")+
  ylab("limma_logFC")+
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
      png(width = 12, height = 12, file=paste0(outFolder,"figures/limmavslm_",unique(v$contrast),":",unique(v$variable),".png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
})
})


# Stacked

library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(data.table)
library(reshape)
library(stringr)
library(data.table)
library(ggrepel)

lapply(split(subvarssig,subvarssig$contrast),function(c){
    lapply(split(c,c$variable),function(v){
cat("running",unique(v$contrast),":",unique(v$variable))

count_df <- v %>%
  group_by(cluster, term, sig) %>%
  summarise(count = n(), .groups = "drop")

fig0 <- ggplot(count_df, aes(x = count, y = cluster, fill = sig)) +
  geom_bar(stat = "identity") +
  #scale_fill_manual(values = group_colors) +
  theme_bw() +
    facet_grid(.~term)+
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("DEGs Overlap Between limma and linear model", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
png(width = 12, height = 7, file=paste0(outFolder,"figures/limmavslm_",unique(v$contrast),":",unique(v$variable),".stackedbar.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(fig0)
dev.off()
    })
})


    subvars <- transform(merged_treat, sig=as.factor(ifelse(adj.P.Val.x<threshold & adj.P.Val.y<threshold, "both", ifelse(adj.P.Val.x<=threshold, "lm_sig", ifelse(adj.P.Val.y<=threshold, "limma_sig", "Not_Sig")))))
    subvarssig <- subset(subvars, !sig=="Not_Sig")
c=split(subvars,subvars$contrast)[[1]]
v=split(c,c$variable)[["IL13_co"]]

listc <- ldply(lapply(split(subvarssig,subvarssig$contrast),function(c){
    list_v <- ldply(lapply(split(c,c$variable),function(v){
        cat("running",unique(v$contrast),":",unique(v$variable))
count_dfd <- ddply(v, c("cluster","term"),plyr::summarize,
    nDEGslimma=sum(adj.P.Val.y<=threshold,na.rm=T),
    nDEGslm=sum(adj.P.Val.x<=threshold,na.rm=T)
)
count_dfdm <- melt(count_dfd)
count_dfw <- reshape2::dcast(count_dfdm,cluster ~ term+variable)
df <- data.frame(contrast=unique(v$contrast),var=unique(v$variable),count_dfw)
df50 <- df %>% filter_at(vars(-c(1:3)), any_vars(. > 49))

dft <- df %>% flextable() %>% split_header(sep="_") %>% set_caption(caption = run) 
#dft <- df50 %>% flextable() %>% split_header(sep="_") %>% set_caption(caption = run) 
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1, j = c(1:length(df50)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(df50)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(seq(3,length(df50),by=2)), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/limmavslm_",unique(df50$contrast),":",unique(df50$var),".table.png"))
return(df50)
}), data.frame)
return(list_v)
}),data.frame)

listc <- subset(listc, interaction_nDEGslimma >49 | interaction_nDEGslm>49)
listc1 <- listc[,-1]
dft <- listc1 %>% flextable() %>% split_header(sep="_") %>% set_caption(caption = run) 
#dft <- df50 %>% flextable() %>% split_header(sep="_") %>% set_caption(caption = run) 
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1, j = c(1:length(listc1)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(2:length(listc1)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(seq(3,length(listc1),by=2)), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)
flextable::save_as_image(
  dft,
  path = paste0(outFolder,"figures/limmavslm_intdeg50.table.png"))







lmcoefs <- ldply(lapply(names(counts_ls),function(cluster){
    #cluster="C0"
        cat("running ", cluster, "\n")
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    df <- ldply(lapply(1:length(contrastdf$control),function(x){
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        ## Normalization of data
        # make edgeR object
        dge <- DGEList(counts=cluster_counts_t)
        #Transform counts to counts per million
        dge <- calcNormFactors(dge)
        cpm <- log2(cpm(dge)+1)
        rownames(cpm) <- rownames(cluster_counts_t)
        samples <- dim(cluster_counts_t)[2]
        design <- model.matrix(~ cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
        model <- lm(t(cpm) ~ design, data = data.frame(t(cpm), cv_d))
        tmodel <- as.data.frame(tidy(model))
        tmodel <- transform(tmodel, cluster=cluster, term=ifelse(grepl("treats", term), contrast, term))
        tmodeltreat <- subset(tmodel, term==contrast)
        return(tmodeltreat)
    }), data.frame)
    return(df)
}), data.frame)

lmcoefspha <- subset(lmcoefs, term=="PHA_vs_CTRL")
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

var="pnsi"
lmcoefs <- ldply(lapply(names(counts_ls),function(cluster){
    #cluster="C0"
        cat("running ", cluster, "\n")
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    df <- ldply(lapply(1:length(contrastdf$control),function(x){
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_metadata_var <- cluster_metadata_var[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        ## Normalization of data
        # make edgeR object
        dge <- DGEList(counts=cluster_counts_t)
        #Transform counts to counts per million
        dge <- calcNormFactors(dge)
        cpm <- log2(cpm(dge)+1)
        rownames(cpm) <- rownames(cluster_counts_t)
        samples <- dim(cluster_counts_t)[2]
        design <- model.matrix(~ as.numeric(cv_d$pnsi) + cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
        model <- lm(t(cpm) ~ design, data = data.frame(t(cpm), cv_d))
        tmodel <- as.data.frame(tidy(model))
        tmodel <- transform(tmodel, cluster=cluster, term=ifelse(grepl("treats", term), contrast, term))
        tmodel <- transform(tmodel, cluster=cluster, term=ifelse(grepl(var, term), var, term))
        tmodeltreat <- subset(tmodel, term==var)
        return(tmodeltreat)
    }), data.frame)
    return(df)
}), data.frame)

#lmcoefspha <- subset(lmcoefs, term=="PHA_vs_CTRL")
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

var="pnsi"
puberty <-c("cgpd","cgpd5","pdpds5","cbpd","ppdpds","pspds")
allvars <- psychvarstorun[c(1:length(psychvarstorun))][!psychvarstorun[c(1:length(psychvarstorun))] %in% puberty]
firstrunvars <- allvars[c(1:10)]
secondrunvars <- c(allvars[c(1:20)],"cage1")

lmcoefs <- ldply(lapply(names(counts_ls),function(cluster){
    #cluster="C0"
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    vars <- ldply(lapply(allvars[c(1:20)],function(var){
        cat("running ", cluster, var, "\n")
    df <- ldply(lapply(1:length(contrastdf$control),function(x){
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_metadata_var <- cluster_metadata_var[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        #cdf <- subset(plyr::count(cluster_metadata_var$Sample_ID),freq<2)
        #cluster_metadata_var <- subset(cluster_metadata_var, !Sample_ID %in% cdf$x)
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        cv_d <- transform(cv_d, treats=factor(treats))
        ## Normalization of data
        # make edgeR object
        dge <- DGEList(counts=cluster_counts_t)
        #Transform counts to counts per million
        dge <- calcNormFactors(dge)
        cpm <- log2(cpm(dge)+1)
        rownames(cpm) <- rownames(cluster_counts_t)
        samples <- dim(cluster_counts_t)[2]
        varcol=cv_d[,var]
        design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
        model <- lm(t(cpm) ~ design, data = data.frame(t(cpm), cv_d))
        tmodel <- as.data.frame(tidy(model))
        tmodel <- transform(tmodel, cluster=cluster, contrast=contrast, variable=var, term=ifelse(grepl("treats", term) & grepl(var, term), paste0(contrast,":",var), ifelse(grepl("treats", term),contrast, ifelse(grepl(var, term), var, term))))
        tmodeltreat <- subset(tmodel, term==paste0(contrast,":",var) | term==contrast | term==var)
        tmodeltreata <- ldply(lapply(unique(tmodeltreat$term),function(t){
            dft <- subset(tmodeltreat, term==t)
            names(dft)[1] <- "identifier"
            dft <- transform(dft, padj=p.adjust(p.value,method="BH"))
            fwrite(dft, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".",var,".lm_",t,".treatvarint.txt"))
            table <- data.frame(variable= var, cluster=cluster,term=t)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts_t))
    table$tested_genes <- paste(nrow(dft))
    table$DEGs_FDR <- paste(length(which(dft$padj<0.05)))
    table$DEGs_FDR_10 <- paste(length(which(dft$padj<0.1)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".",var,".lmstats_",t,".treatvarint.txt"))
        return(dft)
        }),data.frame)
        return(tmodeltreata)
    }), data.frame)
    return(df)
    }), data.frame)
    return(vars)
}), data.frame)
fwrite(lmcoefs, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"all.treatvarint.txt"))



lmcoefspha <- subset(lmcoefs, term=="PHA_vs_CTRL" & variable==var)
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

lmcoefspha <- subset(lmcoefs, term=="PHA_vs_CTRL:pnsi")
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)

lmcoefspha <- subset(lmcoefs, term=="pnsi" & contrast=="PHA_vs_CTRL")
lmcoefsphaccl20 <- subset(lmcoefspha, response=="CCL20" & cluster=="C6")
lmcoefsphaNFKB1 <- subset(lmcoefspha, response=="NFKB1" & cluster=="C6")
lmcoefsphaCAMK1D <- subset(lmcoefspha, response=="CAMK1D" & cluster=="C6")
lmcoefsphaFOS <- subset(lmcoefspha, response=="FOS" & cluster=="C6")
lmcoefsphaLPYPLA1 <- subset(lmcoefspha, response=="LYPLA1" & cluster=="C6")
lmcoefsphaLST3GAL3 <- subset(lmcoefspha, response=="ST3GAL3" & cluster=="C6")
alltested <- rbind(lmcoefsphaccl20,lmcoefsphaNFKB1,lmcoefsphaCAMK1D,lmcoefsphaFOS,lmcoefsphaLPYPLA1,lmcoefsphaLST3GAL3)


    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".",var,".lmstats_",t,".treatvarint.txt"))

myDir <- paste0(outFolder)
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grepl(".lmstats_",filenames)]
data_names <- gsub(".treatvarint.txt", "", filenames)
data_names <- gsub(".lmstats_", "", filenames)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
res <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(res, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"all_stats_lm_treatinteraction.txt"))



    ldply(lapply(names(counts_ls),function(cl){
        ldply(lapply(1:length(contrastdf$control),function(x){
            con=contrastdf[x,]
            co=paste0(con$treatment,"_vs_",con$control)
            df <- subset(lmcoefs, cluster==cl & contrast==co)
            var=unique(df$variable)

    })