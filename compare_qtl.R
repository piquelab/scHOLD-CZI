R
library(data.table)
library(ggplot2)
library(plyr);library(dplyr)
library(tidyverse)

cluster="C4"
treat="CTRL"
baseoutFolder=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,"ALL.0.2.50.DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)

jaxoutFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/results/"
tensoroutFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/"

filenames <- list.files(jaxoutFolder) #file list from directory
jaxfilenames1 <- filenames[grep(".GEPCs.txt", filenames)] #pick specific files from list

filenames <- list.files(tensoroutFolder) #file list from directory
tensorfilenames1 <- filenames[grep(".GEPCs.txt", filenames)] #pick specific files from list

vlist <- ldply(lapply(names(counts_ls),function(c){
jax_pc_signif_pairs <- fread(paste0(jaxoutFolder,jaxfilenames1[grep(paste0(c,".",treat,".best"),jaxfilenames1)]))
jax_pc_signif_pairs <- transform(jax_pc_signif_pairs,method="jaxQTL")
tensor_pc_signif_pairs <- fread(paste0(tensoroutFolder,tensorfilenames1[grep(paste0(c,".",treat,".best"),tensorfilenames1)]))
tensor_pc_signif_pairs <- transform(tensor_pc_signif_pairs,method="tensorQTL")

all_PCs_r_best <- merge(jax_pc_signif_pairs[,c("phenotype_id","variant_id","method","pval_beta","padj")],tensor_pc_signif_pairs[,c("phenotype_id","variant_id","method","pval_beta","qval")],by=c("phenotype_id","variant_id"))
all_PCs_r_best <- transform(all_PCs_r_best, cluster=c)
all_PCs_r_best$pval_beta.x[all_PCs_r_best$pval_beta.x<1e-20] <- 1e-20
all_PCs_r_best$pval_beta.y[all_PCs_r_best$pval_beta.y<1e-20] <- 1e-20

all_PCs_r_best_m <- melt(all_PCs_r_best)

    v <- transform(all_PCs_r_best_m, cluster=as.factor(cluster),variable=as.factor(variable))
    v <- v[!v$variable %in% c("qval","padj"),] %>%
   group_by(cluster,variable)%>%
   arrange(value) %>%
   mutate(observed=-log10(value), expected=-log10(ppoints(length(value))))
   v <- transform(v, variable=ifelse(variable=="pval_beta.x","jaxQTL","tensorQTL"))
   return(v)
}), data.frame)

p0 <- ggplot(vlist, aes(x=expected, y=observed, color=variable))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    facet_wrap(.~cluster, scales="free_y",ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    labs(color = "QTL method")+
    #ggtitle(paste0(unique(v$variable)))+
    theme_bw()
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/figures/jax_vs_tensorqtl_qqplot.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()


#bar plot
filenames <- list.files(jaxoutFolder) #file list from directory
jaxfilenames1 <- filenames[grep("_significant_topeeQTL_pairs.txt", filenames)] #pick specific files from list

filenames <- list.files(tensoroutFolder) #file list from directory
tensorfilenames1 <- filenames[grep("_significant_topeeQTL_pairs.txt", filenames)] #pick specific files from list

vlist <- ldply(lapply(names(counts_ls),function(c){
jax_pc_signif_pairs <- fread(paste0(jaxoutFolder,jaxfilenames1[grep(paste0(c,".",treat),jaxfilenames1)]))
jax_pc_signif_pairs <- transform(jax_pc_signif_pairs,method="jaxQTL",cluster=c)
tensor_pc_signif_pairs <- fread(paste0(tensoroutFolder,tensorfilenames1[grep(paste0(c,".",treat),tensorfilenames1)]))
tensor_pc_signif_pairs <- transform(tensor_pc_signif_pairs,method="tensorQTL",cluster=c)

df <- data.frame(cluster=c,jaxQTL_eGene=dim(jax_pc_signif_pairs)[1],tensorQTL_eGene=dim(tensor_pc_signif_pairs)[1])

   return(df)
}), data.frame)

justyna <- data.frame(cluster=c("C4","C5","C7"),justyna_fastqtl_eGene=c(364,291,253))
vlistm <- merge(vlist,justyna,by="cluster",all=T)
outtablem <- melt(vlistm)
p <- ggplot(outtablem, aes(fill=variable, y=value, x=cluster)) + 
    geom_bar(position="dodge", stat="identity")+
    labs(y="Number SNP:eGene pairs at FDR 10%",fill="QTL method")+
    theme_bw()
    #facet_wrap(.~description,ncol=4)
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/figures/jax_vs_tensorqtl_bar.png")
png(width = 9, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

#compare GxE
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/results/"
jaxint <- fread(paste0(resultsFolder,treat,".GxE_abundance_perclus.txt"))
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/"
tensorint <- fread(paste0(resultsFolder,treat,".GxE_abundance_perclus.txt"))
resultsFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/interaction/"
tensorint_builtin <- fread(paste0(resultsFolder,treat,".GxE_abundance_perclus.txt"))

all_PCs_r_best <- merge(jaxint[,c(1:5,9)],tensorint[,c(1:5,9)],by=c("phenotype_id","variant_id","cluster","treat","variable"))
#all_PCs_r_best <- merge(all_PCs_r_best1,tensorint_builtin[,c(2:6,21)],by=c("phenotype_id","variant_id","cluster","treat","variable")) #forgot it only outputs best not all

all_PCs_r_best$interaction_pval.x[all_PCs_r_best$interaction_pval.x<1e-20] <- 1e-20
all_PCs_r_best$interaction_pval.y[all_PCs_r_best$interaction_pval.y<1e-20] <- 1e-20
#all_PCs_r_best$pval_emt[all_PCs_r_best$pval_emt<1e-20] <- 1e-20

all_PCs_r_best_m <- melt(all_PCs_r_best)

    v <- transform(all_PCs_r_best_m, cluster=as.factor(cluster),QTL_method=as.factor(variable.1))
    v <- v %>%
   group_by(cluster,variable,QTL_method)%>%
   arrange(value) %>%
   mutate(observed=-log10(value), expected=-log10(ppoints(length(value))))
   v <- transform(v, QTL_method=ifelse(QTL_method=="interaction_pval.x","jaxQTL",ifelse(QTL_method=="interaction_pval.y","tensorQTL","tensorQTL_builtin")))

lapply(unique(v$variable),function(var) {
	vsub <- subset(v, variable==var)
    p0 <- ggplot(vsub, aes(x=expected, y=observed, color=QTL_method))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    facet_wrap(.~cluster, scales="free_y",ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    labs(color = "QTL method")+
    ggtitle(paste0(var))+
    theme_bw()
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/figures/jax_vs_tensorqtl_GxE_",var,"_qqplot.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()
})

jaxint <- jaxint %>%
group_by(cluster,variable)%>%
mutate(jaxQTL_eGene=sum(interaction_padj<0.1,na.rm=T)) %>% ungroup()

tensorint <- tensorint %>%
group_by(cluster,variable)%>%
mutate(tensorQTL_eGene=sum(interaction_padj<0.1,na.rm=T)) %>% ungroup()

tensorint_builtin <- tensorint_builtin %>%
group_by(cluster,variable)%>%
mutate(tensorQTL_builtin_eGene=sum(pval_adj_bh<0.1,na.rm=T)) %>% ungroup()

all_PCs_r_best1 <- merge(unique(jaxint[,c("cluster","variable","jaxQTL_eGene")]),unique(tensorint[,c("cluster","variable","tensorQTL_eGene")]))
all_PCs_r_best <- merge(all_PCs_r_best1,unique(tensorint_builtin[,c("cluster","variable","tensorQTL_builtin_eGene")]))

all_PCs_r_best_m <- melt(all_PCs_r_best)
names(all_PCs_r_best_m)[3] <- "QTL_method"
all_PCs_r_best_m <- transform(all_PCs_r_best_m, cluster=as.factor(cluster))
p <- ggplot(all_PCs_r_best_m, aes(fill=QTL_method, y=value, x=cluster)) + 
    geom_bar(position="dodge", stat="identity")+
    labs(y="Number SNP:eGene pairs at FDR 10%",fill="QTL method")+
    theme_bw()
    facet_wrap(.~variable,ncol=4)
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/jaxQTL/output/figures/jax_vs_tensorqtl_GxE_bar.png")
png(width = 9, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

########################################################################

#now comparing ctrl only tensorqtl
cluster_celltype_alltreat <- fread(paste0(base,method,"_pseudobulk_ctrl/",filter,"/",project,".",resset,".",dimset,".cluster_celltype.txt"))
cluster_celltype_alltreat <- transform(cluster_celltype_alltreat, cell_type=ifelse(cluster=="C1","Naive CD4+ T cells",
  ifelse(cluster=="C6","Naive CD4+ T cells",ifelse(cluster=="C7", "Monocytes",ifelse(cluster=="C9","Naive CD4+ T cells",ifelse(cluster=="C10","Memory CD4+ T cells",cell_type))))))
cluster_celltype_alltreat <- rbind(cluster_celltype_alltreat,data.frame(cluster=c("C12"),cell_type=c("B cell")))

cluster_celltype <- fread(paste0(base,method,"_pseudobulk_ctrl/",filter,"/",project,".",resset,".",dimset,".cluster_celltype.txt"))

colnames(cluster_celltype) <- c("cluster.ctrlonly","celltype.ctrlonly")
colnames(cluster_celltype_alltreat) <- c("cluster.alltreat","celltype.alltreat")

#cluster_celltype_alltreat$seuratcluster <- gsub("C","",cluster_celltype_alltreat$cluster.alltreat)

#0/0,1/1,2/3,3/4,4/5,5/7,6/8,7/9,skip8/,/skip2,9/10,10/12,/skip11
celltype <- data.frame(cluster.alltreat=c("C0","C1","C3","C4","C5","C7","C8","C9","C10","C12"),
  cluster.ctrlonly=c("C0","C1","C2","C3","C4","C5","C6","C7","C9","C10"))

celltypeorder <- cbind(cluster_celltype_alltreat[match(celltype$cluster.alltreat,cluster_celltype_alltreat$cluster.alltreat),],cluster_celltype[match(celltype$cluster.ctrlonly,cluster_celltype$cluster.ctrlonly),])

celltypeorder <- fread(paste0(base,"celltype_alltreat_controlonly.txt"))

treat="CTRL"
dexbaseoutFolder=paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/")

ctrltensoroutFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/results/"
tensoroutFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/results/"

filenames <- list.files(ctrltensoroutFolder) #file list from directory
ctrltensorfilenames1 <- filenames[grep(".GEPCs.txt", filenames)] #pick specific files from list

filenames <- list.files(tensoroutFolder) #file list from directory
tensorfilenames1 <- filenames[grep(".GEPCs.txt", filenames)] #pick specific files from list

vlist <- ldply(lapply(celltypeorder$cluster.ctrlonly[!celltypeorder$cluster.ctrlonly %in% c("C9","C10")],function(c){
    cat("running",c,"\n")
    mc <- celltypeorder[celltypeorder$cluster.ctrlonly==c,]
ctrltensor_pc_signif_pairs <- fread(paste0(ctrltensoroutFolder,ctrltensorfilenames1[grep(paste0(c,".",treat,".best"),ctrltensorfilenames1)]))
ctrltensor_pc_signif_pairs <- transform(ctrltensor_pc_signif_pairs,method="ctrl_tensorQTL")
tensor_pc_signif_pairs <- fread(paste0(tensoroutFolder,tensorfilenames1[grep(paste0(mc$cluster.alltreat,".",treat,".best"),tensorfilenames1)]))
tensor_pc_signif_pairs <- transform(tensor_pc_signif_pairs,method="alltreat_tensorQTL")

all_PCs_r_best <- merge(ctrltensor_pc_signif_pairs[,c("phenotype_id","variant_id","method","pval_beta","qval")],tensor_pc_signif_pairs[,c("phenotype_id","variant_id","method","pval_beta","qval")],by=c("phenotype_id","variant_id"))
all_PCs_r_best <- transform(all_PCs_r_best, cluster=c)
all_PCs_r_best$pval_beta.x[all_PCs_r_best$pval_beta.x<1e-20] <- 1e-20
all_PCs_r_best$pval_beta.y[all_PCs_r_best$pval_beta.y<1e-20] <- 1e-20

all_PCs_r_best_m <- melt(all_PCs_r_best)

    v <- transform(all_PCs_r_best_m, cluster=as.factor(cluster),variable=as.factor(variable))
    v <- v[!v$variable %in% c("qval.x","qval.y"),] %>%
   group_by(cluster,variable)%>%
   arrange(value) %>%
   mutate(observed=-log10(value), expected=-log10(ppoints(length(value))))
   v <- transform(v, variable=ifelse(variable=="pval_beta.x","ctrl_tensorQTL","alltreat_tensorQTL"))
   return(v)
}), data.frame)

p0 <- ggplot(vlist, aes(x=expected, y=observed, color=variable))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    facet_wrap(.~cluster, scales="free_y",ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    labs(color = "QTL method")+
    #ggtitle(paste0(unique(v$variable)))+
    theme_bw()
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/ctrltensor_vs_alltreattensorqtl_qqplot.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()


#bar plot
filenames <- list.files(ctrltensoroutFolder) #file list from directory
ctrltensorfilenames1 <- filenames[grep("_significant_topeeQTL_pairs.txt", filenames)] #pick specific files from list

filenames <- list.files(tensoroutFolder) #file list from directory
tensorfilenames1 <- filenames[grep("_significant_topeeQTL_pairs.txt", filenames)] #pick specific files from list

#I did check, the number of eGenes is the same as the number of SNP:eGene pairs
vlist <- ldply(lapply(celltypeorder$cluster.ctrlonly[!celltypeorder$cluster.ctrlonly %in% c("C9","C10")],function(c){
    cat("running",c,"\n")
    mc <- celltypeorder[celltypeorder$cluster.ctrlonly==c,]
ctrltensor_pc_signif_pairs <- fread(paste0(ctrltensoroutFolder,ctrltensorfilenames1[grep(paste0(c,".",treat),ctrltensorfilenames1)]))
ctrltensor_pc_signif_pairs <- transform(ctrltensor_pc_signif_pairs,method="ctrl_tensorQTL",cluster=c)
tensor_pc_signif_pairs <- fread(paste0(tensoroutFolder,tensorfilenames1[grep(paste0(mc$cluster.alltreat,".",treat),tensorfilenames1)]))
tensor_pc_signif_pairs <- transform(tensor_pc_signif_pairs,method="alltreat_tensorQTL",cluster=c)

df <- data.frame(cluster=c,ctrltensorQTL_eGene=length(unique(ctrltensor_pc_signif_pairs$phenotype_id)),
    tensorQTL_eGene=length(unique(tensor_pc_signif_pairs$phenotype_id)))
   return(df)
}), data.frame)

#justyna <- data.frame(cluster=c("C4","C5","C7"),justyna_fastqtl_eGene=c(364,291,253))
#vlistm <- merge(vlist,justyna,by="cluster",all=T)
outtablem <- melt(vlist)
p <- ggplot(outtablem, aes(fill=variable, y=value, x=cluster)) + 
    geom_bar(position="dodge", stat="identity")+
    labs(y="Number eGenes at FDR 10%",fill="QTL method")+
    theme_bw()
    #facet_wrap(.~description,ncol=4)
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output_CTRLonly/figures/ctrltensor_vs_tensorqtl_egene_bar.png")
png(width = 9, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()



data_names <- gsub("_significant_topeeQTL_pairs.txt", "", ctrltensorfilenames1) #remove file ending
for(i in 1:length(ctrltensorfilenames1)) assign(data_names[i], fread(file.path(ctrltensoroutFolder, ctrltensorfilenames1[i]),header = T)) #read in specific files and set the df object names. can dro0p unwanted columns
#combining data
all_PCs <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(all_PCs) <- data_names
all_PCs_r <- ldply(all_PCs, data.frame)

length(unique(all_PCs_r$phenotype_id)) #2974 -> 1926

#checked for slc27a3 (validated asthma gene)
grep("ENSG00000143554",all_PCs_r$phenotype_id) #not eqtl 

data_names <- gsub("_significant_topeeQTL_pairs.txt", "", tensorfilenames1) #remove file ending
for(i in 1:length(tensorfilenames1)) assign(data_names[i], fread(file.path(tensoroutFolder, tensorfilenames1[i]),header = T)) #read in specific files and set the df object names. can dro0p unwanted columns
#combining data
all_PCs <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(all_PCs) <- data_names
all_PCs_r <- ldply(all_PCs, data.frame)

length(unique(all_PCs_r$phenotype_id))
#3439
