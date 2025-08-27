#plotting old vs new deseq
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
baseoutFolder_dex=paste0(base,method,"_pseudobulk_ctrl/adjusted/",0.2,".",11,"/cell20filt/")

run="SES_PCs_sex_age_and_treats_generem"

outFolder_dex=paste0(baseoutFolder_dex,run,"/")
outFolder=paste0(baseoutFolder,run,"/")
celltype_dex <- data.frame(cluster.withdex=c("C0","C1","C2","C3","C4","C5","C6"),celltype=c("CD4_T","CD4_T","CD8_T","NK","Monocyte","B","DC"))
celltype <- data.frame(cluster.withoutdex=c("C0","C1","C2","C3","C4","C5","C6"),celltype=c("CD4_T","CD8_T","NK","Monocyte","B","DC","Monocyte"))

threshold=0.10
#treat="RNA-LPS"
i="RNA-CTRL"
lapply(c("PSS_all_mean","ISEL_Mean"), function(var){
        nodex_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",0.15,".",13,".deseqres_",var,"-",i,".",run,".txt"))
        dex_deseqres <- fread(paste0(outFolder_dex,"deseqres/",project,".",0.2,".",11,".deseqres_",var,"-",i,".",run,".txt"))
        nodex_deseqres <- merge(celltype, nodex_deseqres, by.y="cluster",by.x="cluster.withoutdex")
        dex_deseqres <- merge(celltype_dex, dex_deseqres, by.y="cluster",by.x="cluster.withdex")
#        nodex_deseqres <- transform(nodex_deseqres, model="noDEX")
#        dex_deseqres <- transform(dex_deseqres, model="DEX")
        subvars <- merge(dex_deseqres,nodex_deseqres,by=c("identifier","var","treats","celltype"))

subvars <- transform(subvars, dex_zscore=logFC.x/SE.x, nodex_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "2dex_sig", ifelse(padj.y<=threshold, "3nodex_sig", "1Not_Sig")))))
group_colors <- c("1Not_Sig" = "grey","2dex_sig" = "#E34234", "3nodex_sig" = "turquoise4", "4both" = "#7851A9")
p <- ggplot(subvars[order(subvars$sig),], aes(x=logFC.x, y=logFC.y)) +
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  labs(x="with Dex",y="no Dex")+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".dexvsnodex_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
#per cluster
p <- ggplot(subvars[order(subvars$sig),], aes(x=logFC.x, y=logFC.y)) +
  theme_bw()+
  facet_wrap(.~celltype,ncol=3)+
  geom_point(aes(color=sig))+ #aes(color=sig)
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  labs(x="with Dex",y="no Dex")+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".dexvsnodex_logfc_bycluster.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
})


        pss_nodex_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",0.15,".",13,".deseqres_","PSS_all_mean","-",i,".",run,".txt"))
        isel_nodex_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",0.15,".",13,".deseqres_","ISEL_Mean","-",i,".",run,".txt"))

var="ISEL_Mean"
allvardeseq <- ldply(lapply(c("C0","C1","C2","C3","C4","C5"),function(cluster) {
  cat("running",cluster,"\t")
    opfn <- paste0(outFolder,project,".",0.15,".",13,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
    load(opfn)
    res <- results(dds)
    nodex_deseqres <- data.frame(res@rownames, res$'baseMean', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
    names(nodex_deseqres) <- c('identifier', 'baseMean', 'pvalue', 'logFC','SE')
    nodex_deseqres <- transform(nodex_deseqres, padj=p.adjust(pvalue,method="BH"),var=var,cluster=cluster,treats=i,model="noDEX")
    opfn <- paste0(outFolder_dex,project,".",0.2,".",11,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
    load(opfn)
    res <- results(dds)
    dex_deseqres <- data.frame(res@rownames, res$'baseMean', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
    names(dex_deseqres) <- c('identifier', 'baseMean', 'pvalue', 'logFC','SE')
    dex_deseqres <- transform(dex_deseqres, padj=p.adjust(pvalue,method="BH"),var=var,cluster=cluster,treats=i,model="DEX")

    nodex_deseqres <- merge(nodex_deseqres, celltype, by.x="cluster",by.y="cluster.withoutdex")
    dex_deseqres <- merge(dex_deseqres, celltype_dex, by.x="cluster",by.y="cluster.withdex")

    subvars <- rbind(nodex_deseqres,dex_deseqres)
    fwrite(nodex_deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".noindfilt.txt"))
    fwrite(dex_deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder_dex,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".noindfilt.txt"))
    return(subvars)
    }),data.frame)

    merged_treat <- merge(allvardeseq[allvardeseq$model=="DEX",],allvardeseq[allvardeseq$model=="noDEX",],by=c("identifier","var","treats","celltype"))
    merged_treat <- transform(merged_treat, dex_zscore=logFC.x/SE.x, nodex_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "2dex_sig", ifelse(padj.y<=threshold, "3nodex_sig", "1Not_Sig")))))

merged_treat <- subset(merged_treat, !is.na(sig))
group_colors <- c("1Not_Sig" = "grey","2dex_sig" = "#E34234", "3nodex_sig" = "turquoise4", "4both" = "#7851A9")

p <- ggplot(merged_treat[order(merged_treat$sig),], aes(x=logFC.x, y=logFC.y)) +
  theme_bw()+
  facet_wrap(.~celltype,ncol=3)+
  geom_point(aes(color=sig))+ #aes(color=sig)
  scale_color_manual(values = group_colors) +
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  labs(x="with Dex",y="no Dex")+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".dexvsnodex_logfc_bycluster_noindfilter.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

library(tinyutils)
allclus <- lapply(c("C0","C1","C2","C3","C4","C5"),function(cluster) {
  cat("running",cluster,"\t")
    opfn <- paste0(outFolder,project,".",0.15,".",13,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
    load(opfn)
    p <- plotDispEsts2(dds)
    p <- p+ ggtitle(paste0("no DEX ",cluster))
    return(p)
    })
library(cowplot)
png(width = 12, height = 8, file=paste0(outFolder,"figures/",i,".",var,".",run,".",0.15,".",13,".nodex_dispests.png"), pointsize=12, 
    bg = "transparent", canvas = "white", units = "in", res = 600)
plot_grid(plotlist=allclus, ncol=3)
dev.off()

allclus <- lapply(c("C0","C1","C2","C3","C4","C5"),function(cluster) {
  cat("running",cluster,"\t")
    opfn <- paste0(outFolder_dex,project,".",0.2,".",11,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
    load(opfn)
    p <- plotDispEsts2(dds)
    p <- p+ ggtitle(paste0("DEX ",cluster))
    return(p)
    })
png(width = 12, height = 8, file=paste0(outFolder,"figures/",i,".",var,".",run,".",0.2,".",11,".dex_dispests.png"), pointsize=12, 
    bg = "transparent", canvas = "white", units = "in", res = 600)
plot_grid(plotlist=allclus, ncol=3)
dev.off()


celltype_combo <- data.frame(cluster.withdex=c("C0","C1","C2","C3","C4","C5","C6"),celltype.withdex=c("CD4_T","CD4_T","CD8_T","NK","Monocyte","B","DC"),
  cluster.withoutdex=c("C0","C1","C2","C3","C4","C5","C6"),celltype.withoutdex=c("CD4_T","CD8_T","NK","Monocyte1","B","DC","Monocyte2"))

celltype_combo <- data.frame(cluster.withdex=c("C0","C2","C3","C4","C5","C6","C1"),celltype.withdex=c("CD4_T","CD8_T","NK","Monocyte","B","DC","CD4_T"),
  cluster.withoutdex=c("C0","C1","C2","C3","C4","C5","C6"),celltype.withoutdex=c("CD4_T","CD8_T","NK","Monocyte","B","DC","Monocyte"))


var="ISEL_Mean"
allvardeseq <- ldply(lapply(c("C0","C1","C2","C3","C4"),function(cluster) {
    cellcluster <- subset(celltype_combo, cluster.withoutdex==cluster)
  cat("running",cluster,"\t")
    opfn <- paste0(outFolder,project,".",0.15,".",13,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
    load(opfn)
    nodex_metadata <- colData(dds)

    opfn <- paste0(outFolder_dex,project,".",0.2,".",11,".DESeq_output-",i,"-",var,cellcluster$cluster.withdex,".",run,".RData")
    load(opfn)
    dex_metadata <- colData(dds)

    nodex_not_dex <- tryCatch(data.frame(Sample_ID=nodex_metadata$Sample_ID[!nodex_metadata$Sample_ID %in% dex_metadata$Sample_ID], model="nodex_not_dex"),error=function(e){data.frame(Sample_ID=NA,model="nodex_not_dex")})
    dex_not_nodex <- tryCatch(data.frame(Sample_ID=dex_metadata$Sample_ID[!dex_metadata$Sample_ID %in% nodex_metadata$Sample_ID], model="dex_not_nodex"),error=function(e){data.frame(Sample_ID=NA,model="nodex_not_dex")})
    df <- rbind(nodex_not_dex,dex_not_nodex)
    df <- transform(df, cluster.withoutdex=cellcluster$cluster.withoutdex,cluster.withdex=cellcluster$cluster.withdex,celltype.withoutdex=cellcluster$celltype.withoutdex,celltype.withdex=cellcluster$celltype.withdex,var=var)
    return(df)
    }), data.frame)

allvardeseq_noDC <- subset(allvardeseq, !is.na(Sample_ID))
allvardeseq_noDC[order(allvardeseq_noDC$celltype.withoutdex),]

#which genes
library(ggpubr)
library(scales)

var="ISEL_Mean"
var="PSS_all_mean"
allvardeseq <- ldply(lapply(c("C0","C1","C2","C3","C4"),function(cluster) {
    cellcluster <- subset(celltype_combo, cluster.withoutdex==cluster)
  cat("running",cluster,"\t")
    opfn <- paste0(outFolder,project,".",0.15,".",13,".DESeq_output-",i,"-",var,cluster,".",run,".RData")
    load(opfn)
    res <- results(dds)
    nodex_deseqres <- data.frame(res@rownames, res$'baseMean', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
    names(nodex_deseqres) <- c('identifier', 'baseMean', 'pvalue', 'logFC','SE')
    nodex_deseqres <- transform(nodex_deseqres, padj=p.adjust(pvalue,method="BH"),var=var,cluster=cluster,treats=i,model="noDEX")
nodex_norm <- counts(dds, normalized=TRUE)
nodex_var <- colData(dds)[,var]

    opfn <- paste0(outFolder_dex,project,".",0.2,".",11,".DESeq_output-",i,"-",var,cellcluster$cluster.withdex,".",run,".RData")
    load(opfn)
    res <- results(dds)
    dex_deseqres <- data.frame(res@rownames, res$'baseMean', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
    names(dex_deseqres) <- c('identifier', 'baseMean', 'pvalue', 'logFC','SE')
    dex_deseqres <- transform(dex_deseqres, padj=p.adjust(pvalue,method="BH"),var=var,cluster=cluster,treats=i,model="DEX")
    dex_norm <- counts(dds, normalized=TRUE)
dex_var <- colData(dds)[,var]

#plot x=variable, y=cpm

nodex_norm_c <- transform(nodex_norm,gene_present=ifelse(rownames(nodex_norm) %in% nodex_deseqres$identifier & ! rownames(nodex_norm) %in% dex_deseqres$identifier,"nodex_not_dex",ifelse(rownames(nodex_norm) %in% dex_deseqres$identifier & ! rownames(nodex_norm) %in% nodex_deseqres$identifier,"dex_not_nodex","both")))
nodex_norm_c_i <- transform(nodex_norm_c, identifier=rownames(nodex_norm_c))
nodex_norm_cm <- melt(nodex_norm_c_i)
nodex_norm_cm %>% split(.$gene_present) %>% map(summary)
subset(nodex_norm_cm, !gene_present=="both" & value>2000)

dex_norm_c <- transform(dex_norm,gene_present=ifelse(rownames(dex_norm) %in% nodex_deseqres$identifier & ! rownames(dex_norm) %in% dex_deseqres$identifier,"nodex_not_dex",ifelse(rownames(dex_norm) %in% dex_deseqres$identifier & ! rownames(dex_norm) %in% nodex_deseqres$identifier,"dex_not_nodex","both")))
dex_norm_c_i <- transform(dex_norm_c, identifier=rownames(dex_norm_c))
dex_norm_cm <- melt(dex_norm_c)
dex_norm_cm %>% split(.$gene_present) %>% map(summary)
subset(dex_norm_cm, !gene_present=="both" & value>2000)

dex_plotdf <- data.frame(x=as.vector(dex_var),y=as.vector(dex_norm),gene_present=ifelse(rownames(dex_norm) %in% nodex_not_dex$identifier,"nodex_not_dex",ifelse(rownames(dex_norm) %in% dex_not_nodex$identifier,"dex_not_nodex","both")))

    nodex_not_dex <- tryCatch(data.frame(identifier=nodex_deseqres$identifier[!nodex_deseqres$identifier %in% dex_deseqres$identifier], model="nodex_not_dex"),error=function(e){data.frame(identifier=NA,model="nodex_not_dex")})
    dex_not_nodex <- tryCatch(data.frame(identifier=dex_deseqres$identifier[!dex_deseqres$identifier %in% nodex_deseqres$identifier], model="dex_not_nodex"),error=function(e){data.frame(identifier=NA,model="nodex_not_dex")})
    df <- rbind(nodex_not_dex,dex_not_nodex)
    df <- transform(df, cluster.withoutdex=cellcluster$cluster.withoutdex,cluster.withdex=cellcluster$cluster.withdex,celltype.withoutdex=cellcluster$celltype.withoutdex,celltype.withdex=cellcluster$celltype.withdex,var=var)
    return(df)
    }), data.frame)

allvardeseq_noDC <- subset(allvardeseq, !is.na(identifier))
allvardeseq_noDC[order(allvardeseq_noDC$celltype.withoutdex),]

table(allvardeseq_noDC$celltype.withoutdex, allvardeseq_noDC$model)


nodex_norm <- counts(dds, normalized=TRUE)
nodex_var <- colData(dds)[,var]
#plot x=variable, y=cpm
nodex_plotdf <- data.frame(x=as.vector(nodex_var),y=as.vector(nodex_norm),gene_present=ifelse(rownames(nodex_norm) %in% nodex_not_dex$identifier,"nodex_not_dex",ifelse(rownames(nodex_norm) %in% dex_not_nodex$identifier,"dex_not_nodex","both")))


p <- ggplot(plotdf, aes(x=x, y=y,color=gene_present)) +
  theme_bw()+
  geom_point()+ #aes(color=sig)
    #geom_abline(slope=1,intercept=0) +
  geom_smooth(method = "lm", fill = NA)+
  ylab(paste0(gene," normalized")) +
  xlab(paste0(v)) +
  ggtitle(paste0(c," ",i," ",v,"\n","DESeq logFC=",round(d_info$logFC,2)," DESeq padj=",scientific(d_info$padj, digits = 3)," DESeq pvalue=",scientific(d_info$pvalue, digits = 3),
    "\n","limma logFC=",round(l_info$logFC,2)," limma padj=",scientific(l_info$adj.P.Val, digits = 3)," limma pvalue=",scientific(l_info$P.Value, digits = 3)))+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),plot.title = element_text(size = rel(1.1), face = "bold"),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 6, height = 6, file=paste0(outFolder,"figures/",c,".",i,".",v,".",gene,".png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(p)
dev.off()

})