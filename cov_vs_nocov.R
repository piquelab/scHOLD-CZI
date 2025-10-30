library(Vennerable)

#want to plot isel vs no isel in DESeq model
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
#run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
run="noCombat_DESeq"
treatmentsfirst=c("RNA-CTRL","RNA-LPS")

outFolder=paste0(baseoutFolder,run,"/")
sesrun=paste0(run,"_iselcov")
var="PSS_all_mean"
threshold=0.10

for (i in treatmentsfirst){
lapply(c("PSS_all_mean"),function(var){
    allvardeseq <- ldply(lapply(names(counts_ls),function(cluster) {
      cat("running",cluster,"\t")
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
        subvars <- transform(merged_treat, noisel_zscore=logFC.x/SE.x, isel_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "1both", ifelse(padj.x<=threshold, "3noisel_sig", ifelse(padj.y<=threshold, "2isel_sig", "4Not_Sig")))))
        fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".noindfilt.txt"))
        fwrite(ses_deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",sesrun,".noindfilt.txt"))
        return(subvars)
        }),data.frame)
allvardeseq <- subset(allvardeseq, !is.na(sig))
group_colors <- c("4Not_Sig" = "grey","3noisel_sig" = "#E34234", "2isel_sig" = "turquoise4", "1both" = "#7851A9")

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
png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".iselvsnoisel.noindfilt.png"), pointsize=12, 
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
png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".iselvsnoisel_bycluster.noindfilt.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

noSESressig <- subset(allvardeseq[,c(1:9)],padj.x<0.1)
names(noSESressig)[9] <- "model"
SESressig <- subset(allvardeseq[,c(1:4,10:14)],padj.y<0.1)
names(SESressig)[9] <- "model"

if(dim(noSESressig)[1]>0 & dim(SESressig)[1]>0){
df <- unique(rbind(noSESressig[,c("identifier","model")],SESressig[,c("identifier","model")]))
d_list <- split(df$identifier,df$model)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0(outFolder,"figures/",project,".",var,"-",i,".noiselvsisel_venn.noindfilt.png"), pointsize=12, 
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
png(width = 12, height = 7, file=paste0(outFolder,"figures/noiselvsisel",".",var,"-",i,".stackedbar.noindfilt.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(fig0)
dev.off()
})
}


#want to plot PSS vs no PSS in DESeq model
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
#run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
run="noCombat_DESeq"

outFolder=paste0(baseoutFolder,run,"/")
sesrun=paste0(run,"_PSScov")
var="ISEL_Mean"
threshold=0.10

for (i in treatmentsfirst){
lapply(c("ISEL_Mean"),function(var){
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
        deseqres <- transform(deseqres, model="noPSSCov")
        ses_deseqres <- transform(ses_deseqres, model="PSSCov")
        merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var","treats"))
        subvars <- transform(merged_treat, noPSS_zscore=logFC.x/SE.x, PSS_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "1both", ifelse(padj.x<=threshold, "3noPSS_sig", ifelse(padj.y<=threshold, "2PSS_sig", "4Not_Sig")))))
        fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".noindfilt.txt"))
        fwrite(ses_deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",sesrun,".noindfilt.txt"))  
        return(subvars)
        }),data.frame)
allvardeseq <- subset(allvardeseq, !is.na(sig))
group_colors <- c("4Not_Sig" = "grey","3noPSS_sig" = "#E34234", "2PSS_sig" = "turquoise4", "1both" = "#7851A9")
p <- ggplot(allvardeseq, aes(x=PSS_zscore, y=noPSS_zscore)) +
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
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".PSSvsnoPSS.noindfilt.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
#per cluster
p <- ggplot(allvardeseq, aes(x=PSS_zscore, y=noPSS_zscore)) +
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
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".PSSvsnoPSS_bycluster.noindfilt.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

noSESressig <- subset(allvardeseq[,c(1:9)],padj.x<0.1)
names(noSESressig)[9] <- "model"
SESressig <- subset(allvardeseq[,c(1:4,10:14)],padj.y<0.1)
names(SESressig)[9] <- "model"

if(dim(noSESressig)[1]>0 & dim(SESressig)[1]>0){
df <- unique(rbind(noSESressig[,c("identifier","model")],SESressig[,c("identifier","model")]))
d_list <- split(df$identifier,df$model)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0(outFolder,"figures/",project,".",var,"-",i,".noPSSvsPSS_venn.noindfilt.png"), pointsize=12, 
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

group_colors <- c("3noPSS_sig" = "#E34234", "2PSS_sig" = "turquoise4", "4both" = "#7851A9")

fig0 <- ggplot(count_df, aes(x = count, y = cluster, fill = sig)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group_colors) +
  theme_bw() +
    #facet_grid(.~term)+
  labs(
    x = "Number of DEGs",
    y = "Cell Type",
    fill = "Group",
    title = str_wrap("DEGs Overlap Between PSS and no PSS covariate model", width = 30)
  ) +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, lineheight = 1.1),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  )
# Save the plot
png(width = 12, height = 7, file=paste0(outFolder,"figures/noPSSvsPSS",".",var,"-",i,".stackedbar.noindfilt.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(fig0)
dev.off()
})
}

#want to plot WHR vs no WHR in DESeq model
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
run="noCombat_DESeq"

#run="SES_PCs_sex_age_and_treats_generem"
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
    lapply(c("cytocomp","PSS_all_mean","ISEL_Mean"),function(var){
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

#want to plot SES vs no SES in DESeq model
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
#run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
run="noCombat_DESeq"

outFolder=paste0(baseoutFolder,run,"/")
sesrun=paste0(run,"_SEScov")
var="ISEL_Mean"
threshold=0.10

for (i in treatmentsfirst){
lapply(c("cytocomp","ISEL_Mean","PSS_all_mean"),function(var){
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
        deseqres <- transform(deseqres, model="noSESCov")
        ses_deseqres <- transform(ses_deseqres, model="SESCov")
        merged_treat <- merge(deseqres,ses_deseqres,by=c("identifier","cluster","var","treats"))
        subvars <- transform(merged_treat, noSES_zscore=logFC.x/SE.x, SES_zscore=logFC.y/SE.y, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "1both", ifelse(padj.x<=threshold, "3noSES_sig", ifelse(padj.y<=threshold, "2SES_sig", "4Not_Sig")))))
        fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",run,".noindfilt.txt"))
        fwrite(ses_deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",var,"-",i,".",sesrun,".noindfilt.txt"))  
        return(subvars)
        }),data.frame)
allvardeseq <- subset(allvardeseq, !is.na(sig))
group_colors <- c("4Not_Sig" = "grey","3noSES_sig" = "#E34234", "2SES_sig" = "turquoise4", "1both" = "#7851A9")
p <- ggplot(allvardeseq, aes(x=SES_zscore, y=noSES_zscore)) +
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
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".SESvsnoSES.noindfilt.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()
#per cluster
p <- ggplot(allvardeseq, aes(x=SES_zscore, y=noSES_zscore)) +
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
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",i,".",var,".",run,".SESvsnoSES_bycluster.noindfilt.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

noSESressig <- subset(allvardeseq[,c(1:9)],padj.x<0.1)
names(noSESressig)[9] <- "model"
SESressig <- subset(allvardeseq[,c(1:4,10:14)],padj.y<0.1)
names(SESressig)[9] <- "model"

if(dim(noSESressig)[1]>0 & dim(SESressig)[1]>0){
df <- unique(rbind(noSESressig[,c("identifier","model")],SESressig[,c("identifier","model")]))
d_list <- split(df$identifier,df$model)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0(outFolder,"figures/",project,".",var,"-",i,".noSESvsSES_venn.noindfilt.png"), pointsize=12, 
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

group_colors <- c("3noSES_sig" = "#E34234", "2SES_sig" = "turquoise4", "4both" = "#7851A9")

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
png(width = 12, height = 7, file=paste0(outFolder,"figures/noSESvsSES",".",var,"-",i,".stackedbar.noindfilt.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 500)
print(fig0)
dev.off()
})
}



#combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
#plotting LPS vs CONTROL
#baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/cell20filt/")
#run="treatment_withCOMBAT"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/nodex/")
combat="noCOMBAT"
run=paste0("treatment_",combat)

outFolder=paste0(baseoutFolder,run,"/")
threshold=0.10
treat="RNA-LPS"
control="RNA-CTRL"
#contrastdf <- data.frame(control=c("RNA-CTRL","RNA-LPS"),treatment=c("RNA-LPS","RNA-LPS-DEX"))
contrastdf <- data.frame(control=c("RNA-CTRL"),treatment=c("RNA-LPS"))
contrastdf <- transform(contrastdf, contrast=paste0(treatment,"_vs_",control))
contrast=contrastdf[1,]$contrast

#varrun="SES_PCs_sex_age_and_treats_generem"
varrun="noCombat_DESeq"
varoutFolder=paste0(baseoutFolder,varrun,"/")

library(ggrepel)

lapply(c("PSS_all_mean","ISEL_Mean"), function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_deseqres <- fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",control,".",varrun,".txt"))
        ctrl_deseqres_sig <- subset(ctrl_deseqres, padj<0.1)
        lpsvctrl <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
        lpsvctrl <- transform(lpsvctrl, treatsig=ifelse(padj<0.1 & abs(logFC)>2,"sig","not_sig"),varDEG=ifelse(identifier %in% ctrl_deseqres_sig$identifier,paste0(var,"_sig"),"not_sig"))
        lpsvctrl <- transform(lpsvctrl, bothsig=ifelse(treatsig=="sig" & !varDEG=="not_sig","bothsig", ifelse(!varDEG=="not_sig" & !treatsig=="sig", paste0(var,"_sig"), ifelse(treatsig=="sig" & varDEG=="not_sig", paste0(contrast,"_sig"),"not_sig"))))
        lpsvctrl$bothsig <- factor(lpsvctrl$bothsig, levels = c("bothsig",paste0(var,"_sig"),paste0(contrast,"_sig"),"not_sig"))
        tab <- subset(lpsvctrl, bothsig=="bothsig" & abs(logFC)>2)
p <- ggplot(lpsvctrl[order(lpsvctrl$bothsig),], aes(baseMean, logFC, col=bothsig, label=identifier)) + 
        geom_point() + scale_x_log10() + 
        xlab("mean of normalized counts") + 
        ylab("log fold change") + 
        geom_hline(yintercept=0, col="grey40") + 
        scale_color_manual(values=c("purple", "red", "blue", "grey60") |> `names<-`(c("bothsig", paste0(var,"_sig"),paste0(contrast,"_sig"),"not_sig"))) +
        #geom_point(data=tab, shape=1, size=3, show.legend=FALSE) + 
        geom_label_repel(data=tab, nudge_x = 1, nudge_y = 2*sign(tab$logFC), show.legend=FALSE)

png(width = 6, height = 6, file=paste0(outFolder,"figures/",var,"sig-",contrast,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(p)
dev.off()
})
#single plot, cat=PSS, ISEL, both, none
subvarsls <- lapply(c("PSS_all_mean","ISEL_Mean"), function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_deseqres <- fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",control,".",varrun,".txt"))
        ctrl_deseqres_sig <- subset(ctrl_deseqres, padj<0.1)
        lpsvctrl <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
        lpsvctrl <- transform(lpsvctrl, treatsig=ifelse(padj<0.1 & abs(logFC)>2,"sig","not_sig"),varDEG=ifelse(identifier %in% ctrl_deseqres_sig$identifier,"sig","not_sig"))
        return(lpsvctrl)
        })
names(subvarsls) <- c("PSS_all_mean","ISEL_Mean")
lpsvctrl <- as.data.frame(Reduce(function(x, y) merge(x, y,by=c("identifier","cluster","contrast","treatsig","baseMean","logFC")),subvarsls,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
names(lpsvctrl)[c(10,14)] <- c("PSS_all_mean_sig","ISEL_Mean_sig")
lpsvctrl <- transform(lpsvctrl, bothsig=ifelse(PSS_all_mean_sig=="sig" & ISEL_Mean_sig=="sig","bothsig", ifelse(PSS_all_mean_sig=="sig" & !ISEL_Mean_sig=="sig", "PSS_all_mean_sig", ifelse(!PSS_all_mean_sig=="sig" & ISEL_Mean_sig=="sig", "ISEL_Mean_sig","not_var_sig"))))
lpsvctrl$bothsig <- factor(lpsvctrl$bothsig, levels = c("bothsig","PSS_all_mean_sig","ISEL_Mean_sig","not_var_sig"))
lpsvctrl$bothsig <- factor(lpsvctrl$bothsig, levels = c("not_var_sig","ISEL_Mean_sig","PSS_all_mean_sig","bothsig"))
tab <- subset(lpsvctrl, bothsig=="bothsig" & abs(logFC)>2)

p <- ggplot(lpsvctrl[order(lpsvctrl$bothsig),], aes(baseMean, logFC, col=bothsig, label=identifier)) + 
        geom_point() + scale_x_log10() + 
        xlab("mean of normalized counts") + 
        ylab("log fold change") + 
        geom_hline(yintercept=0, col="grey40") + 
        geom_hline(yintercept=-2, col="grey40",linetype="dashed") + 
        geom_hline(yintercept=2, col="grey40",linetype="dashed") + 
        scale_color_manual(values=c("purple", "red", "blue", "grey60") |> `names<-`(c("bothsig", "PSS_all_mean_sig","ISEL_Mean_sig","not_var_sig"))) +
        #geom_point(data=tab, shape=1, size=3, show.legend=FALSE) + 
        geom_label_repel(data=tab, nudge_x = 1, nudge_y = 2*sign(tab$logFC), show.legend=FALSE)

png(width = 6, height = 6, file=paste0(outFolder,"figures/","pss_isel_comb-",contrast,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(p)
dev.off()

#separate clusters and for the colors maybe we can categorize up vs down regulated and do each variable separately 
lapply(c("PSS_all_mean","ISEL_Mean"), function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_deseqres <- fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",control,".",varrun,".txt"))
        lpsvctrl_a <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
  subvars <- ldply(lapply(names(counts_ls),function(c){
      lpsvctrl <- subset(lpsvctrl_a, cluster==c)
      ctrl_deseqres_sig <- subset(ctrl_deseqres, padj<0.1 & cluster==c)
      ctrl_deseqres_upsig <- subset(ctrl_deseqres_sig, logFC>0 )
      ctrl_deseqres_downsig <- subset(ctrl_deseqres_sig, logFC<0 )
      lpsvctrl <- transform(lpsvctrl, varDEG=ifelse(identifier %in% ctrl_deseqres_upsig$identifier,"upsig",ifelse(identifier %in% ctrl_deseqres_downsig$identifier,"downsig","not_var_sig")))
      lpsvctrl$varDEG <- factor(lpsvctrl$varDEG, levels = c("not_var_sig","upsig","downsig"))
      return(lpsvctrl)
      }), data.frame)

p <- ggplot(subvars[order(subvars$varDEG),], aes(baseMean, logFC, col=varDEG, label=identifier)) + 
        geom_point() + scale_x_log10() + 
        xlab("mean of normalized counts") + 
        ylab("log fold change") + 
        facet_wrap(.~cluster,ncol=3)+
        ggtitle(paste0(var))+
        geom_hline(yintercept=0, col="grey40") + 
        geom_hline(yintercept=-2, col="grey40",linetype="dashed") + 
        geom_hline(yintercept=2, col="grey40",linetype="dashed") + 
        scale_color_manual(values=c("red", "blue", "grey60") |> `names<-`(c("upsig","downsig","not_var_sig"))) 
        #geom_point(data=tab, shape=1, size=3, show.legend=FALSE) + 

png(width = 6, height = 6, file=paste0(outFolder,"figures/",var,"_comb-",contrast,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(p)
dev.off()
})

#MAplot of treatment only no var
lpsvctrl_a <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres.",contrast,".",run,".txt"))
subvars <- ldply(lapply(names(counts_ls),function(c){
    lpsvctrl <- subset(lpsvctrl_a, cluster==c)
    lpsvctrl <- transform(lpsvctrl, treatDEG=ifelse(logFC>0 & padj<0.1,"upsig",ifelse(logFC<0 & padj<0.1,"downsig","not_sig")))
    lpsvctrl$treatDEG <- factor(lpsvctrl$treatDEG, levels = c("not_sig","upsig","downsig"))
    return(lpsvctrl)
    }), data.frame)

p <- ggplot(subvars[order(subvars$treatDEG),], aes(baseMean, logFC, col=treatDEG, label=identifier)) + 
        geom_point() + scale_x_log10() + 
        xlab("mean of normalized counts") + 
        ylab("log fold change") + 
        facet_wrap(.~cluster,ncol=3)+
        ggtitle(paste0(contrast))+
        geom_hline(yintercept=0, col="grey40") + 
        geom_hline(yintercept=-2, col="grey40",linetype="dashed") + 
        geom_hline(yintercept=2, col="grey40",linetype="dashed") + 
        scale_color_manual(values=c("red", "blue", "grey60") |> `names<-`(c("upsig","downsig","not_sig"))) 
        #geom_point(data=tab, shape=1, size=3, show.legend=FALSE) + 

png(width = 6, height = 6, file=paste0(outFolder,"figures/",contrast,".",run,"_MAplot.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(p)
dev.off()



# expression in treatment vs expression in control, and then colored by deg for psychosocial factors. average across samples? and have one dot per gene. log transformed values


cluster="C0"
cpmb_ls <- ldply(lapply(names(counts_ls),function(c){
PSSctrl_deseqres_sig <- subset(fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_","PSS_all_mean","-",control,".",varrun,".txt")),padj<0.1 & cluster==c)
iselctrl_deseqres_sig <- subset(fread(paste0(varoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_","ISEL_Mean","-",control,".",varrun,".txt")),padj<0.1 & cluster==c)

cluster_counts_sce <- counts_ls[[c]]
cluster_counts <- assay(cluster_counts_sce, "counts")
#opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",c,".",combatrun,".RData")
#load(opfn)
cluster_metadata_sce <- metadata_ls[[c]]
cluster_metadata <- data.frame(cluster_metadata_sce)
cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats),rowid=rownames(cluster_metadata))
cluster_metadata_c <- subset(cluster_metadata, treats==control)
cluster_metadata_t <- subset(cluster_metadata, treats==treat)
#cluster_counts_c <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_c))]
#cluster_counts_t <- adjusted_counts[,which(colnames(adjusted_counts) %in% rownames(cluster_metadata_t))]
cluster_counts_c <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_c))]
cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_t))]

dget <- DGEList(counts=cluster_counts_t)
dget <- calcNormFactors(dget)
cpmt <- log2(cpm(dget)+1)
rownames(cpmt) <- rownames(cluster_counts_t)
cpmt_avg <- as.data.frame(rowMeans(cpmt, na.rm = TRUE))
cpmt_avg$identifier <- rownames(cpmt_avg)
dgec <- DGEList(counts=cluster_counts_c)
dgec <- calcNormFactors(dgec)
cpmc <- log2(cpm(dgec)+1)
rownames(cpmc) <- rownames(cluster_counts_c)
cpmc_avg <- as.data.frame(rowMeans(cpmc, na.rm = TRUE))
cpmc_avg$identifier <- rownames(cpmc_avg)

cpmb <- merge(cpmt_avg,cpmc_avg,by="identifier")
names(cpmb)[c(2,3)] <- c(paste0("logcpm_",treat),paste0("logcpm_",control))
cpmb <- transform(cpmb, cluster=c, bothsig=ifelse(identifier %in% PSSctrl_deseqres_sig$identifier & identifier %in% iselctrl_deseqres_sig$identifier,"bothsig", ifelse(identifier %in% PSSctrl_deseqres_sig$identifier, "PSS_all_mean_sig", ifelse(identifier %in% iselctrl_deseqres_sig$identifier, "ISEL_Mean_sig","not_var_sig"))))
cpmb$bothsig <- factor(cpmb$bothsig, levels = c("not_var_sig","ISEL_Mean_sig","PSS_all_mean_sig","bothsig"))
return(cpmb)
}),data.frame)

p <- ggplot(cpmb_ls[order(cpmb_ls$bothsig),], aes(logcpm_RNA.CTRL, logcpm_RNA.LPS, col=bothsig)) + 
    geom_point()+
    facet_wrap(.~cluster,ncol=3)+
    geom_abline(slope=1)+
    scale_color_manual(values=c("purple", "red", "blue", "grey60") |> `names<-`(c("bothsig", "PSS_all_mean_sig","ISEL_Mean_sig","not_var_sig"))) +
png(width = 6, height = 6, file=paste0(outFolder,"figures/","pss_isel_comb-",contrast,"_logcpm.",run,"_scatter.png"), pointsize=12, bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(p)
dev.off()

