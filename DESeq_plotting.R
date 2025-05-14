library(cowplot)
library(EnhancedVolcano)
library(ggcorrplot)
library(ggpubr)
library(data.table)
library(plyr)
library(parallel)
library(DESeq2)


args <- commandArgs(trailingOnly = TRUE)
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_chronic_cond_added_02_07_2024.txt","ALL","fastdemux") #for testing
base <- args[1]
cov_file=fread(args[2]) #t5 is the psych cov file
project=args[3]
method=args[4]
run="SES_PCs_SES_sex_age_and_treats_adjusted"
outFolder=paste0(base,method,"_pseudobulk_ctrl/")
figuredir=paste0(outFolder,"figures/")

filenames <- list.files(outFolder) #file list from directory
filenames1 <- filenames[grep("DESeq_output-", filenames)] #pick specific files from list
x <- na.omit(sapply(strsplit(filenames1,"(?<=.)(?=C[0-9])",perl=TRUE),function(y)y[2]))
clusters <- unique(sapply(strsplit(x,"[.]"),function(z)z[1]))
eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
eigenvec2 <- merge(eigenvec2_o,cov_file,by.x="Sample_ID",by.y="dbgap.ID",all.x=T)
notrun_var <- c("DSES_01","DSES_03","PWaist","PHip")
colnumuotovar <- grep("czi_exp",colnames(eigenvec2))+1
psychvarstorun <- eigenvec2[,colnumuotovar:length(colnames(eigenvec2))]
psychvarstorun <- colnames(psychvarstorun)[!colnames(psychvarstorun) %in% notrun_var]

cluster=clusters[1]
 var="LPSvsCTRL"
opfn <- paste0(outFolder,project,".DESeq_countlists.RData")
load(opfn)

load(paste0(outFolder,project,".DESeq_output-",var,cluster,".RData"))
annotation <- metadata_ls[[cluster]][match(rownames(colData(dds)), rownames(metadata_ls[[cluster]])),]

#mcols(dds) <- cbind(mcols(dds), metadata_ls[[cluster]])
colData(dds) <- cbind(colData(dds), annotation)

#reorder_idx = match(rownames(colData(dds)),rownames(metadata_ls[[cluster]])  )
#reordered <- metadata_ls[[cluster]][reorder_idx,]  

    vsd <- vst(dds, blind=FALSE)



    # Plot PCA
    fig0 <- DESeq2::plotPCA(vsd, intgroup = c("treats"))+ggtitle("treats")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig1 <- DESeq2::plotPCA(vsd, intgroup = c("sex_alph"))+ggtitle("sex")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig2 <- DESeq2::plotPCA(vsd, intgroup = c("age"))+ggtitle("age")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig3 <- DESeq2::plotPCA(vsd, intgroup = c("BATCH"))+ggtitle("BATCH")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig <-plot_grid(fig0, fig1,fig2,fig3, nrow=4, ncol=1, align="h")
    png(paste0(figuredir,project,".deseqPCA-treats",var,"-",cluster,".png"), width=1000, height=1000, res=120)
    print(fig)
    dev.off()

pca <-DESeq2::plotPCA(vsd, intgroup = c("treats"),returnData=TRUE)
outlier <- pca[pca$treats=="RNA-CTRL" & pca$PC1<(-20) & pca$PC2<0,]


for (cluster in clusters){
    #cluster=clusters[1]
    cat("running ", cluster, "\n")
    var="LPSvsCTRL"

    sub.table <- fread(paste0(outFolder,project,".deseqres_",var,".txt"))
    colnames(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','var','Cluster')
    sub.table <- subset(sub.table, Cluster==cluster)
    sub.table <- sub.table[order(sub.table$padj,-abs(sub.table$logFC)), ]

    topp1 <- min(subset(sub.table, padj > quantile(padj, prob = 1 -99/100,na.rm=T))$pvalue)
    tolab <- c(unique(head(sub.table,n=20)$identifier))
    png(width = 8, height = 8, file=paste0(figuredir,project,".dge_volcano-",var,"-",cluster,".png"), pointsize=12, 
          bg = "transparent", units = "in", res = 1200)
    p <- EnhancedVolcano(sub.table,
      lab = sub.table$identifier,
      x = 'logFC',
      y = 'pvalue',
      title=paste0(var,"-",cluster),
      subtitle = NULL,
      xlim = c(-10, 10),
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
                #run covariates separately for each treatment condition
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
        #i <- "RNA-CTRL"
        # Transform counts for data visualization
        lapply(c("sex","age","sex_age_int",psychvarstorun,"factor_HS_CRP"),function(var){
        #mclapply(psychvarstorun,function(var){
        #var="sex"
        #load(paste0(outFolder,"DESeq_output-",i,"-",var,cluster,".RData"))
        if(isTRUE(file.size(paste0(outFolder,project,".deseqres_",var,"-",i,".",run,".txt")) > 0)){

        cat("running ",var,i)
        sub.table <- fread(paste0(outFolder,project,".deseqres_",var,"-",i,".",run,".txt"))
        sub.table <- subset(sub.table, cluster==cluster)
        sub.table <- sub.table[order(sub.table$padj,-abs(sub.table$logFC)), ]

        topp1 <- min(subset(sub.table, padj > quantile(padj, prob = 1 -99/100,na.rm=T))$pvalue)
        tolab <- c(unique(head(sub.table,n=20)$identifier))
        png(width = 8, height = 8, file=paste0(figuredir,project,".dge_volcano-",var,"-",cluster,"-",i,".",run,".png"), pointsize=12, 
              bg = "transparent", units = "in", res = 1200)
        p <- EnhancedVolcano(sub.table,
          lab = sub.table$identifier,
          x = 'logFC',
          y = 'pvalue',
          title=paste0(var,"-",cluster,"-",i),
          subtitle = NULL,
          xlim = c(-10, 10),
          pCutoffCol= 'padj',
          pCutoff = 0.1,
          FCcutoff = 0.5,
            #pointSize = 1.5,
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
        })
}}


var="LPSvsCTRL"
dge <- fread(paste0(outFolder,"deseqres_",var,".txt"))
colnames(dge) <- c('gene_symbol', 'padj', 'pvalue', 'logFC','var','Cluster')

ddf <- split(dge,dge$Cluster)
ddf10 <- lapply(ddf, function(i) {
  i <- i[order(i$padj,-abs(i$logFC)),]
  return(head(i,n=5))
  })
top10 <- ldply(ddf10,data.frame)[,-1]

#top10 <- transform(top10, gene_symbol=ifelse(is.na(gene_symbol),ensembl_gene_id,gene_symbol),contrast=as.factor(contrast))
top10m <- reshape2::dcast(gene_symbol ~ Cluster, data=top10, value.var="logFC",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$gene_symbol
myMat <- top10m[,-1]

top10m <- reshape2::dcast(gene_symbol ~ Cluster, data=top10, value.var="padj",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$gene_symbol
pMat <- top10m[,-1]

myMat[is.na(pMat)] <- NA
pMat[is.na(myMat)] <- NA

png(width = 17, height = 8, file=paste0(figuredir,var,".dge_top_heatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2, insig = "blank",p.mat = pMat, sig.level=0.1,digits = 1) +
scale_fill_gradient2(low = "blue", high =  "red", mid = "white", midpoint = 0)
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
dev.off()

#taking into account direction for shared / not shared
dge_matching <- ldply(lapply(split(dge,dge$gene_symbol),function(i){
  #i <- split(dge,dge$gene_symbol)[[1]]
  logfc=as.vector(i$logFC)
  matching_dir=ifelse(all(logfc>0)|all(logfc<0),TRUE,FALSE)
  matching=ifelse(matching_dir==TRUE & length(logfc)==length(clusters),TRUE,
    ifelse(matching_dir==FALSE & length(logfc)==length(clusters),FALSE, NA))
  i <- transform(i, matching_direction=matching)
  return(i)
  }),data.frame)

dge_sub <- subset(dge_matching, matching_direction==TRUE)

ddf <- split(dge_sub,dge_sub$Cluster)
ddf10 <- lapply(ddf, function(i) {
  i <- i[order(i$padj,-abs(i$logFC)),]
  return(head(i,n=10))
  })
top10 <- ldply(ddf10,data.frame)[,-1]
dge_top <- dge_sub[dge_sub$gene_symbol %in% top10$gene_symbol, ]

top10m <- reshape2::dcast(gene_symbol ~ Cluster, data=dge_top, value.var="logFC",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$gene_symbol
myMat <- top10m[,-1]

top10m <- reshape2::dcast(gene_symbol ~ Cluster, data=dge_top, value.var="padj",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$gene_symbol
pMat <- top10m[,-1]

myMat[is.na(pMat)] <- NA
pMat[is.na(myMat)] <- NA

png(width = 17, height = 8, file=paste0(figuredir,var,".dge_topshared_heatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2, insig = "blank",p.mat = pMat, sig.level=0.1,digits = 1) +
scale_fill_gradient2(low = "blue", high =  "red", mid = "white", midpoint = 0)
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
dev.off()

dge_sub <- subset(dge_matching, matching_direction==FALSE)

ddf <- split(dge_sub,dge_sub$Cluster)
ddf10 <- lapply(ddf, function(i) {
  i <- i[order(i$padj,-abs(i$logFC)),]
  return(head(i,n=10))
  })
top10 <- ldply(ddf10,data.frame)[,-1]
dge_top <- dge_sub[dge_sub$gene_symbol %in% top10$gene_symbol, ]

top10m <- reshape2::dcast(gene_symbol ~ Cluster, data=dge_top, value.var="logFC",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$gene_symbol
myMat <- top10m[,-1]

top10m <- reshape2::dcast(gene_symbol ~ Cluster, data=dge_top, value.var="padj",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$gene_symbol
pMat <- top10m[,-1]

myMat[is.na(pMat)] <- NA
pMat[is.na(myMat)] <- NA

png(width = 17, height = 8, file=paste0(figuredir,var,".dge_topnotshared_heatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2, insig = "blank",p.mat = pMat, sig.level=0.1,digits = 1) +
scale_fill_gradient2(low = "blue", high =  "red", mid = "white", midpoint = 0)
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
dev.off()


#plot genotype PC by discrimination -- does genotype removal from the model remove the discrimination
library(ggpubr)

p <- ggplot(eigenvec2, aes(x=PC1, y=EDS_mean)) +
  theme_bw()+
  geom_point()+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
  ylab("Discrimination") +
  xlab("Genotype PC1") +
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))

png(width = 12, height = 12, file=paste0(figuredir,project,".PC1_vs_EDS.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

#CRP values , highlight HO-003, HO-123, HO-183
eigenvec2$all <- "All"
p <- ggplot(eigenvec2, aes(x=all,y=HS_CRP)) +
  theme_bw()+
  geom_boxplot()+ #aes(color=sig)
  geom_point(data=subset(eigenvec2, Sample_ID %in% c("HO-003", "HO-123", "HO-183")),aes(color="red"))+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))

png(width = 12, height = 12, file=paste0(figuredir,project,".CRP_box.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

#why are only higher clusters (fewer cells) showing DEGs
#zscore all clusters for PHA cdres
run="income_PCs_sex_age_and_treats_adjusted_withWave"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
outFolder=paste0(baseoutFolder,run,"/")
figuredir=paste0(outFolder,"figures/")
i="PHA"
var="cdres"
deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",run,".txt"))
deseqres <- transform(deseqres, zscore=logFC/SE,sig=ifelse(padj<0.1,"sig","not_sig"))
deseqres <- deseqres[order(deseqres$padj,-abs(deseqres$logFC)),]
top20genes <- head(deseqres,n=20)
dge_top <- deseqres[deseqres$identifier %in% top20genes$identifier, ]

top10m <- reshape2::dcast(identifier ~ cluster, data=dge_top, value.var="zscore",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$identifier
myMat <- top10m[,-1]
top10p <- reshape2::dcast(identifier ~ cluster, data=dge_top, value.var="padj",fun.aggregate=mean,na.rm=T)
rownames(top10p) <- top10p$identifier
pMat <- top10p[,-1]

myMat[is.na(pMat)] <- NA
pMat[is.na(myMat)] <- NA

png(width = 17, height = 8, file=paste0(figuredir,project,".",resset,".",dimset,".",var,"-",i,".",run,".zscore_heatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2, insig = "pch",pch=4,p.mat = pMat, sig.level=0.1,digits = 1) +
scale_fill_gradient2(low = "blue", high =  "red", mid = "white", midpoint = 0)
dev.off()

# Function to get asteriks
labs.function = function(x){
  case_when(x >= 0.1 ~ "",
            x < 0.1 & x >= 0.05 ~ "*",
            x < 0.05 & x >= 0.01 ~ "**",
            x < 0.01 ~ "***")
}

# Get asteriks matrix based on p-values
p.labs = pMat  %>%                      
  mutate_all(labs.function)

# Reshaping asteriks matrix to match ggcorrplot data output
p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

# Initial ggcorrplot
cor.plot = ggcorrplot(myMat, lab = TRUE, lab_size=2,digits = 1)+scale_fill_gradient2(low = "blue", high =  "red", mid = "white", midpoint = 0)

# Subsetting asteriks matrix to only those rows within ggcorrplot data
p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2), 
                                  paste0(cor.plot[["data"]]$Var1, cor.plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

# Add asteriks to ggcorrplot
cor.plot.labs = cor.plot + 
  geom_text(aes(x = p.labs$Var1, 
                y = p.labs$Var2), 
            label = p.labs$lab, 
            nudge_y = 0.25, 
            size = 5)

png(width = 17, height = 8, file=paste0(figuredir,project,".",resset,".",dimset,".",var,"-",i,".",run,".zscore_heatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
print(cor.plot.labs)
dev.off()

#, p.mat = pMat, sig.level=0.1
p <- ggplot(deseqres, aes(x=cluster, y=zscore)) +
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_hline(yintercept = 0)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 12, height = 12, file=paste0(figuredir,project,".",resset,".",dimset,".",var,"-",i,".",run,".zscore_scatter.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

####################################### not updated
##################################################

for (i in unique(dds$treats)){
var="sex"
dge <- fread(paste0(outFolder,paste0(outFolder,"deseqres_",var,"-",i,".txt")))










dge <- dge[order(dge$padj,-abs(dge$logFC)),]

dge %>% 
    group_by(Cluster) %>%
    head(n=40) -> top10
top10 <- as.data.frame(top10)
