library(cowplot)
library(EnhancedVolcano)
library(ggcorrplot)
library(data.table)
library(plyr)
library(parallel)
library(DESeq2)


args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/") #for testing
base <- args[1]
outFolder=paste0(base,"pseudobulk_ctrl/")
figuredir=paste0(outFolder,"figures/")

filenames <- list.files(outFolder) #file list from directory
filenames1 <- filenames[grep("DESeq_output-", filenames)] #pick specific files from list
x <- na.omit(sapply(strsplit(filenames1,"(?<=.)(?=C[0-9])",perl=TRUE),function(y)y[2]))
clusters <- unique(sapply(strsplit(x,"[.]"),function(z)z[1]))

for (cluster in clusters){
    #cluster=clusters[1]
    cat("running ", cluster, "\n")
    var="LPSvsCTRL"
    load(paste0(outFolder,"DESeq_output-",var,cluster,".RData"))

    vsd <- vst(dds, blind=FALSE)
    # Plot PCA
    fig0 <- DESeq2::plotPCA(vsd, intgroup = c("treats"))+ggtitle("treats")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig1 <- DESeq2::plotPCA(vsd, intgroup = c("BATCH"))+ggtitle("BATCH")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig2 <- DESeq2::plotPCA(vsd, intgroup = c("sex_alph"))+ggtitle("sex")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig3 <- DESeq2::plotPCA(vsd, intgroup = c("age"))+ggtitle("age")+theme(plot.margin=margin(b=-0.8,unit="cm"))
    fig <-plot_grid(fig0, fig1, fig2,fig3, nrow=4, ncol=1, align="h")
    png(paste0(figuredir,"deseqPCA-",var,"-",cluster,".png"), width=1000, height=1000, res=120)
    print(fig)
    dev.off()

    sub.table <- fread(paste0(outFolder,"deseqres_",var,".txt"))
    colnames(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','var','Cluster')
    sub.table <- subset(sub.table, Cluster==cluster)
    sub.table <- sub.table[order(sub.table$padj,-abs(sub.table$logFC)), ]

    topp1 <- min(subset(sub.table, padj > quantile(padj, prob = 1 -99/100,na.rm=T))$pvalue)
    tolab <- c(unique(head(sub.table,n=20)$identifier))
    png(width = 8, height = 8, file=paste0(figuredir,"dge_volcano-",var,"-",cluster,".png"), pointsize=12, 
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
    for (i in unique(dds$treats)){
        # Transform counts for data visualization
        mclapply(c("sex","age","SES","DSES_07","DSES_09","HVS_mean","EDS_mean"),function(var){
        #mclapply(c("SES","DSES_07","DSES_09","HVS_mean","EDS_mean"),function(var){
        #var="sex"
        #load(paste0(outFolder,"DESeq_output-",i,"-",var,cluster,".RData"))
        cat("running ",var,i)
        sub.table <- fread(paste0(outFolder,"deseqres_",var,"-",i,".txt"))
            colnames(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','var','Cluster','treats')
        sub.table <- subset(sub.table, Cluster==cluster & treats==i)
        sub.table <- sub.table[order(sub.table$padj,-abs(sub.table$logFC)), ]

        topp1 <- min(subset(sub.table, padj > quantile(padj, prob = 1 -99/100,na.rm=T))$pvalue)
        tolab <- c(unique(head(sub.table,n=20)$identifier))
        png(width = 8, height = 8, file=paste0(figuredir,"dge_volcano-",var,"-",cluster,"-",i,".png"), pointsize=12, 
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
