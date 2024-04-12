maleGenes <- data.frame(ENSG=c('ENSG00000129824', 'ENSG00000198692', 'ENSG00000067048', 'ENSG00000012817'),
        symbol=c("RPS4Y1","EIF1AY","DDX3Y","KDM5D"))
femaleGenes <- data.frame(ENSG=c('ENSG00000229807'),symbol=c("XIST"))
males <- subset(x = sc, subset = RPS4Y1 > 2 | EIF1AY>2 | DDX3Y>2 | KDM5D>2)
female <- subset(x = sc, subset = XIST > 2 )

sc@meta.data$sex_inf <- ifelse(sc@meta.data$Sample_ID %in% males@meta.data$Sample_ID,"Male",
    ifelse(sc@meta.data$Sample_ID %in% female@meta.data$Sample_ID, "Female","unknown"))

table(sc@meta.data$sex_inf,useNA = "always")
length(which(!sc@meta.data$sex_inf==sc@meta.data$sex_alph))
head(sc@meta.data[which(!sc@meta.data$sex_inf==sc@meta.data$sex_alph),])

        var="age"
        cat("running deseq ",var,cluster," \n")
        #cluster_metadata_t <- transform(cluster_metadata_t, age_cscale=scale(age, scale = FALSE))
        #cluster_metadata_t_dt <- as.data.table(cluster_metadata_t)[ , agequartile := cut(cluster_metadata_t$age,
        #  breaks = quantile(cluster_metadata_t$age, probs = 0:4/4), labels = 1:4, right = T,include.lowest = TRUE)]
        #cluster_metadata_t_dt <- transform(cluster_metadata_t_dt, agequartile=relevel(as.factor(agequartile), ref = "2"))
        #cluster_metadata_t_dt <- as.data.frame(cluster_metadata_t_dt)
        #rownames(cluster_metadata_t_dt) <- rownames(cluster_metadata_t)
        #cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_t_dt))]
        #cluster_counts_t <- cluster_counts[,which(colnames(cluster_counts) %in% rownames(cluster_metadata_t))]
        #all(colnames(cluster_counts_t) == rownames(cluster_metadata_t_dt))
        #res <- lapply(c("1","3","4"),function(t){
        #  res <- results(dds, contrast=c("agequartile",t,"2"))
        #  sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
        #  names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC')
        #  sub.table <- sub.table[!is.na(sub.table$padj), ]
        #  cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
        #  sub.table 
        #}); names(res) <- c("1","3","4")
        #res_n <- ldply(res,data.frame)

        res_n$var=var
        res_n$cluster=cluster
        res_n$treats =i
        fwrite(res_n, sep='\t', quote=F, row.names=F, col.names=F, paste0(outFolder,"deseqres_",var,"-",i,".txt"),append=T)
        sigDEGs <- subset(res_n,padj<0.1)
        sigDEGs <- sigDEGs[order(sigDEGs$padj), ]
        fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=F, paste0(outFolder,"sigDEGs_",var,"-",i,".txt"),append=T)
        table <- data.frame(symb=var, variable= var, cluster=cluster)
        table$CTRL_numb <- paste(nrow(cluster_metadata))
        table$CTRL_uniq_ID <- paste(length(unique(cluster_metadata$Sample_ID)))
        table$CTRL_gene_numb <- paste(nrow(cluster_counts))
        table$CTRL_DEGs <- paste(nrow(sigDEGs))
        fwrite(table, sep='\t', quote=F, row.names=F, col.names=F, paste0(outFolder,"stats_all_cell_types-",var,"-",i,".txt"),append=T)    

################################
################################ compare cluster run -- was not turning up any enrichments and abandonded

### GO
library("AnnotationHub")
library(tidyverse)

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

### background gene list
BgDf <-geneIDs[geneIDs$symbol %in% dge$gene_symbol,] %>% dplyr::filter(grepl("protein_coding", biotype))

### geneCluster for enrichment results
res <- dge%>%drop_na(logFC,padj)%>%dplyr::filter( padj<0.1) #abs(logFC)>0.5,
geneCluster <- res %>%mutate(direction=ifelse(logFC>0, "Up", "Down"))

eg=bitr(geneCluster$gene_symbol, "SYMBOL", "ENTREZID", "org.Hs.eg.db")[, "ENTREZID"]


cg <- compareCluster(gene_symbol~Cluster+direction,
                     data=geneCluster,  
                     universe=unique(BgDf$symbol),
                     fun="enrichGO", 
                     OrgDb="org.Hs.eg.db",
                     pvalueCutoff=1,
                     qvalueCutoff=1, 
                     ont="ALL",
                     minGSSize=0,
                     maxGSSize=1000)

cg <- cg%>%mutate(Cluster1=Cluster)
opfn1 <- paste(outdir, prefix, "_enrichGO.rds", sep="")
write_rds(cg, opfn1)
                          
fig1 <- enrichplot::dotplot(diff_analysis, x=~Cluster + direction, showCategory=5)+
        scale_x_discrete(labels=lab2)+
        theme(axis.text.x=element_text(angle=60, hjust=1, size=12),
              axis.text.y=element_text(size=10))
        
figfn <- paste(outdir, prefix, "_Figure2.1_GO.png", sep="")
png(figfn, width=3000, height=2000, res=180)
print(fig1)
dev.off()