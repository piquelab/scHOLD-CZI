# Enrichr
library(ggcorrplot)
library(data.table)
library(plyr);library(dplyr)
library(msigdbr)
  require(splitstackshape)
library(clusterProfiler)
m_df <- bind_rows(msigdbr(species = "Homo sapiens", category = "H"),
                  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS"))
msigdbr_t2g = m_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
run="treatment_withCOMBAT"
dge <- fread(paste0(baseoutFolder,run,"/deseqres/ALL.0.2.11.deseqres.RNA-LPS_vs_RNA-CTRL.treatment_withCOMBAT.txt"))
figuredir <- paste0(baseoutFolder,run,"/figures/")

ddf <- split(dge,dge$cluster)
ddf10 <- lapply(ddf, function(i) {
  i <- i[order(i$padj,-abs(i$logFC)),]
  return(head(i,n=5))
  })
dge <- dge[order(dge$padj,-abs(dge$logFC)),]

dge %>% 
    group_by(cluster) %>%
    head(n=40) -> top10
top10 <- as.data.frame(top10)

top10 <- ldply(ddf10,data.frame)
top10m <- reshape2::dcast(identifier ~ cluster, data=top10, value.var="logFC",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$identifier
myMat <- top10m[,-1]

top10m <- reshape2::dcast(identifier ~ cluster, data=top10, value.var="padj",fun.aggregate=mean,na.rm=T)
rownames(top10m) <- top10m$gene_symbol
pMat <- top10m[,-1]

myMat[is.na(pMat)] <- NA
pMat[is.na(myMat)] <- NA

png(width = 17, height = 8, file=paste0(figuredir,"top_dge_heatmap.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
ggcorrplot(myMat, method = "square", outline.col = "grey", ggtheme = ggplot2::theme_bw(), lab = TRUE, lab_size=2, insig = "blank",p.mat = pMat, sig.level=0.1,digits = 1) +
scale_fill_gradient2(low = "blue", high =  "red", mid = "white", midpoint = 0)
    #theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = a))
dev.off()

lapply(ddf,function(f){
#f <- ddf[[1]]
  n <- unique(as.character(f$cluster))
message(paste0("running ",n))
f_df <- subset(f,padj<0.1) #& abs(logFC)>0.5)
f_df_up <- subset(f_df, logFC>=0)
f_df_up <- transform(f_df_up, direction="up")
f_df_down <- subset(f_df, logFC<0)
f_df_down <- transform(f_df_down, direction="down")
me <- lapply(list(down=f_df_down,up=f_df_up),function(u){
  #u <-f_df_down
  nn <- ifelse(mean(u$logFC,na.rm=T)<0,"downregulated","upregulated")
  if(!dim(u)[1]==0){
em <- enricher(u$identifier, TERM2GENE=msigdbr_t2g,pvalueCutoff =0.1)
  if(!is.null(em)){
res <- filter(em@result, p.adjust < 0.05)
res <- res[order(res$p.adjust),]
   fwrite(res, sep='\t', quote=F, row.names=F, file=paste0(baseoutFolder,run,"/enrichr/",n,"_",nn,"_gse_enrichr_fullres.txt"))
  if(!dim(res)[1]==0){
   s <- cSplit(as.data.table(res)[c(1:10),], "geneID", "/",direction = "long")
   s_dge <- merge(s, dge[dge$cluster==n,c("identifier","logFC","padj")], by.x="geneID",by.y="identifier")
   s_dge <- s_dge[order(s_dge$padj,-abs(s_dge$logFC)),]
   names(s_dge)[c(10:11)] <- c("dge_log2FoldChange","dge_padj")
   #s <- transform(s, Z=qnorm(pvalue)) #remember direction is nonsensical
   fwrite(s_dge, sep='\t', quote=F, row.names=F, file=paste0(baseoutFolder,run,"/enrichr/",n,"_",nn,"_gse_enrichr.txt"))
ss <- ldply(lapply(split(s_dge,s_dge$Description),function(i){
      ds <- paste(i$geneID[c(1:5)],collapse=';')
      df <- cbind(i[1,c("Description","qvalue","pvalue","p.adjust")],data.frame(geneID=ds))
      return(df)
      }),data.frame)[,-1]
   fwrite(ss, sep='\t', quote=F, row.names=F, file=paste0(baseoutFolder,run,"/enrichr/",n,"_",nn,"_gse_enrichrtop.txt"))
  if(nn=="upregulated"){ 
      hc <- "indianred1"
      lc <- "darkred"
      title <- "Upregulated DEGs"
  } else {
      hc <- "lightblue"
      lc <- "darkblue"
    title <- "Downregulated DEGs"
  }
  em <- dplyr::mutate(em, log10.p.adjust = -log10(as.numeric(p.adjust)))
 #png(width = 8, height = 10, file=paste0("figure/",n,"_",nn,"_barplot_enrichr.png"), pointsize=12, 
  #      bg = "transparent", units = "in", res = 1200)
  #p <-barplot(em, x="log10.p.adjust", showCategory=10) + scale_color_gradient(low = lc, high = hc, na.value = "darkgrey") + ggtitle(title)
  #print(p)
  #dev.off()
 #png(width = 8, height = 10, file=paste0("figure/",n,"_",nn,"_dotplot_enrichr.png"), pointsize=12, 
 #       bg = "transparent", units = "in", res = 1200)
  p <-dotplot(em, showCategory=10) + scale_color_gradient(low = lc, high = hc, na.value = "darkgrey") + ggtitle(title)
 # print(p)
 # dev.off()
}}}
if(exists("p")){
  return(p)
}
  })
#png(width = 13, height = 18, file=paste0("figure/",n,"_barplot_enrichr.png"), pointsize=12, 
#        bg = "transparent", units = "in", res = 1200)
#  ppp <- cowplot::plot_grid(plotlist =me, align = "v", nrow=3,ncol=2)
#  print(ppp)
#  dev.off()
png(width = 12, height = 9, file=paste0(baseoutFolder,run,"/enrichr/",n,"_dotplot_enrichr.png"), pointsize=12, 
        bg = "transparent", units = "in", res = 1200)
  ppp <- cowplot::plot_grid(plotlist =me, align = "v", nrow=1,ncol=2)
  print(ppp)
  dev.off()
})

myDir <- getwd() #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames1 <- filenames[grep("upregulated_gse_enrichrtop.txt", filenames)] #pick specific files from list
filenames2 <- filenames[grep("downregulated_gse_enrichrtop.txt", filenames)] #pick specific files from list
filenames <- c(filenames1,filenames2)
data_names <- gsub("_gse_enrichrtop.txt", "", filenames) #remove file ending
shortnames <- data_names
for(i in 1:length(filenames)) assign(shortnames[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')[,analysis:=shortnames[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
#combining data
ddf <- lapply(shortnames, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(shortnames)
rm(list=ls(pattern="vs"))
diff_analysis <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
diff_analysis <- transform(diff_analysis, contrast=gsub('_[^_]*$',"",analysis),direction=gsub(".*_","",analysis),p.adjust=round(p.adjust,3))
diff_analysis <- subset(diff_analysis, !contrast=="")
diff_analysis <- diff_analysis[,c("contrast","Description","geneID","p.adjust","direction")]
library(gridExtra)
library(condformat)

lapply(split(diff_analysis,diff_analysis$contrast), function(i){
  #i <- split(diff_analysis,diff_analysis$contrast)[[1]]
  n <- unique(i$contrast)
ss <- condformat(i) %>%
  rule_text_color(geneID, ifelse(direction=="downregulated","darkblue", ifelse(direction=="upregulated", "red4", "black"))) %>%
  condformat2grob()
png(width = 16, height = 8, file=paste0("figure/enrichrtop_table.",n,".png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
grid.arrange(ss)
dev.off()
})
