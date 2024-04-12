library(plyr);library(dplyr)
library(msigdbr)
library(clusterProfiler)
  require(splitstackshape)
  library(gridExtra)
library(condformat)
library(data.table)
library(ggplot2)
library(cowplot)

m_df <- bind_rows(msigdbr(species = "Homo sapiens", category = "H"),
                  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS"))
msigdbr_t2g = m_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_HOLD01-HOLD14_dbgap.ID_cziexp_03_04_2024.txt","CZ1_group.txt") #for testing
base <- args[1]
outFolder=paste0(base,"enrichr/")
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)

figuredir=paste0(outFolder,"figures/")
if (!file.exists(figuredir)) dir.create(figuredir, showWarnings=F)


var="LPSvsCTRL"
dge <- fread(paste0(base,"pseudobulk_ctrl/deseqres_",var,".txt"))
colnames(dge) <- c('gene_symbol', 'padj', 'pvalue', 'logFC','var','Cluster')

ddf <- split(dge,dge$Cluster)

list_em <- lapply(ddf,function(f){
#f <- ddf[[1]]
  n <- unique(as.character(f$Cluster))
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
em <- enricher(u$gene_symbol, TERM2GENE=msigdbr_t2g,pvalueCutoff =0.1)
  if(!is.null(em)){
res <- dplyr::filter(em@result, p.adjust < 0.05)
res <- res[order(res$p.adjust),]
   fwrite(res, sep='\t', quote=F, row.names=F, file=paste0(outFolder,n,"_",nn,"_gse_enrichr_fullres.txt"))
  if(!dim(res)[1]==0){
   s <- cSplit(as.data.table(res)[c(1:10),], "geneID", "/",direction = "long")
   s_dge <- merge(s, dge[dge$Cluster==n,c("gene_symbol","logFC","padj")], by.x="geneID",by.y="gene_symbol")
   s_dge <- s_dge[order(s_dge$padj,-abs(s_dge$logFC)),]
   names(s_dge)[c(10:11)] <- c("dge_log2FoldChange","dge_padj")
   #s <- transform(s, Z=qnorm(pvalue)) #remember direction is nonsensical
   fwrite(s_dge, sep='\t', quote=F, row.names=F, file=paste0(outFolder,n,"_",nn,"_gse_enrichr.txt"))
ss <- ldply(lapply(split(s_dge,s_dge$Description),function(i){
      ds <- paste(i$geneID[c(1:5)],collapse=';')
      df <- cbind(i[1,c("Description","qvalue","pvalue","p.adjust")],data.frame(geneID=ds))
      return(df)
      }),data.frame)[,-1]
   fwrite(ss, sep='\t', quote=F, row.names=F, file=paste0(outFolder,n,"_",nn,"_gse_enrichrtop.txt"))
  if(nn=="upregulated"){ 
      hc <- "indianred1"
      lc <- "darkred"
      title <- "Upregulated DEGs"
  } else {
      hc <- "lightblue"
      lc <- "darkblue"
    title <- "Downregulated DEGs"
  }
  em <- dplyr::mutate(em, direction=nn)
 #png(width = 8, height = 10, file=paste0("figure/",n,"_",nn,"_dotplot_enrichr.png"), pointsize=12, 
 #       bg = "transparent", units = "in", res = 1200)
  #p <-dotplot(em, showCategory=5) + scale_color_gradient(low = lc, high = hc, na.value = "darkgrey") + ggtitle(title)
 # print(p)
 # dev.off()
}}}
if(exists("em")){
  return(em)
}
  })
#fig <-plot_grid(plotlist=me, nrow=1, ncol=2, align="h")
#figfn <- paste(figuredir, n,".dotplot.png", sep="")
#png(figfn, width=3000, height=2000, res=180)
#print(fig)
#dev.off()
})

combined_gsea<-merge_result(unlist(list_em,recursive=F)) #must be a single level named list
figfn <- paste(figuredir,"me.dotplot.png", sep="")
png(figfn, width=3000, height=2000, res=180)
dotplot(combined_gsea, showCategory=5)
dev.off()

myDir <- outFolder #directory to load from
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
rm(list=ls(pattern="C[0-9]"))
diff_analysis <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
diff_analysis <- transform(diff_analysis, Cluster=gsub('_[^_]*$',"",analysis),direction=gsub(".*_","",analysis),p.adjust=round(p.adjust,3))
diff_analysis <- subset(diff_analysis, !Cluster=="")
diff_analysis <- diff_analysis[,c("Cluster","Description","geneID","p.adjust","direction")]

lapply(split(diff_analysis,diff_analysis$Cluster), function(i){
  #i <- split(diff_analysis,diff_analysis$Cluster)[[1]]
  n <- unique(i$Cluster)
ss <- condformat(i) %>%
  rule_text_color(geneID, ifelse(direction=="downregulated","darkblue", ifelse(direction=="upregulated", "red4", "black"))) %>%
  condformat2grob()
png(width = 16, height = 8, file=paste0(figuredir,"enrichrtop_table.",n,".png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
grid.arrange(ss)
dev.off()
})

