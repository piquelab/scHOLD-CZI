corr <- cor(TMP0)[Neworder, Neworder]
mycol <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100) 
#mycol <- viridisLite::viridis(100)

x <- str_split(colnames(corr), "_", simplify=T)
tmp_column <- data.frame(celltype=x[,1], treatment=x[,2])
rownames(tmp_column) <- colnames(corr)
tmp_colors <- list(celltype=col2, treatment=col1)

fig2 <- pheatmap(corr, col=mycol, scale="none", border_color="NA",
                 cluster_rows=F, cluster_cols=F,
                 annotation_col=tmp_column, annotation_colors=tmp_colors, annotation_legend =T,
                 show_colnames=T, show_rownames=F, na_col="white")

figfn <- "./6_DEG.CelltypeNew_output/Filter2/Figure3.2_corr.beta.png"
png(figfn, width=1000, height=1000,res=150)
print(fig2)
dev.off()

