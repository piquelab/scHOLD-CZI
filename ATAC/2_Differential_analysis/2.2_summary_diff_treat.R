####
####
library(tidyverse)
library(Matrix)
## library(DESeq2)
## library(biobroom)
library(data.table)

##
library(cowplot)
library(RColorBrewer)
library(scales)
library(viridis)
library(circlize)
library(ComplexHeatmap)
library(openxlsx)
library(ggrastr)

###

rm(list=ls())

option <- "option_nFeature15K_cluster_res0.12"
outdir2 <- paste("./2_motif.outs/summary_", option, "_treat/", sep="")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)
 



###
###
outdir <- paste("./2_motif.outs/", option, "_th20_treat/", sep="")
fn_ls <- list.files(outdir, "_motifs.results.txt.gz$")
 
summ <- map_dfr(fn_ls, function(ii){
    ###
    cluster0 <- gsub("_motifs.*", "", ii)
    cat(ii, cluster0, "\n")
    
    fn0 <- paste(outdir, ii, sep="")
    res0 <- fread(fn0, header=T, data.table = FALSE)
    th0 <- quantile(abs(res0$beta), 0.9)
    
    summ0 <- data.frame(cluster=cluster0, nsig_fdr0.01=sum(res0$padj<0.01), nsig_fdr0.05=sum(res0$padj<0.05), 
           nsig_fdr0.1=sum(res0$padj<0.1), nsig_th0.1=sum(abs(res0$beta)>th0),
           nsig_fdr0.1_th0.1 = sum(res0$padj<0.1 & abs(res0$beta)>th0),
           th0.1=th0, nind=res0$nind[1])
    summ0
})

summ2 <- summ%>%mutate(cl_val=as.numeric(gsub("^C", "", cluster)))%>%arrange(cl_val)
                       
###
### 
opfn <- paste(outdir2, "1_summary.sigs.xlsx", sep="")
write.xlsx(summ2, file=opfn)


###
### response DAMs
res_sig <- map_dfr(fn_ls, function(ii){
    ###
    cluster0 <- gsub("_motifs.*", "", ii)
    
    fn0 <- paste(outdir, ii, sep="")
    res0 <- fread(fn0, header=T, data.table = FALSE)
    th0 <- quantile(abs(res0$beta), 0.9)

    ### save
    res2 <- res0%>%filter(padj<0.1, abs(beta)>th0)

    cat(ii, cluster0, th0, nrow(res2), "\n")    
    res2
})

### 143 unique DAMs
opfn <- paste(outdir2, "2_response_DAMs.comb.txt.gz", sep="")
fwrite(res_sig, file=opfn, sep="\t", quote=F, na=NA)
 




##########################
### QQ plots 
########################

 
###
### combine results
res_all <- NULL
for ( ii in fn_ls){
    ###
    fn0 <- paste(outdir, ii, sep="")
    res <- fread(fn0, header=T, data.table=F)
    ntest <- nrow(res)
    
    res <- res%>%arrange(pval)%>%
       mutate(observed=-log10(pval), expected=-log10(ppoints(ntest)))
       ###
       ### 
    res_all <- rbind(res_all, res)    
}

opfn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
fwrite(res_all, file=opfn, sep="\t", quote=F, na=NA)


###
### plots
col2 <- hue_pal()(11) ### default ggplot color
names(col2) <- paste("C", 0:10, sep="")
 
infn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
plotDF <- fread(infn, header=T, data.table=F)

cluster_sel <- paste("C", 0:7, sep="")
plotDF2 <- plotDF%>%filter(Cluster%in%cluster_sel)

## motifs <- plotDF2%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()

###
p0 <- ggplot(plotDF2, aes(x=expected, y=observed, color=Cluster))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=col2[1:8],
        guide=guide_legend(override.aes=list(size=3)))+
    ## facet_wrap(~psycho_variable, scales="free_y", ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          axis.title=element_text(size=12),
          axis.text=element_text(size=10))

          ## strip.text=element_text(size=12))

###
###
figfn <- paste(outdir2, "Figure2.1_treat.qq.pdf", sep="")
ggsave(figfn, p0, width=4.5, height=4)


## x <- plotDF%>%filter(Cluster=="C0", psycho_variable==varSel[4])


 







