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
type <- "psycho" ##"PFAS"
outdir2 <- paste("./2_motif.outs/summary_", option, "_", type, "_CTRL/", sep="")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)

### correct t test


###
outdir <- paste("./2_motif.outs/", option, "_th20_psycho/", sep="")
fn_ls <- list.files(outdir, "_motifs.summ.tsv$")
fn_ls <- fn_ls[grepl("CTRL", fn_ls)]
fn_ls <- fn_ls[grepl(type, fn_ls)]

summ <- map_dfr(fn_ls, function(ii){
    ###
    cat(ii, "\n") 
    fn0 <- paste(outdir, ii, sep="")
    df0 <- read_tsv(fn0, show_col_types=F)
    df0
})

###
###
opfn <- paste(outdir2, "1_", type, "_summary.sigs.tsv", sep="")
write_tsv(summ, file=opfn)


###
### show wide tables 
## fn <- paste(outdir2, "1_th", th0, "_summary.sigs.tsv", sep="")
## summ_df <- read_tsv(fn, show_col_types=F)%>%as.data.frame()


## ###
## summ2 <- summ_df%>%
##     dplyr::select(variable, Cluster, nind, nmotif=ngene_test, nsig_z=nsig_fdr0.1, nsig_t=nsig_fdr0.1_t)%>%
##     arrange(variable, .locale="en")
## summ2 <- summ2%>%
##     pivot_wider(id_cols=variable, names_from="Cluster",
##         values_from=c(nsig_z, nsig_t), values_fill=0,
##         names_glue = "{Cluster}_{.value}")%>%as.data.frame()

## colSel <- c("variable", sort(colnames(summ2)[-1]))
## summ2 <- summ2[, colSel]
 
## varSel <- summ2$variable[c(3,4,6:10)]

## ###
## ### save
## opfn <- paste(outdir2, "1_th20_summary_mat.compare.xlsx", sep="")
## write.xlsx(summ2%>%filter(variable%in%varSel), file=opfn)
 


## ###
## ### summarize sigs
## summ2 <- summ_df%>%filter(variable%in%varSel)%>%arrange(variable, .locale="en")%>%
##     pivot_wider(id_cols=variable, names_from="Cluster", values_from=nsig_fdr0.1_t, values_fill=0)
## ###
## opfn2 <- paste(outdir2, "1.2_th", th0, "_mat.sigs.xlsx", sep="")
## write.xlsx(summ2, file=opfn2, overwrite=T)
    

## ###
## ### nind and motifs

## summ2 <- summ_df%>%filter(variable%in%varSel)%>%
##     dplyr::select(variable, Cluster, nind, nmotif=ngene_test, nsig=nsig_fdr0.1_t)%>%
##     arrange(variable, .locale="en")
## summ2 <- summ2%>%
##     pivot_wider(id_cols=variable, names_from="Cluster",
##         values_from=c(nind, nsig), values_fill=0,
##         names_glue = "{Cluster}_{.value}")%>%as.data.frame()

## co lSel <- c("variable", sort(colnames(summ2)[-1]))
## summ2 <- summ2[, colSel]

## opfn2 <- paste(outdir2, "1.2_th", th0, "_mat.xlsx", sep="")
## w rite.xlsx(summ2, file=opfn2, overwrite=T)

 

####
### summary PFAS
##var0 <- "cytokines"
type <- "psycho"
fn0 <- paste(outdir2, "1_", type, "_summary.sigs.tsv", sep="")
summ_df <- read_tsv(fn0, show_col_types=F)

var_all <- sort(unique(summ_df$variable))
varSel <- var_all[c(3,4,7:10)] ##varSel ##var_all[c(3, 4, 7:10)]

summ2 <- summ_df%>%dplyr::filter(variable%in%varSel)%>%
    arrange(variable, .locale="en")%>%
    pivot_wider(id_cols=variable, names_from="Cluster", values_from=nsig_fdr0.1_t, values_fill=0)
###
opfn2 <- paste(outdir2, "1.2_", type, "_mat.sigs.xlsx", sep="")
write.xlsx(summ2, file=opfn2, overwrite=T)


### full

summ2 <- summ_df%>%dplyr::filter(variable%in%varSel)%>%
    dplyr::select(variable, Cluster, nind, nmotif=ngene_test, nsig=nsig_fdr0.1_t)%>%
    arrange(variable, .locale="en")

summ2 <- summ2%>%
    pivot_wider(id_cols=variable, names_from="Cluster",
        values_from=c(nind, nmotif, nsig), values_fill=0,
        names_glue = "{Cluster}_{.value}")%>%as.data.frame()

colSel <- c("variable", sort(colnames(summ2)[-1]))
summ2 <- summ2[, colSel]

opfn2 <- paste(outdir2, "1.2_", type, "_mat.full.xlsx", sep="")
write.xlsx(summ2, file=opfn2, overwrite=T)



## ##########################
## ### QQ plots 
## ########################

## fn <- "./psycho_var_dir/PFAS_vars.txt"
## var2 <- sort(read.table(fn, header=F)$V1)
## ## varSel <- sort(var2)[c(3, 4, 6:9)]


 
## outdir <- "./2_motif.outs/correct_excludeX/Cluster_res0.07_th20_psycho/"
## fn_ls <- list.files(outdir, ".results.txt.gz$")
## fn_ls <- fn_ls[grepl("CTRL", fn_ls)]

## fn_ls <- fn_ls[grepl("PFAS", fn_ls)]
   
## res_all <- NULL
## for ( ii in fn_ls){
##     ###

##     fn0 <- paste(outdir, ii, sep="")
##     res0 <- fread(fn0, header=T, data.table=F)

##     for ( var0 in var2){

##        res <- res0%>%filter(psycho_variable==var0) 
##        ntest <- nrow(res)

##        cat(var0, res$nind[1], ntest, "\n")
        
##        res <- res%>%
##            arrange(pval_t)%>%
##            mutate(observed=-log10(pval_t), expected=-log10(ppoints(ntest)))
##        ###
##        ### 
##        res_all <- rbind(res_all, res)
##    }     
## }
 
## opfn <- paste(outdir2, "2_PFAS_plotData.comb.txt.gz", sep="")
## fwrite(res_all, file=opfn, sep="\t", quote=F, na=NA)


## ###
## ### plots
## ## outdir2 <- "./2_motif.outs/results_old_includeX/summary_Cluster_res0.07_psycho/"
## col2 <- hue_pal()(9) ### default ggplot color

## infn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
## plotDF <- fread(infn, header=T, data.table=F)

## var2 <- sort(unique(plotDF$psycho_variable))
## varSel2 <- var2[c(3, 4, 7:10)]
## plotDF2 <- plotDF%>%filter(psycho_variable%in%varSel2)

## ## motifs <- plotDF2%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()

## ###
## p0 <- ggplot(plotDF2, aes(x=expected, y=observed, color=Cluster))+
##     rasterise(geom_point(size=0.06), dpi=300)+
##     geom_abline(color="grey")+
##     scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
##         "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C9"="#FF61C3"),
##         guide=guide_legend(override.aes=list(size=3)))+
##     facet_wrap(~psycho_variable, scales="free_y", ncol=3)+
##     xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
##     ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           legend.text=element_text(size=9),
##           legend.key.size=grid::unit(0.4, "cm"),
##           axis.title=element_text(size=12),
##           axis.text=element_text(size=10),
##           strip.text=element_text(size=12))

## ###
## ###
## figfn <- paste(outdir2, "Figure2.1_th20_psycho_ctrl.qq.png", sep="")
## ggsave(figfn, p0, width=800, height=500, units="px", dpi=120)





###############################################
### qq plots for psycho, cytokines, PFAS 
#################################################



fn <- "./psycho_var_dir/psychosocial_top10.txt"
###fn <- "./psycho_var_dir/cytokines_var11.txt"
###fn <- "./psycho_var_dir/PFAS_vars.txt"
var2 <- sort(read.table(fn, header=F)$V1)
var2 <- var2[c(3, 4, 7:10)]
 
####
option <- "option_nFeature15K_cluster_res0.12"
type <- "psycho"
outdir <- paste("./2_motif.outs/", option, "_th20_psycho/", sep="")
fn_ls <- list.files(outdir, ".results.txt.gz$")
fn_ls <- fn_ls[grepl("CTRL", fn_ls)]
fn_ls <- fn_ls[grepl(type, fn_ls)]

   
res_all <- NULL
for ( ii in fn_ls){
    ###

    fn0 <- paste(outdir, ii, sep="")
    res0 <- fread(fn0, header=T, data.table=F)

    for ( var0 in var2){

       res <- res0%>%filter(psycho_variable==var0) 
       ntest <- nrow(res)

       cat(var0, res$nind[1], ntest, "\n")
        
       res <- res%>%
           arrange(pval_t)%>%
           mutate(observed=-log10(pval_t), expected=-log10(ppoints(ntest)))
       ###
       ### 
       res_all <- rbind(res_all, res)
   }     
}
 
opfn <- paste(outdir2, "2_", type, "_plotData.comb.txt.gz", sep="")
fwrite(res_all, file=opfn, sep="\t", quote=F, na=NA)


###
### plots
## outdir2 <- "./2_motif.outs/results_old_includeX/summary_Cluster_res0.07_psycho/"
col2 <- hue_pal()(11) ### default ggplot color
names(col2) <- paste("C", 0:10, sep="")

infn <- paste(outdir2, "2_", type, "_plotData.comb.txt.gz", sep="")
plotDF <- fread(infn, header=T, data.table=F)

cluster_sel <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
plotDF2 <- plotDF%>%filter(Cluster%in%cluster_sel)

###
###
p0 <- ggplot(plotDF2, aes(x=expected, y=observed, color=Cluster))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=col2[1:8],
        guide=guide_legend(override.aes=list(size=2)))+
    facet_wrap(~psycho_variable, scales="free_y", ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          ## legend.position="inside",
          ## legend.position.inside=c(0.9, 0.12),
          axis.title=element_text(size=12),
          axis.text=element_text(size=10),
          strip.text=element_text(size=10))

###
###
figfn <- paste(outdir2, "Figure2.1_", type, "_ctrl.qq.pdf", sep="")
ggsave(figfn, p0, width=8, height=5)
##ggsave(figfn, p0, width=900, height=600, units="px", dpi=100)
### ggsave(figfn, p0, width=800, height=600, units="px", dpi=100)
### ggsave(figfn, p0, width=800, height=500, units="px", dpi=120)  









#################################################################################################
### compare motif activity new (correct cellranger output) and old ( old cellranger output) 
##################################################################################################

type <- "psycho"


cl_df <- data.frame(old=paste("C", 0:4, sep=""), new=paste("C", c(0:2, 4, 6), sep=""),
                    MCls=c("CD4+", "CD8+", "NK", "Mono", "Bcell"))
cl_df <- cl_df%>%mutate(comb=paste(gsub("C", "New", new), gsub("C", "Old", old), MCls, sep="_"))                    


###
### old-with old cellranger output
fn0 <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/old/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_psycho_plotData.comb.txt.gz"
###
res0 <- fread(fn0, header=T, sep="\t", data.table=F)
res0 <- res0%>%filter(Cluster%in%cl_df$old)%>%
    left_join(cl_df%>%dplyr::select(Cluster=old, Cluster2=comb), by="Cluster")
res0 <- res0%>%mutate(comb=paste(Cluster2, psycho_variable, sep="_"))

###
### new with up-to-date cellranger output
fn2 <- paste(outdir2, "2_", type, "_plotData.comb.txt.gz", sep="")
res2 <- fread(fn2, header=T, sep="\t", data.table=F)
res2 <- res2%>%filter(Cluster%in%cl_df$new)%>%
    left_join(cl_df%>%dplyr::select(Cluster=new, Cluster2=comb), by="Cluster")
res2 <- res2%>%mutate(comb=paste(Cluster2, psycho_variable, sep="_"))
 
###
cl_var <- sort(unique(res2$comb))
plotDF <- map_dfr(cl_var, function(ii){
    ##
    x0 <- res0%>%filter(comb==ii)%>%
        dplyr::select(gene, beta_x=estimate, se_x=stderror, 
                      zscore_x=statistic, pt_x=pval_t, fdr_x=padj_t)
    ##
    x1 <- res2%>%filter(comb==ii)%>%
        dplyr::select(comb, Cluster2, psycho_variable, gene, beta_y=estimate, se_y=stderror,
                      zscore_y=statistic, pt_y=pval_t, fdr_y=padj_t)
    ##
    xx <- x0%>%inner_join(x1, by="gene")
    cat(ii, nrow(x0), nrow(x1), "\n")
    xx
})

varSel <- c("ISEL_Mean", "PSS_all_mean")
plotDF2 <- plotDF%>%filter(psycho_variable%in%varSel)%>%  
    mutate(sig_gr=case_when(fdr_x<0.1&fdr_y>0.1 ~ "sig1",
            fdr_x>0.1&fdr_y<0.1 ~ "sig2",
            fdr_x<0.1&fdr_y<0.1 ~ "sig3",
            TRUE ~ "sig0"))

## x <- plotDF2%>%filter(comb=="New4_Old3_Mono_ISEL_Mean", sig_gr!="sig0")%>%arrange(desc(zscore_x))

## x2 <- x%>%select(gene, beta_x, se_x, zscore_x, pt_x, fdr_x,
##                  beta_y, se_y, zscore_y, pt_y, fdr_y)



###
### position 
mypos <- function(x, a=0.1){
    ##
    R <- max(x)-min(x)
    pos <- min(x)+a*R
    pos
}

###
### eq
feq <- function(x){
   ## 
   r <- round(as.numeric(x$estimate),digits=3)
   p <- x$p.value
   if ( p<0.001) symb <- "***"
   if ( p>=0.001 & p<0.01) symb <- "**"
   if ( p>=0.01 & p<0.05) symb <- "*"
   if ( p>0.05) symb <- "NS"

   eq <- bquote(italic(R)==.(r)) ##~","~.(symb))
   eq
}


###
### annotation text 
anno_df2 <- plotDF2%>%
   group_by(Cluster2, psycho_variable)%>%
   ## group_by(comb)%>%   
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$beta_x, (.x)$beta_y, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=mypos(plotDF2$beta_x, 0.25),
          ypos=mypos(plotDF2$beta_y, 0.9))%>%
       dplyr::select(-data,-corr)

 
       

###
### plots 
pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=beta_x, y=beta_y))+ 
    geom_point(aes(color=sig_gr), size=0.8)+
    scale_color_manual(values=c("sig0"="grey", "sig1"="blue", "sig2"="red", "sig3"="green"), 
        labels=c("sig0" = "Not", "sig1"="DAM(old)", "sig2"="DAM(new)", "sig3"="Both"),
        breaks=c("sig0", "sig1", "sig2", "sig3"),
        guide=guide_legend(override.aes=list(size=2)))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    facet_grid(psycho_variable~Cluster2, scales="fixed")+
    ## facet_wrap(~comb, nrow=2, ncol=5, dir="v", scales="free")+
    scale_x_continuous(bquote(beta~"on motif (old)"), expand=expansion(mult=0.1))+
    scale_y_continuous(bquote(beta~"on motif (new)"), expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=10),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_blank(),
          legend.key.size=grid::unit(0.4, "cm"))
 
figfn <- paste(outdir2, "Figure3.1_", type, "_beta.scatter.pdf", sep="")
## ggsave(figfn, pcomb, width=1200, height=1000, units="px", dpi=120)
ggsave(figfn, pcomb, width=12, height=5)


################################
### scatter plots of zscore
#################################
###
### annotation text  
anno_df2 <- plotDF2%>%
   group_by(Cluster2, psycho_variable)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$zscore_x, (.x)$zscore_y, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=mypos(plotDF2$zscore_x, 0.25),
          ypos=mypos(plotDF2$zscore_y, 0.9))%>%
       dplyr::select(-data,-corr)

       

###
### plots  
pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=zscore_x, y=zscore_y))+ 
    geom_point(aes(color=sig_gr), size=0.8)+
    scale_color_manual(values=c("sig0"="grey", "sig1"="blue", "sig2"="red", "sig3"="green"), 
        labels=c("sig0" = "Not", "sig1"="DAM (old)", "sig2"="DAM(new)", "sig3"="Both"),
        breaks=c("sig0", "sig1", "sig2", "sig3"),
        guide=guide_legend(override.aes=list(size=2)))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    facet_grid(psycho_variable~Cluster2, scales="fixed")+
    scale_x_continuous(bquote(italic(Z)~"-score on motif (old)"), expand=expansion(mult=0.1))+
    scale_y_continuous(bquote(italic(Z)~"-score on motif (new)"), expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=10),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_blank(),
          legend.key.size=grid::unit(0.4, "cm"))
 
figfn <- paste(outdir2, "Figure3.1_", type, "_zscore.scatter.pdf", sep="")
##ggsave(figfn, pcomb, width=1200, height=1000, units="px", dpi=120)
ggsave(figfn, pcomb, width=12, height=5)



####
####

############################################################
### get differential gene expression results from Ali and Cindy 
############################################################

###
### Before we have replicates results for each var::celltype because the file ALL.0.15.13.deseqres has include all cluster results. Now I corrected it.

option <- "option_nFeature15K_cluster_res0.12"
var_type <- "psycho"
outdir2 <- paste("./2_motif.outs/summary_", option, "_psycho_CTRL/", sep="")
varSel <- sort(unique(read.table("./psycho_var_dir/psychosocial_top10.txt")$V1))
varSel2 <- varSel[c(4, 7, 10)]
 
 
###
### get differential gene results from Ali DESeq2 results 
dir_deg <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/deseqres/"

fns <- list.files(dir_deg, "^ALL.0.1.13")
df_fn <- data.frame(fn_name=fns,
     psycho_var=gsub(".*.deseqres_|-RNA-.*", "", fns),
     treat=gsub(".*-RNA-|\\..*", "", fns))

## df_fn%>%filter(treat=="CTRL", psycho_var%in%varSel[1], grepl("_DESeq.txt$", fn_name))

df_fn2 <- df_fn%>%
    dplyr::filter(treat=="CTRL", psycho_var%in%varSel2, grepl("^ALL.0.1.13.deseqres", fn_name),
                  grepl("_DESeq.txt$", fn_name))

 
res <- map_dfr(df_fn2$fn_name, function(ii){
    ##
    fn0 <- paste(dir_deg, ii, sep="")
    x <- read.table(fn0, header=T, sep="\t")
    x
})    
 
opfn <- paste(outdir2, "Ali_RNA_", var_type, ".diff.DESeq.txt.gz", sep="")
fwrite(res, file=opfn, sep="\t", quote=F, na=NA)

## ##
DEGs <- res%>%dplyr::filter(padj<0.1)%>%pull(identifier)%>%unique()  ## 2756 DEGs
opfn2 <- paste(outdir2, "Ali_RNA_", var_type, ".DEGs.txt", sep="")
write.table(DEGs, file=opfn2, quote=F, row.names=F, col.names=F)



################################
##### DEG list and DAM list 
#################################


outdir2 <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"

varSel <- sort(unique(read.table("./psycho_var_dir/psychosocial_top10.txt")$V1))
varSel <- varSel[c(7, 10)]

###
### RNA results 
fn <- paste(outdir2, "Ali_RNA_psycho.diff.DESeq.txt.gz", sep="")
res_rna <- fread(fn, header=T, data.table=F)
 
cl_rna <- paste("C", 0:4, sep="")
summ_rna <- res_rna%>%filter(cluster%in%cl_rna, var%in%varSel)%>%
    group_by(cluster, var)%>%
    summarize(nDEG=sum(padj<0.1, na.rm=T), .groups="drop")%>%ungroup()
summ_rna <- summ_rna%>%mutate(rna_cluster=gsub("^C", "R", cluster))%>%
    dplyr::select(rna_cluster, var, nDEG)%>%
    mutate(comb=paste(rna_cluster, var, sep="_"))


###
### DAMs results
fn <- paste(outdir2, "2_psycho_plotData.comb.txt.gz", sep="")
res_motif <- fread(fn, header=T, data.table=F)
cl_atac <- paste("C", 0:7, sep="")
summ_motif <- res_motif%>%filter(Cluster%in%cl_atac, psycho_variable%in%varSel)%>%
    group_by(Cluster, psycho_variable)%>%
    summarize(nTF=sum(padj_t<0.1, na.rm=T), .groups="drop")%>%ungroup()


###
### cluster infor
fn <- "../1.2_ArchR_process/2_Integrate_output/harmony_default/TableS_summ.tsv"
cluster_df <- read_tsv(fn)
cl_rna <- cluster_df$RNA_cluster
names(cl_rna) <- cluster_df$ATAC_cluster

MCl_name <- cluster_df$MCl_comb
names(MCl_name) <- cluster_df$ATAC_cluster

      
summ_motif <- summ_motif%>%
    mutate(atac_cluster=gsub("C", "A", Cluster), rna_cluster=cl_rna[atac_cluster],
           comb=paste(rna_cluster, psycho_variable, sep="_"))
summ_motif <- summ_motif%>%dplyr::select(comb, atac_cluster, psycho_variable, nTF)

###
### combine

summ2 <- summ_rna%>%full_join(summ_motif, by="comb")


### add cell type
summ2 <- summ2%>%mutate(MCls=MCl_name[atac_cluster])

summ3 <- summ2%>%dplyr::select(MCls, rna_cluster, atac_cluster, psycho_variable, nDEG, nTF)

opfn <- paste(outdir2, "list_cluster_var_nDEGs_nTF_new.txt", sep="")
write.table(summ3, file=opfn, row.names=F, col.names=T, sep="\t", quote=F)
 
opfn2 <- paste(outdir2, "list_cluster_var_nDEGs_nTF_new.xlsx", sep="")
write.xlsx(summ3, file=opfn2)
###
## opfn <- paste(outdir2, "list"


##############
### DEGs
##############

varSel <- c("ISEL_Mean", "PSS_all_mean")

###
outdir2 <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"
fn <- paste(outdir2, "Ali_RNA_psycho.diff.DESeq.txt.gz", sep="")
res <- fread(fn, header=T, data.table=F)
res <- res%>%filter(var%in%varSel)

summ <- res%>%group_by(cluster, var)%>%
    summarize(ngene=n(), nsig=sum(padj<0.1, na.rm=T), .groups="drop")%>%ungroup()
summ <- summ%>%mutate(comb=paste(cluster, var, sep="_"))
 
###
### old
fn2 <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/old/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/Ali_RNA_psycho.diff.DESeq.txt.gz" 
old <- fread(fn2, header=T, data.table=F)
old <- old%>%filter(var%in%varSel)

summ_old <- old%>%group_by(cluster, var)%>%
    summarize(ngene=n(), nsig=sum(padj<0.1, na.rm=T), .groups="drop")%>%ungroup()
summ_old <- summ_old%>%mutate(comb=paste(cluster, var, sep="_"))%>%
    dplyr::select(comb, ngene_old=ngene, nsig_old=nsig)

summ2 <- summ%>%left_join(summ_old, by="comb")
opfn <- paste(outdir2, "0_summ_statistics.xlsx", sep="")
write.xlsx(summ2, file=opfn)

###
### save xlsx file 
## fn <- paste(outdir2, "list_cluster_var_nDEGs20_nTF1_new.txt", sep="")
## x <- read.table(fn, sep="\t", header=T)

## ###
## ### save
## x2 <- x%>%dplyr::select(rna_cluster, atac_cluster, celltype=MCls, var, comb, nDEG, nTF)
## opfn <- paste(outdir2, "TableS_res_summ.xlsx", sep="")
## write.xlsx(x2, file=opfn)

###
### END








##############################
### stratified QQ plots 
###############################

## ##
## fn_motif <- "../1.2_ArchR_process/4_motif_output/Motif_jaspar2022.match.mat.rds"
## motif_match <- read_rds(fn_motif)

## df_row <- rowRanges(motif_match)
## anno <- data.frame(chr=


## mat2 <- assay(motif_match)


## fn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
## res <- fread(fn, header=T)





## option <- "Cluster_res0.07"
## outdir2 <- paste("./2_motif.outs/summary_", option, "_psycho/", sep="")
## if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


## ##
## prefix <- "../1.2_ArchR_process/4_motif_output/motif_anno/"
## fn0 <- paste(prefix, "1_motif.list.txt", sep="")
## motifList <- read.table(fn0)$V1


## ###
## df2 <- map_dfr(motifList, function(ii){
##     ##
##     fn <- paste(prefix, "2_motif_", ii, ".annotation.txt", sep="")
##     x <- read.table(fn, header=T, sep="\t")
##     ##
##     df0 <- data.frame(motif_name = ii, ngene = length(unique(x$gene)), npeak = length(unique(x$peak)))
##     df0
## })


## ###
## ###

## p1 <- ggplot(df2, aes(x=ngene))+
##     geom_histogram(fill="white", color="grey50")+
##     xlab("#target genes") + ylab("#Motifs")+
##     ggtitle("#target genes per motif")+
##     theme_bw()+
##     theme(plot.title=element_text(hjust=0.5, size=10))

## ##
## p2 <- ggplot(df2, aes(x = npeak))+
##     geom_histogram(fill = "white", color = "grey50")+
##     xlab("#peaks") + ylab("#Motifs")+
##     ggtitle("#peaks per motif")+
##     theme_bw()+
##     theme(plot.title = element_text(hjust=0.5, size=10))


## ##
## figfn <- paste(outdir2, "Figure3.0_motifs.hist.png", sep="")
## pcomb <- plot_grid(p1, p2, nrow=2)
## ggsave(figfn, pcomb, width=380, height=600, units="px", dpi=120)

## ### save data
## data_fn <- paste(outdir2, "Figure3.0_plotdata.txt", sep="")
## write.table(df2, data_fn, quote=F, row.names=F, sep="\t")








