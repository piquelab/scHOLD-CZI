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
outdir2 <- paste("./2_motif.outs/summary_", option, "_psycho_LPS/", sep="")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)




###
###
var_type <- "psycho"
outdir <- paste("./2_motif.outs/", option, "_th20_psycho/", sep="")
fn_ls <- list.files(outdir, "_motifs.summ.tsv$")
fn_ls <- fn_ls[grepl("LPS", fn_ls)]
fn_ls <- fn_ls[grepl(var_type, fn_ls)]

summ <- map_dfr(fn_ls, function(ii){
    ###
    cat(ii, "\n") 
    fn0 <- paste(outdir, ii, sep="")
    df0 <- read_tsv(fn0, show_col_types=F)
    df0
})



###
### save
opfn <- paste(outdir2, "1_", var_type, "_summary.sigs.tsv", sep="")
write_tsv(summ, file=opfn)

## ###
## ### show wide tables 
## fn <- paste(outdir2, "1_psychosocial_top10_summary.sigs.tsv", sep="")
## summ_df <- read_tsv(fn, show_col_types=F)%>%as.data.frame()


## varSel <- sort(unique(read.table("./psycho_var_dir/psychosocial_top10.txt")$V1))
## varSel <- varSel[c(3, 4, 7:10)]


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
 
## ## varSel <- summ2$variable[c(3,4,6:10)]

## ###
## ### save
## opfn <- paste(outdir2, "1_psychosocial_summary_mat.compare.xlsx", sep="")
## write.xlsx(summ2%>%filter(variable%in%varSel), file=opfn)
 


## ###
## ### summarize sigs
## summ2 <- summ_df%>%filter(variable%in%varSel)%>%arrange(variable, .locale="en")%>%
##     pivot_wider(id_cols=variable, names_from="Cluster", values_from=nsig_fdr0.1_t, values_fill=0)
## ###
## opfn2 <- paste(outdir2, "1.2_psychosocial_mat.sigs.xlsx", sep="")
## write.xlsx(summ2, file=opfn2, overwrite=T)
    

###
### nind and motifs

## summ2 <- summ_df%>%filter(variable%in%varSel)%>%
##     dplyr::select(variable, Cluster, nind, nmotif=ngene_test, nsig=nsig_fdr0.1_t)%>%
##     arrange(variable, .locale="en")
## summ2 <- summ2%>%
##     pivot_wider(id_cols=variable, names_from="Cluster",
##         values_from=c(nind, nsig), values_fill=0,
##         names_glue = "{Cluster}_{.value}")%>%as.data.frame()

## colSel <- c("variable", sort(colnames(summ2)[-1]))
## summ2 <- summ2[, colSel]

## opfn2 <- paste(outdir2, "1.2_psychosocial_mat.xlsx", sep="")
## write.xlsx(summ2, file=opfn2, overwrite=T)



###
### full
var_type <- "psycho"
fn0 <- paste(outdir2, "1_", var_type, "_summary.sigs.tsv", sep="")
summ_df <- read_tsv(fn0, show_col_types=F)

var_all <- sort(unique(summ_df$variable))
varSel <- var_all[c(3, 4, 7:10)]

###
### summarize sigs
summ2 <- summ_df%>%filter(variable%in%varSel)%>%arrange(variable, .locale="en")%>%
    pivot_wider(id_cols=variable, names_from="Cluster", values_from=nsig_fdr0.1_t, values_fill=0)
###
opfn2 <- paste(outdir2, "1.2_", var_type, "_mat.sigs.xlsx", sep="")
write.xlsx(summ2, file=opfn2, overwrite=T)



###
### full table 

summ2 <- summ_df%>%filter(variable%in%varSel)%>%
    dplyr::select(variable, Cluster, nind, nmotif=ngene_test, nsig=nsig_fdr0.1_t)%>%
    arrange(variable, .locale="en")

summ2 <- summ2%>%
    pivot_wider(id_cols=variable, names_from="Cluster",
        values_from=c(nind, nmotif, nsig), values_fill=0,
        names_glue = "{Cluster}_{.value}")%>%as.data.frame()

colSel <- c("variable", sort(colnames(summ2)[-1]))
summ2 <- summ2[, colSel]
 
opfn2 <- paste(outdir2, "1.2_", var_type, "_mat.full.xlsx", sep="")
write.xlsx(summ2, file=opfn2, overwrite=T)







##########################
### QQ plots 
########################

 
outdir2 <- paste("./2_motif.outs/summary_", option, "_psycho_LPS/", sep="")
var_type <- "psycho"

###
outdir <- paste("./2_motif.outs/", option, "_th20_psycho/", sep="")
fn_ls <- list.files(outdir, "results.txt.gz$")
fn_ls <- fn_ls[grepl("LPS", fn_ls)]
fn_ls <- fn_ls[grepl(var_type, fn_ls)]


### variable
var2 <- sort(unique(read.table("./psycho_var_dir/psychosocial_top10.txt")$V1))
varSel <- var2[c(3, 4, 7:10)]



###
### combine results
res_all <- NULL
for ( ii in fn_ls){
    ###
    fn0 <- paste(outdir, ii, sep="")
    res0 <- fread(fn0, header=T, data.table=F)

    for ( var0 in varSel){

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

opfn <- paste(outdir2, "2_", var_type, "_plotData.comb.txt.gz", sep="")
fwrite(res_all, file=opfn, sep="\t", quote=F, na=NA)


###
### plots
col2 <- hue_pal()(11) ### default ggplot color
names(col2) <- paste("C", 0:10, sep="")

infn <- paste(outdir2, "2_", var_type, "_plotData.comb.txt.gz", sep="")
plotDF <- fread(infn, header=T, data.table=F)

cluster_sel <- paste("C", 0:7, sep="")
plotDF2 <- plotDF%>%filter(Cluster%in%cluster_sel)

## plotDF2 <- plotDF%>%filter(psycho_variable%in%varSel) ##, !Cluster%in%c("C9"))

## motifs <- plotDF2%>%filter(p.adjusted<0.1)%>%pull(gene)%>%unique()

###
p0 <- ggplot(plotDF2, aes(x=expected, y=observed, color=Cluster))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=col2[1:8],
        guide=guide_legend(override.aes=list(size=3)))+
    facet_wrap(~psycho_variable, scales="free_y", ncol=3)+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          axis.title=element_text(size=12),
          axis.text=element_text(size=10),
          strip.text=element_text(size=12))

###
###
figfn <- paste(outdir2, "Figure2.1_", var_type, "_LPS.qq.pdf", sep="")
ggsave(figfn, p0, width=8, height=5)

## x <- plotDF%>%filter(Cluster=="C0", psycho_variable==varSel[4])



#################################################
### compare isel (old) and  ISEL_mean (new)
##################################################

## infn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
## df0 <- fread(infn, header=T, data.table=F)

## ###
## plotDF <- map_dfr(sort(unique(df0$Cluster)), function(ii){
##     ###
##     df_old <- df0%>%filter(Cluster==ii, psycho_variable=="isel")%>%
##         dplyr::select(Cluster, gene, beta_x=estimate, zscore_x=statistic,  padj_x=padj_t)
##     df_new <- df0%>%filter(Cluster==ii, psycho_variable=="ISEL_Mean")%>%
##         dplyr::select(gene, beta_y=estimate, zscore_y=statistic, padj_y=padj_t)

##     df2 <- df_old%>%left_join(df_new, by="gene")
##     df2
## })    



## ###
## ###
## plotDF2 <- plotDF%>%
##     mutate(sig_gr=case_when(padj_x<0.1&padj_y>0.1 ~ "sig1",
##             padj_x>0.1&padj_y<0.1 ~ "sig2",
##             padj_x<0.1&padj_y<0.1 ~ "sig3",
##             TRUE ~ "sig0"))
 

## ###
## ### position 
## mypos <- function(x, a=0.1){
##     ##
##     R <- max(x)-min(x)
##     pos <- min(x)+a*R
##     pos
## }

## ###
## ### eq
## feq <- function(x){
##    ## 
##    r <- round(as.numeric(x$estimate),digits=3)
##    p <- x$p.value
##    if ( p<0.001) symb <- "***"
##    if ( p>=0.001 & p<0.01) symb <- "**"
##    if ( p>=0.01 & p<0.05) symb <- "*"
##    if ( p>0.05) symb <- "NS"

##    eq <- bquote(italic(R)==.(r)) ##~","~.(symb))
##    eq
## }

 
## ###
## ### annotation text 
## anno_df2 <- plotDF2%>%
##    group_by(Cluster)%>%
##    nest()%>%
##    mutate(corr=map(data, ~cor.test((.x)$beta_x, (.x)$beta_y, method="pearson")),
##           eq=map(corr,feq),
##           rr=map_dbl(corr, ~(.x)$estimate),
##           xpos=mypos(plotDF2$beta_x, 0.25),
##           ypos=mypos(plotDF2$beta_x, 0.9))%>%
##        dplyr::select(-data,-corr)

 
       

## ###
## ### plots 
## pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=beta_x, y=beta_y))+ 
##     geom_point(aes(color=sig_gr), size=0.8)+
##     scale_color_manual(values=c("sig0"="grey", "sig1"="red", "sig2"="blue", "sig3"="green"), 
##         labels=c("sig0" = "Not", "sig1"="DAM-old", "sig2"="DAM-new", "sig3"="Both"),
##         breaks=c("sig0", "sig1", "sig2", "sig3"),
##         guide=guide_legend(override.aes=list(size=2)))+
##     geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
##     geom_abline(slope=1, intercept=0, linewidth=0.3)+
##     facet_wrap(~Cluster, scales="fixed", nrow=3)+
##     scale_x_continuous(bquote(beta~"on motif in isel"), expand=expansion(mult=0.1))+
##     scale_y_continuous(bquote(beta~"on motif in ISEL_Mean"), expand=expansion(mult=0.1))+
##     theme_bw()+
##     theme(strip.text=element_text(size=12),
##           axis.title=element_text(size=10),
##           axis.text=element_text(size=10),
##           legend.title=element_blank(),
##           legend.key.size=grid::unit(0.4, "cm"))
 
## figfn <- paste(outdir2, "Figure0.1_beta_isel.scatter.png", sep="")
## ggsave(figfn, pcomb, width=900, height=800, units="px", dpi=120)


## ###
## ### compare z-score 
## ###


## ### annotation text 
## anno_df2 <- plotDF2%>%
##    group_by(Cluster)%>%
##    nest()%>%
##    mutate(corr=map(data, ~cor.test((.x)$zscore_x, (.x)$zscore_y, method="pearson")),
##           eq=map(corr,feq),
##           rr=map_dbl(corr, ~(.x)$estimate),
##           xpos=mypos(plotDF2$zscore_x, 0.25),
##           ypos=mypos(plotDF2$zscore_x, 0.9))%>%
##        dplyr::select(-data,-corr)
        

## ###
## ### plots 
## pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=zscore_x, y=zscore_y))+ 
##     geom_point(aes(color=sig_gr), size=0.8)+
##     scale_color_manual(values=c("sig0"="grey", "sig1"="red", "sig2"="blue", "sig3"="green"), 
##         labels=c("sig0" = "Not", "sig1"="DAM-old", "sig2"="DAM-new", "sig3"="Both"),
##         breaks=c("sig0", "sig1", "sig2", "sig3"),
##         guide=guide_legend(override.aes=list(size=2)))+
##     geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
##     geom_abline(slope=1, intercept=0, linewidth=0.3)+
##     facet_wrap(~Cluster, scales="fixed", nrow=3)+
##     scale_x_continuous(bquote(italic(Z)~"-score"~"on motif in isel"), expand=expansion(mult=0.1))+
##     scale_y_continuous(bquote(italic(Z)~"-score"~"on motif in ISEL_Mean"), expand=expansion(mult=0.1))+
##     theme_bw()+
##     theme(strip.text=element_text(size=12),
##           axis.title=element_text(size=10),
##           axis.text=element_text(size=10),
##           legend.title=element_blank(),
##           legend.key.size=grid::unit(0.4, "cm"))
 
## figfn <- paste(outdir2, "Figure0.1_zscore_isel.scatter.png", sep="")
## ggsave(figfn, pcomb, width=900, height=800, units="px", dpi=120)





#####################################
### compare CTRL and LPS 
#####################################

var2 <- sort(read.table("./psycho_var_dir/psychosocial_top10.txt")$V1)
varSel <- var2[c(3, 4, 7:10)]


var_type <- "psycho"

###
fn0 <- paste("./2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_",
    var_type, "_plotData.comb.txt.gz", sep="")

res0 <- fread(fn0, header=T, sep="\t", data.table=F)
res0 <- res0%>%mutate(comb=paste(Cluster, psycho_variable, sep="_"))

###
fn2 <- paste("./2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho_LPS/2_",
    var_type, "_plotData.comb.txt.gz", sep="")
res2 <- fread(fn2, header=T, sep="\t", data.table=F)
res2 <- res2%>%mutate(comb=paste(Cluster, psycho_variable, sep="_"))  ##%>%filter(!Cluster%in%c("C9"))

###
cl_var <- sort(unique(res2$comb))
plotDF <- map_dfr(cl_var, function(ii){
    ##
    x0 <- res0%>%filter(comb==ii)%>%
        dplyr::select(Cluster, psycho_variable, gene, beta_CTRL=estimate,
                      zscore_CTRL=statistic, fdr_CTRL=padj_t)
    ##
    x1 <- res2%>%filter(comb==ii)%>%
        dplyr::select(gene, beta_LPS=estimate,
                      zscore_LPS=statistic, fdr_LPS=padj_t)
    ##
    xx <- x0%>%inner_join(x1, by="gene")
    cat(ii, nrow(x0), nrow(x1), "\n")
    xx
})


plotDF2 <- plotDF%>%filter(psycho_variable%in%varSel)%>%
    mutate(sig_gr=case_when(fdr_CTRL<0.1&fdr_LPS>0.1 ~ "sig1",
            fdr_CTRL>0.1&fdr_LPS<0.1 ~ "sig2",
            fdr_CTRL<0.1&fdr_LPS<0.1 ~ "sig3",
            TRUE ~ "sig0"))

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
   group_by(Cluster, psycho_variable)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$beta_CTRL, (.x)$beta_LPS, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=mypos(plotDF2$beta_CTRL, 0.25),
          ypos=mypos(plotDF2$beta_LPS, 0.9))%>%
       dplyr::select(-data,-corr)

 
       

###
### plots 
pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=beta_CTRL, y=beta_LPS))+ 
    geom_point(aes(color=sig_gr), size=0.8)+
    scale_color_manual(values=c("sig0"="grey", "sig1"="red", "sig2"="blue", "sig3"="green"), 
        labels=c("sig0" = "Not", "sig1"="DAM-ctrl", "sig2"="DAM-LPS", "sig3"="Both"),
        breaks=c("sig0", "sig1", "sig2", "sig3"),
        guide=guide_legend(override.aes=list(size=2)))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    facet_grid(Cluster~psycho_variable, scales="fixed")+
    scale_x_continuous(bquote(beta~"on motif in CTRL"), expand=expansion(mult=0.1))+
    scale_y_continuous(bquote(beta~"on motif in LPS"), expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_blank(),
          legend.key.size=grid::unit(0.4, "cm"))
 
figfn <- paste(outdir2, "Figure3.1_LPSvsCTRL_", var_type, "_beta.scatter.png", sep="")
ggsave(figfn, pcomb, width=1200, height=1000, units="px", dpi=120)



################################
### scatter plots of zscore
#################################
###
### annotation text 
anno_df2 <- plotDF2%>%
   group_by(psycho_variable)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$zscore_CTRL, (.x)$zscore_LPS, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=mypos(plotDF2$zscore_CTRL, 0.25),
          ypos=mypos(plotDF2$zscore_LPS, 0.9))%>%
       dplyr::select(-data,-corr)

       

###
### plots 
pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=zscore_CTRL, y=zscore_LPS))+ 
    geom_point(aes(color=sig_gr), size=0.8)+
    scale_color_manual(values=c("sig0"="grey", "sig1"="red", "sig2"="blue", "sig3"="green"), 
        labels=c("sig0" = "Not", "sig1"="DAM-ctrl", "sig2"="DAM-LPS", "sig3"="Both"),
        breaks=c("sig0", "sig1", "sig2", "sig3"),
        guide=guide_legend(override.aes=list(size=2)))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    facet_grid(Cluster~psycho_variable, scales="fixed")+
    scale_x_continuous(bquote(italic(Z)~"-score on motif in CTRL"), expand=expansion(mult=0.1))+
    scale_y_continuous(bquote(italic(Z)~"-score on motif in LPS"), expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_blank(),
          legend.key.size=grid::unit(0.4, "cm"))
 
figfn <- paste(outdir2, "Figure3.1_LPSvsCTRL_", var_type, "_zscore.scatter.png", sep="")
ggsave(figfn, pcomb, width=1200, height=1000, units="px", dpi=120)





###################################################################################
#### compare LPS vs CTRL across clusters and split psychosocial variables
###################################################################################


###
### annotation text 
anno_df2 <- plotDF2%>%
   group_by(psycho_variable)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$beta_CTRL, (.x)$beta_LPS, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=map_dbl(data, ~mypos((.x)$beta_CTRL, 0.25)),
          ypos=map_dbl(data, ~mypos((.x)$beta_LPS, 0.9)))%>%
       dplyr::select(-data,-corr)
 

###
### plots 
pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=beta_CTRL, y=beta_LPS))+ 
    geom_point(aes(color=sig_gr), size=0.8)+
    scale_color_manual(values=c("sig0"="grey", "sig1"="red", "sig2"="blue", "sig3"="green"), 
        labels=c("sig0" = "Not", "sig1"="DAM-ctrl", "sig2"="DAM-LPS", "sig3"="Both"),
        breaks=c("sig0", "sig1", "sig2", "sig3"),
        guide=guide_legend(override.aes=list(size=2)))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    geom_vline(xintercept=0, linewidth=0.3)+
    geom_hline(yintercept=0, linewidth=0.3)+
    facet_wrap(~psycho_variable, ncol=3, scales="free")+
    scale_x_continuous(bquote(beta~"on motif in CTRL"), expand=expansion(mult=0.1))+
    scale_y_continuous(bquote(beta~"on motif in LPS"), expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          legend.key.size=grid::unit(0.4, "cm"))
   
figfn <- paste(outdir2, "Figure3.2_beta_psycho_cluster.scatter.png", sep="")
ggsave(figfn, pcomb, width=1000, height=600, units="px", dpi=120)



###
### scatter plots of z-score 

anno_df2 <- plotDF2%>%
   group_by(psycho_variable)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$zscore_CTRL, (.x)$zscore_LPS, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=map_dbl(data, ~mypos((.x)$zscore_CTRL, 0.25)),
          ypos=map_dbl(data, ~mypos((.x)$zscore_LPS, 0.9)))%>%
       dplyr::select(-data,-corr)
 

###
### plots 
pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=zscore_CTRL, y=zscore_LPS))+ 
    geom_point(aes(color=sig_gr), size=0.8)+
    scale_color_manual(values=c("sig0"="grey", "sig1"="red", "sig2"="blue", "sig3"="green"), 
        labels=c("sig0" = "Not", "sig1"="DAM-ctrl", "sig2"="DAM-LPS", "sig3"="Both"),
        breaks=c("sig0", "sig1", "sig2", "sig3"),
        guide=guide_legend(override.aes=list(size=2)))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    geom_vline(xintercept=0, linewidth=0.3)+
    geom_hline(yintercept=0, linewidth=0.3)+
    facet_wrap(~psycho_variable, ncol=3, scales="free")+
    scale_x_continuous(bquote(italic(Z)~"-score on motif in CTRL"), expand=expansion(mult=0.1))+
    scale_y_continuous(bquote(italic(Z)~"-score on motif in LPS"), expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          legend.key.size=grid::unit(0.4, "cm"))
   
figfn <- paste(outdir2, "Figure3.2_zscore_psycho_cluster.scatter.png", sep="")
ggsave(figfn, pcomb, width=1000, height=600, units="px", dpi=120)










###########################################################################
### compare cluster current and old different cellranger-output
#############################################################################


type <- "psycho"

cl_df <- data.frame(old=paste("C", 0:4, sep=""), new=paste("C", c(0:2, 4, 6), sep=""),
                    MCls=c("CD4+", "CD8+", "NK", "Mono", "Bcell"))
cl_df <- cl_df%>%mutate(comb=paste(gsub("C", "New", new), gsub("C", "Old", old), MCls, sep="_"))                    


###
### old-with old cellranger output
fn0 <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/old/analyses_new_2025_08_07/2_Differential_analysis/2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho_LPS/2_psycho_plotData.comb.txt.gz"
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
        dplyr::select(gene, beta_x=estimate,
                      zscore_x=statistic, fdr_x=padj_t)
    ##
    x1 <- res2%>%filter(comb==ii)%>%
        dplyr::select(comb, Cluster2, psycho_variable, gene, beta_y=estimate,
                      zscore_y=statistic, fdr_y=padj_t)
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



