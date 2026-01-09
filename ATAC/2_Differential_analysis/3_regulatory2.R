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

library(dendsort, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")

## install.packages("dendsort", lib="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")

###
###
rm(list=ls())


###############################
### use cell type peaks 
##################################

outdir <- "./3_regulatory.outs/" 
if ( !file.exists(outdir)) dir.create(outdir, showWarnings = F, recursive = T)

indir <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/"

###
### input-1, differential gene results
fn_rna <- paste(indir, "Ali_RNA_psycho.diff.DESeq.txt.gz", sep="")
res_rna <- fread(fn_rna, header=T, data.table=F) 

## x <- res_rna%>%filter(cluster=="C1", var=="PSS_all_mean")
## fn <- paste(outdir2, "Ali_RNA.DEGs.txt", sep="") 
## DEGs <- read.table(fn, header=F)$V1

###
### input-2, differential TF motifs results   
fn_diff_motif <- paste(indir, "2_psycho_plotData.comb.txt.gz", sep="")
res_motif <- fread(fn_diff_motif, header=T, data.table=F)


###
### input-3, summary cluster from Ali
### cp /rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/DEGs_TF/list_var_clust_nDEGs20_nTF1.txt
###varSel <- c("ISEL_Mean", "PSS_all_mean")
fn_cl <- "./3_regulatory.outs/list_cluster_var_nDEGs_nTF_new.txt"
df_cl <- read.table(fn_cl, header=T, sep="\t")
df_cl <- df_cl%>%filter(nDEG>0, nTF>0)





## fn2 <- "./3_regulatory.outs/summary_nDEGs_nDAMs_per_cluster_noDEX.txt"
## df_ali <- read.table(fn2, header=T)
## df_ali <- df_ali%>%drop_na(nTF)



###
### input-4, motif annotation in genes_
dir_annot <- "../1.2_ArchR_process/4b_motif_output/motif_anno/"
## opfn_annot <- "../1.2_ArchR_process/4_motif_output/motif_anno/2_comb_motif.annot.txt.gz" 
## anno_motif <- fread(fn_annot, header=T, data.table=F)
fn_ls <- list.files(dir_annot, "2_comb_cluster")

anno_df <- data.frame(cluster=gsub("2_comb_cluster_|.motif.*", "", fn_ls), fn_name=fn_ls)



###
### input-5, motif list
fn_motif <- "./3_regulatory.outs/motif.list.txt"
motif_list <- unique(read.table(fn_motif)$V1)

###
### overall enrichment test 

################################################
### Test if DEGs are in enriched in DAM 
#################################################

###
###
df_all <- NULL
for ( i in 1:nrow(df_cl)){
   ##
   var0 <- df_cl$psycho_variable[i]
   cl_rna <- gsub("R", "C", df_cl$rna_cluster[i])
   cl_atac <- gsub("A", "C", df_cl$atac_cluster[i])

   fn0 <- anno_df%>%filter(cluster==cl_atac)%>%pull(fn_name) 
   anno_fn <- paste(dir_annot, fn0, sep="")
   anno_motif <- fread(anno_fn, header=T, sep="\t") 

   cat("rna:", cl_rna, "atac:", cl_atac, var0, "\n")  
 
   res0_rna <- res_rna%>%dplyr::filter(cluster==cl_rna, var==var0)
   res0_motif <- res_motif%>%dplyr::filter(Cluster==cl_atac, psycho_variable == var0)

   gene_test <- unique(res0_rna$identifier)
   ##intersect(unique(res0_rna$identifier), anno_motif$gene)

   DAM <- res0_motif%>%filter(padj_t<0.1)%>%pull(gene)%>%unique()
   gene_DAM <- anno_motif%>%dplyr::filter(motif_name%in%DAM)%>%pull(gene)%>%unique()
   gene_DAM <- intersect(gene_DAM, gene_test) 
    
   DEG <- res0_rna%>%dplyr::filter(padj<0.1)%>%pull(identifier)%>%unique()
   notDEG <- setdiff(gene_test, DEG) 

   ## 
   interest.in <- length(intersect(DEG, gene_DAM))
   interest.not <- length(setdiff(DEG, gene_DAM))
    
   not.interest.in <- length(intersect(notDEG, gene_DAM))
   not.interest.not <- length(setdiff(notDEG, gene_DAM))
    
   dmat <- matrix(c(interest.in, interest.not, not.interest.in, not.interest.not), 2, 2)
   enrich <- fisher.test(dmat)  ## two side test
   enrich2 <- fisher.test(dmat, alternative="greater") ## one side test

   df0 <- data.frame("interest.in" = interest.in, "interest.not" = interest.not,
       "not.interest.in" = not.interest.in, "not.interest.not" = not.interest.not,
       "odds" = enrich$estimate, "pval" = enrich2$p.value,
       "lower" = enrich$conf.int[1], "upper" = enrich$conf.int[2])

   df_all <- rbind(df_all, df0) 
   ### 
}
 
rownames(df_all) <- NULL

res_all <- cbind(df_cl, df_all)

res_all <- res_all%>%
    mutate(log_odds = log(odds), log_lower = log(lower), log_upper = log(upper),
           se = abs(log_odds - log_lower)/1.96)

####
#### 
opfn <- paste(outdir, "1_enrich_DAM.fisher.txt", sep="")
write.table(res_all, file = opfn, quote=F, row.names=F, col.names=T, sep="\t")



######################################
### compare old and new results 
########################################

outdir <- "./3_regulatory.outs/"
 
fn <- paste(outdir, "1_enrich_DAM.fisher.txt", sep="")
res <- read.table(fn, header = T, sep="\t")

###
###
fn <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/old/analyses_new_2025_08_07/2_Differential_analysis/3_regulatory.outs/correct_excludeX/1_enrich_DAM.fisher.txt"
old <- read.table(fn, sep="\t", header=T)

###
###
MCls_old <- c("A0_R0_T CD4+"="A0_R0_CD4+", "A1_R1_T CD8+"="A1_R1_CD8+", "A2_R2_NK"="A2_R2_NK",
              "A4_R3_Monocyte"="A3_R3_Monocyte")
MCl_sel <- names(MCls_old)

res2 <- res%>%filter(MCls%in%MCl_sel)%>%
    mutate(old_MCls=MCls_old[MCls], old_comb=paste(old_MCls, psycho_variable, sep="_"))%>%
    dplyr::select(MCls, rna_cluster, atac_cluster, psycho_variable, nDEG, nTF,
                  interest.in, interest.not, not.interest.in, not.interest.not,
                  odds, lower, upper, old_comb)

old2 <- old%>%filter(MCls%in%MCls_old)%>%
    mutate(comb=paste(MCls, var, sep="_"))%>%
    dplyr::select(comb, rna_cluster, atac_cluster, nDEG, nTF,
                  interest.in, interest.not, not.interest.in, not.interest.not, odds, lower, upper)
names(old2) <- paste("old_", names(old2), sep="")

dd <- res2%>%inner_join(old2, by="old_comb")%>%
    arrange(desc(psycho_variable), desc(MCls))

###
### save 
opfn <- paste(outdir, "2_compare_new_old.xlsx", sep="")
write.xlsx(dd, file=opfn)







#########################
### Forest plots
#########################
 
## outdir <- "./3_regulatory.outs/"
 
## fn <- paste(outdir, "1_enrich_DAM.fisher.txt", sep="")
## res <- read.table(fn, header = T, sep="\t")

## plotDF <- res%>%mutate(var_clust=paste(var, MCls, sep="_"),
##                        var_clust2 = fct_reorder(var_clust, nDEG))
 
## ## x <- plotDF%>%arrange(desc(nTFs), desc(var_clust))
## ## x <- x%>%dplyr::select(var_clust, varname, nDEGs, nTFs, interest.in, interest.not, not.interest.in, not.interest.not, odds, pval)

## ## ##
## ## opfn <- paste(outdir, "1.0_short.xlsx", sep="")
## ## write.xlsx(x, file=opfn)

## varSel <- sort(unique(plotDF$var))
## varSel2 <- varSel[c(2:3)]


## col2 <- hue_pal()(10)
## names(col2) <- paste("A", 0:9, sep="")

## ## col1 <- c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
## ##         "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7" = "#DB72FB")
##         ## guide=guide_legend(override.aes=list(size=2)))+
 
## plotDF2 <- plotDF%>%filter(var%in%varSel2)
## p0 <- ggplot(plotDF2, aes(x = log_odds, y = var_clust2, color = factor(atac_cluster)))+
##     geom_errorbarh(aes(xmax = log_upper, xmin = log_lower), linewidth = 0.5, height = 0.2)+
##     geom_point(shape = 19, size = 1.5)+
##     scale_colour_manual(values = col2,
##         labels = c("A0" = "A0_R0_CD4+", "A3" = "A3_R3_Mono", "A5" = "A5_R0_CD4"))+     
##     geom_vline(aes(xintercept = 0), linewidth = 0.25, linetype = "dashed")+
##     xlab("log odds ratio")+xlim(-2, 4)+
##     theme_bw()+
##     theme(axis.title.y = element_blank(),
##           axis.title.x = element_text(size=10),
##           axis.text.x = element_text(size=10),
##           axis.text.y = element_text(size=10),
##           legend.title = element_blank(),
##           legend.text = element_text(size=10),
##           legend.key.size = unit(0.4, "cm"))
 
## ###
## ###
## figfn <- paste(outdir, "Figure1_enrich_DAM.forest.png", sep="")
## ggsave(figfn, p0, width = 680, height = 480, units = "px", dpi=120)




## ###
## ### barplots of DAM

## p2 <- ggplot(plotDF2, aes(x=var_clust2, y=nTF, fill=factor(atac_cluster)))+
##     geom_bar(stat="identity")+
##     coord_flip()+
##     scale_fill_manual(values=col2, guide="none")+
##     scale_y_continuous("#DAMs")+
##     theme_bw()+
##     theme(legend.position = "none",
##           axis.title.y = element_blank(),
##           axis.text.y =  element_blank(),
##           axis.ticks.y = element_blank(),
##           #axis.title.x = element_text(size=12),
##           axis.title.x = element_text(size=12),
##           axis.text = element_text(size=10))

## ##
## figfn <- paste(outdir, "Figure1.2_DAM.bar.png", sep="")
## ggsave(figfn, p2, width = 280, height = 480, units = "px", dpi=120)
 

## ###
## ### barplots of DEGs 

## p3 <- ggplot(plotDF2, aes(x = var_clust2, y = nDEG, fill = factor(atac_cluster)))+
##     geom_bar(stat="identity")+
##     coord_flip()+
##     scale_fill_manual(values=col2, guide="none")+
##     scale_y_continuous("#DEGs")+
##     theme_bw()+
##     theme(legend.position = "none",
##           axis.title.y = element_blank(),
##           axis.text.y = element_blank(),
##           axis.ticks.y = element_blank(),
##           #axis.title.x = element_text(size=12),
##           axis.title.x = element_text(size=12),
##           axis.text = element_text(size=10))

## ##
## figfn <- paste(outdir, "Figure1.3_DEG.bar.png", sep="")
## ggsave(figfn, p3, width = 280, height = 480, units = "px", dpi=120)



####
#### histogram

## outdir <- "./3_regulatory.outs/Celltype/"

## dir_annot <- "../1.2_ArchR_process/4_motif_output/motif_anno/"
## ## opfn_annot <- "../1.2_ArchR_process/4_motif_output/motif_anno/2_comb_motif.annot.txt.gz" 
## ## anno_motif <- fread(fn_annot, header=T, data.table=F)
## fn_ls <- list.files(dir_annot, "2_cluster")

## anno_df <- data.frame(cluster=gsub("2_cluster_|_motif.*", "", fn_ls), fn_name=fn_ls)

## ncl <- nrow(anno_df)
## summ_df <- map_dfr(1:ncl, function(i){
##     ###
##     cl <- anno_df$cluster[i]
##     fn0 <- anno_df$fn_name[i]
##     fn_full <- paste(dir_annot, fn0, sep="")
##     x <- fread(fn_full, header=T, sep="\t")
##     cat(cl, fn0, "\n")
##     ##
##     df0 <- as.data.frame(table(x$motif_name))
##     df0$Cluster <- cl
##     df0
## })



## ###
## ###
 
## p0 <- ggplot(summ_df, aes(x=Freq))+
##    geom_histogram(fill="white", color="grey50")+
##    xlab("#genes by motif")+ylab("#Motifs")+ 
##    facet_wrap(~Cluster, scales="free", nrow=3)+
##    theme_bw()+
##    theme(axis.title=element_text(size=9),
##          axis.text=element_text(size=6))
         

## ### save
## figfn <- paste(outdir, "Figure0_motif.hist.png", sep="")
## ggsave(figfn, p0, width=600, height=600, units="px", dpi=120)


###
### 
## x <- read.xlsx("../1.2_ArchR_process/3_reCallPeaks_output/Peak_matrix_cluster/1_summ.xlsx")

## summ <- summ_df%>%group_by(Cluster)%>%summarize(nmotif_median=median(Freq), .groups="drop")
## summ <- x%>%left_join(summ, by="Cluster")
## opfn <- paste(outdir, "1.0_summ_num.xlsx", sep="")
## write.xlsx(summ, file=opfn)



###############################################
### 2, test which motifs are enriched in DEGs 
###############################################

### 3.0_submit_enrich.sh and 3.0_enrich_motif.R

outdir2 <- "./3_regulatory.outs/enrich_motif/" 
fn_ls <- list.files(outdir2, "enrich.motif.txt$")

varSel <- c("ISEL_Mean", "PSS_all_mean")

###
res_all <- map_dfr(fn_ls, function(ii){
    ##
    fn <- paste(outdir2, ii, sep="")
    res <- read.table(fn, header = T, sep = "\t")
    res <- res%>%
        mutate(comb = paste(cluster_rna, psycho_variable, sep="_"), padj = p.adjust(pval, "BH"))%>%
        arrange(pval)
    n0 <- nrow(res)
    res <- res%>%mutate(observed = -log10(pval), expected = -log10(ppoints(n0)))
    ###  
    res    
})


res_all <- res_all%>%filter(psycho_variable%in%varSel)
 

opfn <- paste(outdir2, "2_comb_enrich.motif.txt.gz", sep="")
fwrite(res_all, file=opfn, quote=F, na=NA, sep="\t")

## opfn <- paste(outdir, "2_comb_enrich.motif.txt.gz", sep="")
## fwrite(res, file=opfn, na=NA, sep="\t")

## ###
## ###
## fn_enrich <- paste(outdir, "2_comb_enrich.motif.txt.gz", sep="")
## res <- fread(fn_enrich, header=T, sep="\t")

## summ_df <- res%>%mutate(is_sig=padj<0.1)%>%
##     group_by(comb)%>%summarize(nsig=sum(is_sig), .groups="drop")

## ###
## fn <- "./3_regulatory.outs/list_var_clust_nDEGs20_nTF1.txt"
## summ <- read.table(fn, header=T, sep="\t")
## summ <- summ%>%mutate(comb=paste(variable, "_rna", Cluster_scRNA, sep=""))
## summ <- summ%>%left_join(summ_df, by="comb")

## ###
## opfn <- paste(outdir, "2_enrich.summ.tsv", sep="")
## write_tsv(summ, file=opfn)

## fn <- "./3_regulatory.outs/list_var_clust_nDEGs20_nTF1.txt"
## summ <- read.table(fn, header=T, sep="\t")




##################
### qq plots
##################

outdir2 <- "./3_regulatory.outs/enrich_motif/"

fn_enrich <- paste(outdir2, "2_comb_enrich.motif.txt.gz", sep="")
plotDF <- fread(fn_enrich, header=T, sep="\t")

var2 <- sort(unique(plotDF$psycho_variable))
varSel <- var2 ##var2[2:3]  ##var2[c(2, 3, 5:7)]
plotDF2 <- plotDF%>%filter(psycho_variable%in%varSel)%>%
    mutate(comb_new=paste(gsub("C", "R", cluster_rna),
                          gsub("C", "A", cluster_atac),  psycho_variable, sep="_"))

###
df_cl <- read_tsv("./3_regulatory.outs/TableS_nDEGs_nTFs_enriched.tsv")
combSel <- paste(df_cl$rna_cluster, df_cl$atac_cluster, df_cl$psycho_variable, sep="_")


col2 <- hue_pal()(7)
names(col2) <- paste("C", 0:6, sep="")

unique(plotDF2$comb_new)

plotDF2 <- plotDF2%>%filter(comb_new%in%combSel)

###
###  
p0 <- ggplot(plotDF2, aes(x = expected, y = observed, color=cluster_rna))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=col2,
         guide=guide_legend(override.aes=list(size=2)))+
    facet_wrap(~comb_new, scales="free_y", ncol=4)+
    ## facet_wrap(~comb, scales="free_y", ncol=8)+  ## ncol=7
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          legend.position = "inside",
          legend.position.inside=c(0.8, 0.25),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          strip.text=element_text(size=9))

 
###
### 
figfn <- paste(outdir2, "Figure1.0_enrich.qq.pdf", sep="")
ggsave(figfn, p0, width=9, height=4.8)



##############################
### compare old and new 
################################

outdir2 <- "./3_regulatory.outs/enrich_motif/"

fn_enrich <- "./3_regulatory.outs/enrich_motif/2_comb_enrich.motif.txt.gz"
res <- fread(fn_enrich, header=T, sep="\t", data.table=F)
res <- res%>%mutate(MCls=paste(gsub("C", "A", cluster_atac), gsub("C", "R", cluster_rna), sep="_"))

MCl_name <- c("A0_R0"="A0_R0_T CD4+", "A1_R1"="A1_R1_T CD8+", "A2_R2"="A2_R2_NK", "A4_R3"="A4_R3_Monocyte")
res <- res%>%mutate(MCls=MCl_name[MCls])
 
MCls_old <- c("A0_R0_T CD4+"="A0_R0_T CD4+", "A1_R1_T CD8+"="A1_R1_T CD8+", "A2_R2_NK"="A2_R2_NK",
              "A4_R3_Monocyte"="A3_R3_Monocyte")
MCl_sel <- names(MCls_old)

res <- res%>%mutate(old_MCls=MCls_old[MCls], comb=paste(old_MCls, psycho_variable, sep="_"))%>%
    filter(MCls%in%MCl_sel)
 

fn <- "/rs/rs_grp_scatac/schold/ATAC/sc-atac-cziHOLD/old/analyses_new_2025_08_07/2_Differential_analysis/3_regulatory.outs/correct_excludeX/enrich_motif/2_comb_enrich.motif.txt.gz"
old <- fread(fn, sep="\t", header=T, data.table=F)
old <- old%>%mutate(MCls=paste(gsub("C", "A", cluster_atac), gsub("C", "R", cluster_rna), sep="_"))
MCl_name <- c("A0_R0"="A0_R0_T CD4+", "A1_R1"="A1_R1_T CD8+", "A2_R2"="A2_R2_NK", "A3_R3"="A3_R3_Monocyte")
old <- old%>%mutate(MCls=MCl_name[MCls])%>%
    mutate(comb=paste(MCls, psycho_variable, sep="_"))

comb2 <- sort(unique(res$comb))

comb_sel <- unique(res$comb)[c(1, 2, 4, 7:8)]
###
###
plotDF <- map_dfr(comb_sel, function(ii){
    ###
    res0 <- res%>%filter(comb==ii)%>%mutate(log10p=-log10(pval))%>%
        dplyr::select(comb, MCls, old_MCls, psycho_variable, motif_name, log_odds_y=log_odds,
                      log10p_y=log10p)
    ###
    old0 <- old%>%filter(comb==ii)%>%mutate(log10p=-log10(pval))%>%
       dplyr::select(motif_name, log_odds_x=log_odds, log10p_x=log10p)
    ##
    res00 <- res0%>%left_join(old0, by="motif_name")
    res00
 })   



###
###
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

plotDF <- plotDF%>%filter(is.finite(log_odds_y), is.finite(log_odds_x), is.finite(log10p_y), is.finite(log10p_x))

###
### annotation text 
anno_df2 <- plotDF%>%
   group_by(comb)%>%
   ## group_by(comb)%>%   
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$log_odds_x, (.x)$log_odds_y, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=mypos(plotDF$log_odds_x, 0.25),
          ypos=mypos(plotDF$log_odds_y, 0.9))%>%
       dplyr::select(-data,-corr)

 
       

###
### plots
ncomb <- length(unique(plotDF$comb))
pcomb <- ggplot(plotDF, aes(x=log_odds_x, y=log_odds_y))+ 
    geom_point(color="grey", size=0.8)+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    facet_wrap(~comb, scales="fixed", nrow=2)+
    ## facet_wrap(~comb, nrow=2, ncol=5, dir="v", scales="free")+
    scale_x_continuous("log odds from old", expand=expansion(mult=0.1))+
    scale_y_continuous("log odds from new", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=10),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10))
   
 
figfn <- paste(outdir2, "Figure2.1_enrich_log_odds.scatter.pdf", sep="")
## ggsave(figfn, pcomb, width=1200, height=1000, units="px", dpi=120)
ggsave(figfn, pcomb, width=10, height=5)



###
###
p2 <- ggplot(plotDF, aes(x=log10p_x, y=log10p_y))+ 
    geom_point(color="grey", size=0.8)+
    ###geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    facet_wrap(~comb, scales="fixed", nrow=2)+
    ## facet_wrap(~comb, nrow=2, ncol=5, dir="v", scales="free")+
    scale_x_continuous("log10p from old", expand=expansion(mult=0.1))+
    scale_y_continuous("log10p from new", expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=10),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10))
   
 
figfn <- paste(outdir2, "Figure2.1_enrich_log10p.scatter.pdf", sep="")
## ggsave(figfn, pcomb, width=1200, height=1000, units="px", dpi=120)
ggsave(figfn, p2, width=10, height=5)






###################################
### ISEL_Mean and PSS_all_mean
####################################

## var2 <- sort(unique(plotDF$psycho_variable))
## varSel <- var2[c(2, 4)]
## var_val <- c("ISEL_Mean"=1, "PSS_all_mean"=2)
## plotDF2 <- plotDF%>%filter(psycho_variable%in%varSel)%>%
##     mutate(comb=paste(gsub("^C", "R", comb), gsub("^C", "A", cluster_atac), sep="_"), 
##            var_value = as.numeric(var_val[as.character(psycho_variable)]), 
##            comb2=fct_reorder(comb, var_value)) 
 
## col_cl <- c("C0"="#E60A0A",  "C1"="#F9A602", "C2"="#A65192", "C3"="#0570B0", "C4"="#00AF72", "C5"="#87CEFA")
## ###
## ###
## p0 <- ggplot(plotDF2, aes(x = expected, y = observed, color=cluster_rna))+
##     rasterise(geom_point(size=0.06), dpi=300)+
##     geom_abline(color="grey")+
##     scale_color_manual(values=col_cl,
##          guide=guide_legend(override.aes=list(size=2)))+
##     facet_wrap(~comb2, scales="free_y", ncol=4)+
##     ## facet_wrap(~comb, scales="free_y", ncol=8)+  ## ncol=7
##     xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
##     ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           legend.text=element_text(size=9),
##           legend.key.size=grid::unit(0.4, "cm"),
##           ## legend.position = "inside",
##           ## legend.position.inside=c(0.9, 0.25),
##           axis.title=element_text(size=10),
##           axis.text=element_text(size=10),
##           strip.text=element_text(size=9))

 
## ###
## ### 
## figfn <- paste(outdir2, "Figure1.1_enrich.qq.png", sep="")
## ggsave(figfn, p0, width=950, height=420, units="px", dpi=120)




####################
### summary sigs
####################

###
varSel <- c("ISEL_Mean", "PSS_all_mean")

##
fn0 <- "./3_regulatory.outs/list_cluster_var_nDEGs_nTF_new.txt" 
df_cl <- read.table(fn0, header=T, sep="\t")
df_cl <- df_cl%>%filter(nDEG>0, nTF>0, psycho_variable%in%varSel)


###
fn <- "./3_regulatory.outs/enrich_motif/2_comb_enrich.motif.txt.gz"
res_enrich <- fread(fn, header=T, data.table=F, sep="\t")

###
fn2 <- "./2_motif.outs/summary_option_nFeature15K_cluster_res0.12_psycho_CTRL/2_psycho_plotData.comb.txt.gz"
res <- fread(fn2, header=T, data.table=F, sep="\t")


for (i in 1:nrow(df_cl)){
    ##
    cl_atac <- gsub("A", "C", df_cl$atac_cluster[i])
    cl_rna <- gsub("R", "C", df_cl$rna_cluster[i])
    var0 <- df_cl$psycho_variable[i]
    ##
    DAM <- res%>%filter(padj_t<0.1, Cluster==cl_atac, psycho_variable==var0)%>%pull(gene)%>%unique()
    motif_enrich <- res_enrich%>%
        filter(cluster_atac==cl_atac,
               psycho_variable==var0, cluster_rna==cl_rna, odds>1, padj<0.1)%>%pull(motif_name)%>%unique()

    df_cl$nsig_enrich[i] <- length(motif_enrich)
    df_cl$nolap[i] <- length(intersect(DAM, motif_enrich))
}    


df_cl <- df_cl%>%
    dplyr::select(rna_cluster, atac_cluster, MCls, psycho_variable, nDEG, nTF, nsig_enrich, nolap)

df_cl <- df_cl%>%arrange(rna_cluster, psycho_variable)

opfn <- "./3_regulatory.outs/TableS_nDEGs_nTFs_enriched.tsv"
write_tsv(df_cl, file=opfn)

opfn2 <- "./3_regulatory.outs/TableS_nDEGs_nTFs_enriched.xlsx"
write.xlsx(df_cl, file=opfn2)




####################################
### forest plots 
####################################

## outdir <- "./3_regulatory.outs/Celltype/"


## varSel <- c("ISEL_Mean", "PSS_all_mean")
## var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Perceived Stress")

## ###
## ###
## fn <- "./3_regulatory.outs/list_cluster_var_nDEGs20_nTF1_new.txt"
## summ <- read.table(fn, header=T, sep="\t")
## summ <- summ%>%filter(var%in%varSel)%>%arrange(var)

## ##
## ### combine from RNA cluster and psycho_variable 
## fn_enrich <- "./3_regulatory.outs/Celltype/enrich_motif/2_comb_enrich.motif.txt.gz"
## res_enrich <- fread(fn_enrich, header=T, sep="\t")
## res_enrich <- res_enrich%>%filter(psycho_variable%in%varSel)


## ###
## ### diff motif, combine from ATAC cluster and psycho_variable 
## fn_motif <- "./2_motif.outs/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz"
## res_motif <- fread(fn_motif, header=T, data.table=F)
## res_motif <- res_motif%>%
##     mutate(comb=paste(Cluster, psycho_variable, sep="_"))

## res_motif <- res_motif%>%filter(psycho_variable%in%varSel)


## ###
## ###
## for ( i in 1:nrow(summ)){    
##    ##  
##    var0 <- summ$var[i]
    
##     MCl <- summ$MCls[i]    
##     cl_rna <- gsub("^R", "C", summ$rna_cluster[i])
##     cl_atac <- gsub("^A", "C", summ$atac_cluster[i])
##     ii <- paste(cl_rna, var0, sep="_")
##     ii2 <- paste(cl_atac, var0, sep="_")

##     ##
    
##     ###
##     ### motif activity results
##     res0 <- res_motif%>%filter(comb==ii2)%>%
##         dplyr::select(motif_name=gene, estimate, stderror, padj_activity=padj_t)

##     df0 <- res_enrich%>%filter(comb==ii)%>%
##         dplyr::select(motif_name, log_odds, log_lower, log_upper, pval_enrich=pval, padj_enrich=padj)
    
##     df0 <- df0%>%inner_join(res0, by="motif_name")  ##%>%filter(padj_activity<0.1, padj<0.01)
    
##     df0$cluster_RNA <- cl_rna
##     df0$comb <- ii    
##     df0$cluster_atac <- cl_atac
##     df0$celltype <- MCl
##     df0$psycho_variable <- var0

##     #####
##     #####
##     plotDF <- df0%>%filter(padj_activity<0.1)%>%    
##          mutate(motif2=gsub("_.*", "", motif_name), direct=ifelse(estimate>0, "gr1", "gr2"),
##                direct=ifelse(padj_enrich < 0.1, direct, "gr0")) 
##     plotDF <- plotDF%>%distinct(motif2, .keep_all=T)    
## var_name <- var_full[var0]
## p0 <- ggplot(plotDF, aes(x = log_odds, y = fct_reorder(motif2, log_odds), color=factor(direct)  ))+
##     geom_errorbarh(aes(xmax = log_upper, xmin = log_lower), linewidth = 0.5, height = 0.2)+
##     geom_point(shape = 19, size = 1.5)+
##     scale_color_manual(values=c("gr0"="grey", "gr1"="#F8766D", "gr2"="#00A9FF"),
##         labels=c("gr0"="Not enrich", "gr1"="Increase", "gr2"="Decrease"),
##         guide=guide_legend(override.aes=list(size=2)))+    
##     geom_vline(aes(xintercept = 0), linewidth = 0.25, linetype = "dashed")+
##     ggtitle(paste(var_name , " (", MCl, ")", sep=""))+
##     xlab("log odds ratio")+xlim(-6, 4)+
##     theme_bw()+
##     theme(axis.title.y = element_blank(),
##           axis.title.x = element_text(size=10),
##           axis.text.x = element_text(size=10),
##           axis.text.y = element_text(size=9),
##           plot.title = element_text(size=12, hjust=0.5), 
##           legend.title = element_blank(),
##           legend.text = element_text(size=10),
##           legend.key.size = unit(0.35, "cm"))

## ###
## ###
## nm <- nrow(plotDF)
## cat(var0, cl_rna, nm, "motifs", "\n")    
## if ( nm > 40){
##   ##  
##   figfn <- paste(outdir2, "Figure0_", var0, "_", cl_rna, "_enrich.forest.png", sep="")
##   ggsave(figfn, p0, width = 620, height = 650, units = "px", dpi=120)
##   ##  
## } else if (nm > 20){
##   ##  
##   figfn <- paste(outdir2, "Figure0_", var0, "_", cl_rna, "_enrich.forest.png", sep="")
##   ggsave(figfn, p0, width = 620, height = 600, units = "px", dpi=120)
## } else{
##   ###
##   figfn <- paste(outdir2, "Figure0_", var0, "_", cl_rna, "_enrich.forest.png", sep="")
##   ggsave(figfn, p0, width = 500, height = 300, units = "px", dpi=120)    
## }


## }    
###
### END 


##############################
### Heatmap for plots 
##############################

### Script will be used final publication
### I will copy these in the file '888_pub_figure.R' later
### Modified Jun-26-2025 By JW

###
### read data

outdir <- "./3_regulatory.outs/correct_excludeX/"


varSel <- c("ISEL_Mean", "PSS_all_mean")
var_full <- c("ISEL_Mean"="Social Support", "PSS_all_mean"="Perceived Stress")

###
### cluster 
fn0 <- "./3_regulatory.outs/list_cluster_var_nDEGs20_nTF1_new.txt" 
df_cl <- read.table(fn0, header=T, sep="\t")
df_cl <- df_cl%>%filter(var%in%varSel, nDEG > 50)%>%arrange(var)


###
### 
fn_enrich <- "./3_regulatory.outs/correct_excludeX/enrich_motif/2_comb_enrich.motif.txt.gz"
res_enrich <- fread(fn_enrich, header=T, sep="\t")
## res_enrich <- res_enrich%>%filter(psycho_variable%in%varSel)

###
fn2 <- "./2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz"
res <- fread(fn2, header=T, data.table=F, sep="\t")
res <- res%>%filter(Cluster!="C9")
res <- res%>%filter(psycho_variable%in%varSel) ##, Cluster!="C9")


###
####
motif_sel <- lapply(1:nrow(df_cl), function(i){
   ## 
   cl_rna <- gsub("^R", "C", df_cl$rna_cluster[i])
   cl_atac <- gsub("^A", "C", df_cl$atac_cluster[i])
   var0 <- df_cl$var[i]

   ### 
   enrich_motif <- res_enrich%>%
       filter(cluster_rna==cl_rna, psycho_variable==var0, odds>1, padj<0.1)%>%pull(motif_name)%>%unique()    
   DAMs <- res%>%filter(Cluster==cl_atac, psycho_variable==var0, padj_t<0.1)%>%pull(gene)%>%unique()

   motif2 <- intersect(enrich_motif, DAMs)
   motif2
})    

motif_sel <- unique(unlist(motif_sel))  ### 59 motifs 

## motif_sel <- motif_sel[grepl("IRF|STAT|REL", motif_sel)] ### 14 motif
 

###
### make 59 motif table
## res_all <- map_dfr(1:nrow(df_cl), function(i){
##     ##
##     cl_rna <- gsub("^R", "C", df_cl$rna_cluster[i])
##     cl_atac <- gsub("^A", "C", df_cl$atac_cluster[i])
##     var0 <- df_cl$var[i]

##    ### 
##    res0_enrich <- res_enrich%>%
##        filter(cluster_rna==cl_rna, cluster_atac==cl_atac, psycho_variable==var0, motif_name%in%motif_sel)%>%
##        dplyr::select(cluster_rna, psycho_variable, motif_name, odds,
##                     pval_odds=pval, padj_odds=padj, log_odds, log_lower, log_upper)
##    ### 
##    res0 <- res%>%filter(Cluster==cl_atac, psycho_variable==var0, gene%in%motif_sel)%>%
##         dplyr::select(cluster_atac=Cluster,
##                       motif_name=gene, beta_diff=estimate, std_diff=stderror, zval_diff=statistic,
##                       pval_diff=pval_t, padj_diff=padj_t)
##    res0 <- res0_enrich%>%inner_join(res0, by="motif_name") 
## })

## ###
## ### save 
## opfn <- paste(outdir, "TableS_motif.list.txt", sep="")
## write.table(res_all, file=opfn, quote=F, sep="\t", row.names=F)




#########################################
### all cluster and variable in df_cl
#################################
## res_all <- map_dfr(1:nrow(df_cl), function(i){
##     ##
##     cl_rna <- gsub("^R", "C", df_cl$rna_cluster[i])
##     cl_atac <- gsub("^A", "C", df_cl$atac_cluster[i])
##     var0 <- df_cl$var[i]

##    ### 
##    res0_enrich <- res_enrich%>%
##        filter(cluster_rna==cl_rna, cluster_atac==cl_atac, psycho_variable==var0)%>%
##        dplyr::select(cluster_rna, psycho_variable, motif_name, odds,
##                     pval_odds=pval, padj_odds=padj, log_odds, log_lower, log_upper)
##    ### 
##    res0 <- res%>%filter(Cluster==cl_atac, psycho_variable==var0)%>%
##         dplyr::select(cluster_atac=Cluster,
##                       motif_name=gene, beta_diff=estimate, std_diff=stderror, zval_diff=statistic,
##                       pval_diff=pval_t, padj_diff=padj_t)
##    res0 <- res0_enrich%>%inner_join(res0, by="motif_name")
##    res0 <- res0%>%filter(padj_odds<0.1, padj_diff<0.1)
##    res0 
## })
   
## ###
## ### save 
## opfn <- paste(outdir, "TableS_motif_all.infor.txt", sep="")
## write.table(res_all, file=opfn, quote=F, sep="\t", row.names=F)


###
###


#################
### dot plots
#################



plotDF <- map_dfr(1:nrow(df_cl), function(i){
   ###
   ###
   cl_rna <- gsub("^R", "C", df_cl$rna_cluster[i])
   cl_atac <- gsub("^A", "C", df_cl$atac_cluster[i])
   var0 <- df_cl$var[i]

   ### 
   df0 <- res_enrich%>%
       filter(cluster_rna==cl_rna, cluster_atac==cl_atac, psycho_variable==var0, motif_name%in%motif_sel)%>%
       dplyr::select(motif_name, log_odds, log_odds_se=se, pval_enrich=pval, padj_enrich=padj)

   df2 <- res%>%filter(Cluster==cl_atac, psycho_variable==var0, gene%in%motif_sel)%>%
       dplyr::select(motif_name=gene, beta_eff=estimate, beta_se=stderror, zval=statistic, pval=pval_t, padj=padj_t)

   ##
   ## gene2 <- gsub(motif_sel 
   ## df3 <- res%>%filter(cluster==cl_rna, var==var0, identifier%in% 

   ##
   df0 <- df0%>%left_join(df2, by="motif_name")
   df0$rna_cluster <- df_cl$rna_cluster[i]
   df0$atac_cluster <- df_cl$atac_cluster[i]
   df0$MCls <- df_cl$MCls[i] 
   df0$psycho_var <- var0

   ###
   ###
   df0
})    
 


#############################################
### Heatmap for Effects on motif activity
###############################################


var_full <- c("ISEL_Mean" = "Social Support", "PSS_all_mean" = "Perceived Stress")
plotDF2 <- plotDF%>%
    mutate(is_sig=as.integer(padj<0.1), zval2 = zval * is_sig,
           var2 = var_full[psycho_var], comb = paste(var2, MCls, sep = "_"))



###
### version-1, heatmap show all values 
mat <- plotDF2%>%
    pivot_wider(id_cols=motif_name, names_from=comb, values_from=zval, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()

rnSel <- gsub("_.*", "", rownames(mat))
rownames(mat) <- rnSel

rn2 <- sort(rnSel)
comb <- sort(colnames(mat))
mat2 <-  mat[rn2, comb]


###
### setting color
max_val <- max(abs(max(mat)), abs(min(mat)))
max_val <- round(max_val, 2)
min_val <- -max_val
mycol <- colorRamp2(seq(min_val, max_val, length.out=20), colorRampPalette(c("blue", "white", "red"))(20))

###
### column annotation
col_var <- c("Social Support"="#29465B", "Perceived Stress"="#E34234")
col_MCls <- c("R0"="#E60A0A", "R1"="#F9A602", "R2"="#A65192", "R3"="#0570B0", "R4"="#00AF72", "R5"="#87CEFA")
anno_df <- data.frame(variable=gsub("_A.*", "", colnames(mat2)),
    cluster=gsub("_.*", "", gsub(".*_R", "R", colnames(mat2))) )

col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, cluster=col_MCls),
    show_annotation_name = F, show_legend = F)                  




###
### row annotation

fn <- paste(outdir, "Figure2.1_row_cluster.tsv", sep="")
df_cluster <- read_tsv(fn, show_col_types = F)

anno_row <- data.frame(motif_name = rownames(mat2))%>%
    left_join(df_cluster, by="motif_name")%>%
    dplyr::select(cluster)
col_cluster <- c("cluster1"="#1b9e77", "cluster2"="#d95f02", "cluster3"="#e7298a", "cluster4"="#7570b3")
row_ha <- rowAnnotation(df = anno_row, col=list(cluster=col_cluster),
    show_annotation_name=F, show_legend=F)

###
### make raw heatmap and use all zval 
p0 <- Heatmap(mat2, col=mycol, cluster_rows=T, row_split = 4,
   cluster_columns=F,
   column_names_rot = -45, column_names_gp = gpar(fontsize=9),
   row_names_gp = gpar(fontsize=9),
   top_annotation = col_ha, right_annotation = row_ha,
   heatmap_legend_param = list(title = bquote(~italic(Z)~"-score"),
      at = round(seq(min_val, max_val, length.out=4), 1),
      grid_width = grid::unit(0.38, "cm"), legend_height = grid::unit(4, "cm")))

###
### save
figfn <- paste(outdir, "Figure2.1_zscore.heatmap0.png", sep="")
png(figfn, width = 600, height=950, res = 120)
set.seed(0)
p0 <- draw(p0)
dev.off()
  


###
### get clusters and it needs run in the first time. Only run once.  
## hmap <- Heatmap(mat2, cluster_rows = T, cluster_columns = F, row_split = 4)
## set.seed(0)
## hmap <- draw(hmap)
## cl <- row_order(hmap)


## df_cluster <- NULL
## for ( i in 1:length(cl)){
##     ##
##     df0 <- data.frame(cluster = paste("cluster", i, sep=""),
##         motif_name = rownames(mat2)[cl[[i]]])
##     df_cluster <- rbind(df_cluster, df0)
## }

## opfn <- paste(outdir, "Figure2.1_row_cluster.tsv", sep="")
## write_tsv(df_cluster, file = opfn)


###########################################################
### version2, heatmap show significant values 10% fdr
###########################################################



###
fn <- paste(outdir, "Figure2.1_row_cluster.tsv", sep="")
df_cluster <- read_tsv(fn, show_col_types = F)


###
mat3 <- plotDF2%>%
    pivot_wider(id_cols=motif_name, names_from=comb, values_from=zval2, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()

rnSel <- gsub("_.*", "", rownames(mat3))
rownames(mat3) <- rnSel

rn2 <- sort(rnSel)
comb <- sort(colnames(mat3))
mat3 <-  mat3[rn2, comb]

mat3 <- mat3[df_cluster$motif_name, ]



###
### column annotation
col_var <- c("Social Support"="#29465B", "Perceived Stress"="#E34234")
col_MCls <- c("R0"="#E60A0A", "R1"="#F9A602", "R2"="#A65192", "R3"="#0570B0", "R4"="#00AF72", "R5"="#87CEFA")
anno_df <- data.frame(variable=gsub("_A.*", "", colnames(mat3)),
    cluster=gsub("_.*", "", gsub(".*_R", "R", colnames(mat3))) )

col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, cluster=col_MCls),
    show_annotation_name = F, show_legend = F)                  



##
## hight row names  
row_hlight <- c("BCL6", "IRF2", "IRF4", "IRF7", "STAT1", "STAT1..STAT2")
geneSel <- read.table("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem/figures/DEGs_TF/motifs/TF-genes_expressed_in_RNA.txt", header=T)
row_hlight <- gsub("_.*", "", geneSel$motif_name)


iwhich <- which(rownames(mat3)%in%row_hlight)
fontcolors <- rep("black", nrow(mat3))
fontcolors[iwhich] <- "black"

fontfaces <- rep("plain", nrow(mat3))
fontfaces[iwhich] <- "bold"


####
### dendrograms 
row_dend <- as.dendrogram(hclust(dist(mat2[rownames(mat3),])))


###
### row annotation 
anno_row <- data.frame(motif_name = rownames(mat3))%>%
    left_join(df_cluster, by="motif_name")%>%
    dplyr::select(cluster)
col_cluster <- c("cluster1"="#1b9e77", "cluster2"="#d95f02", "cluster3"="#e7298a", "cluster4"="#7570b3")
row_ha <- rowAnnotation(rows = anno_text(rownames(mat3), gp=gpar(fontsize=9, fontface=fontfaces, col = fontcolors)),
    df = anno_row, col=list(cluster=col_cluster),
    show_annotation_name=F, show_legend=F)



###
### setting color
max_val <- max(abs(max(mat3)), abs(min(mat3)))
max_val <- round(max_val, 2)
min_val <- -max_val
mycol <- colorRamp2(seq(min_val, max_val, length.out=20), colorRampPalette(c("blue", "white", "red"))(20))


###
### make heatmap 
p1 <- Heatmap(mat3, col=mycol,
   cluster_rows=row_dend, row_order=rownames(mat3), row_dend_reorder = FALSE, 
   cluster_columns=F,
   column_names_rot = -45, column_names_gp = gpar(fontsize=10),
   show_row_names = FALSE, ##row_names_gp = gpar(fontsize=9),
   top_annotation = col_ha, right_annotation = row_ha,  
   heatmap_legend_param = list(title = bquote(~italic(Z)~"-score"),
      at = round(seq(min_val, max_val, length.out=4), 1),
      grid_width = grid::unit(0.38, "cm"), legend_height = grid::unit(4, "cm")))

###
### save
figfn <- paste(outdir, "Figure2.1_zscore.heatmap.png", sep="")
png(figfn, width = 580, height=1050, res = 120)
set.seed(0)
p1 <- draw(p1)
dev.off()

## x <- read.table("./3_regulatory.outs/correct_excludeX/motif_list_from_Ali.txt", header=T, sep="\t")



###########################
### log odds ratio
###########################
 
var_full <- c("ISEL_Mean" = "Social Support", "PSS_all_mean" = "Perceived Stress")
plotDF2 <- plotDF%>%
    mutate(is_sig=as.integer(padj_enrich<0.1), log_odds=ifelse(is_sig, log_odds, 0), 
           var2 = var_full[psycho_var], comb = paste(var2, MCls, sep = "_"))

mat <- plotDF2%>%
    pivot_wider(id_cols=motif_name, names_from=comb, values_from=log_odds, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()

rnSel <- gsub("_.*", "", rownames(mat))
rownames(mat) <- rnSel

rn2 <- sort(rnSel)
comb <- sort(colnames(mat))
mat2 <- mat[rn2, comb]


###
### get the cluster 
fn <- "./3_regulatory.outs/correct_excludeX/Figure2.1_row_cluster.tsv"
df_cluster <- read_tsv(fn, show_col_types = F)

mat2 <- mat2[df_cluster$motif_name, ]

###
### row annotation 
anno_row <- data.frame(motif_name = rownames(mat2))%>%
    left_join(df_cluster, by="motif_name")%>%
    dplyr::select(cluster)
col_cluster <- c("cluster1"="#1b9e77", "cluster2"="#d95f02", "cluster3"="#e7298a", "cluster4"="#7570b3")
row_ha <- rowAnnotation(df = anno_row, col=list(cluster=col_cluster),
    show_annotation_name=F, show_legend=F)


###
### setting color
max_odds <- max(plotDF$log_odds)
mycol <- colorRamp2(seq(0, max_odds, length.out=20),
                    colorRampPalette(c("white", brewer.pal(n=7, name="RdPu")[2:7]))(20))

###
### column annotation
col_var <- c("Social Support"="#29465B", "Perceived Stress"="#E34234")
col_cl <- c("R0"="#E60A0A", "R1"="#F9A602", "R2"="#A65192", "R3"="#0570B0", "R4"="#00AF72", "R5"="#87CEFA")
anno_df <- data.frame(variable=gsub("_A.*", "", colnames(mat2)),
    cluster=gsub("_.*", "", gsub(".*_R", "R", colnames(mat2))) )

col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, cluster=col_cl),
    show_annotation_name = F, show_legend = F)                        


###
### make a heatmap  
p2 <- Heatmap(mat2, col=mycol, cluster_rows=F, cluster_columns=F, na_col="white",
   column_names_rot = -45, column_names_gp = gpar(fontsize=9),
   row_names_gp = gpar(fontsize=9),
   top_annotation = col_ha, right_annotation = row_ha, 
   heatmap_legend_param = list(title = bquote(log~"odds ratio"),
      at = round(seq(0, max_odds, length.out=3)),
      grid_width = grid::unit(0.38, "cm"), legend_height = grid::unit(4, "cm")))

###
### save
figfn <- paste(outdir, "Figure2.2_odds.heatmap.png", sep="")
png(figfn, width = 580, height=950, res = 120)
p2 <- draw(p2)
dev.off()


#####################################################
### version-3, Heatmap across cell types using 59 motifs 
######################################################

var_full <- c("ISEL_Mean" = "Social Support", "PSS_all_mean" = "Perceived Stress")
MCls_name <- c("C0" = "A0_R0_CD4+", "C1" = "A1_R2_CD8+", "C2" = "A2_R3_NK", "C3"="A3_R4_Monocyte",
    "C4" = "A4_R5_B", "C5" = "A5_R1_CD4+", "C6"="A6_R1_CD4+", "C7"="A7_R1_CD4+", "C8"="A8_R6_DC")

res2 <- res%>%filter(gene%in%motif_sel)%>%
    mutate(MCls = MCls_name[Cluster], var2 = var_full[psycho_variable], comb = paste(var2, MCls, sep = "_"),
           motif_name = gsub("_.*", "", gene),
           is_sig = ifelse(pval_t<0.05, 1, 0), zval2 = statistic * is_sig)
 
##
## make plot data 
mat <- res2%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=statistic, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()

### sort column
## comb <- sort(colnames(mat))
col_df <- data.frame(comb=sort(colnames(mat)))%>%
    mutate(MCls=gsub(".*_A", "A", comb),
           MCl2=gsub(".*_R", "R", comb), psycho_var=gsub("_A.*", "", comb))

## col_df <- col_df%>%arrange(MCl2, psycho_var) 
col_df <- col_df%>%arrange(psycho_var, MCl2) 

    
mat <- mat[, col_df$comb]



###
### setting color
max_val <- max(abs(max(mat)), abs(min(mat)))
max_val <- round(max_val, 2)
min_val <- -max_val
mycol <- colorRamp2(seq(min_val, max_val, length.out=20), colorRampPalette(c("blue", "white", "red"))(20))

###
### column annotation
col_var <- c("Social Support"="#29465B", "Perceived Stress"="#E34234")
col_MCls <- c("R0"="#E60A0A", "R1"="#F9A602", "R2"="#A65192", "R3"="#0570B0", "R4"="#00AF72", "R5"="#87CEFA")
anno_df <- data.frame(variable=gsub("_A.*", "", colnames(mat)),
    cluster=gsub("_.*", "", gsub(".*_R", "R", colnames(mat))) )

col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, cluster=col_MCls),
    show_annotation_name = F, show_legend = F)                  



###
### row annotation

## fn <- paste(outdir, "Figure2.3_row_cluster.tsv", sep="")
## df_cluster <- read_tsv(fn, show_col_types = F)
 
## anno_row <- data.frame(motif_name = rownames(mat))%>%
##     left_join(df_cluster, by="motif_name")%>%
##     dplyr::select(cluster)
## col_cluster <- c("cluster1"="#1b9e77", "cluster2"="#d95f02", "cluster3"="#e7298a", "cluster4"="#7570b3")
## row_ha <- rowAnnotation(df = anno_row, col=list(cluster=col_cluster),
##     show_annotation_name=F, show_legend=F)


###
row_dend <- dendsort(hclust(dist(mat)))


###
### make raw heatmap and use all zval  
p0 <- Heatmap(mat, col=mycol, cluster_rows=row_dend, row_split = 4,
   cluster_columns=F,
   column_names_rot = -45, column_names_gp = gpar(fontsize=9),
   row_names_gp = gpar(fontsize=9),
   top_annotation = col_ha, ###right_annotation = row_ha,
   heatmap_legend_param = list(title = bquote(~italic(Z)~"-score"),
      at = round(seq(min_val, max_val, length.out=4), 1),
      grid_width = grid::unit(0.38, "cm"), legend_height = grid::unit(4, "cm")))

###
### save
figfn <- paste(outdir, "Figure2.4_zscore.heatmap0.png", sep="")
png(figfn, width = 600, height=950, res = 120)
set.seed(0)
p0 <- draw(p0)
dev.off()
  


###
### get clusters and it needs run in the first time. Only run once.  
hmap <- Heatmap(mat, cluster_rows = row_dend, cluster_columns = F, row_split = 4)
set.seed(0)
hmap <- draw(hmap)
cl <- row_order(hmap)


df_cluster <- NULL
for ( i in 1:length(cl)){
    ##
    df0 <- data.frame(cluster = paste("cluster", i, sep=""),
        motif_name = rownames(mat)[cl[[i]]])
    df_cluster <- rbind(df_cluster, df0)
}

opfn <- paste(outdir, "Figure2.4_row_cluster.tsv", sep="")
write_tsv(df_cluster, file = opfn)



##############################
### use nominal p values 
##############################

 
##
## make plot data 
mat2 <- res2%>%pivot_wider(id_cols=motif_name, names_from=comb, values_from=zval2, values_fill=0)%>%
    column_to_rownames(var="motif_name")%>%
    as.matrix()


###
### sort column 
col_df <- data.frame(comb=sort(colnames(mat2)))%>%
    mutate(MCls=gsub(".*_A", "A", comb),
           MCl2=gsub(".*_R", "R", comb), psycho_var=gsub("_A.*", "", comb))

## col_df <- col_df%>%arrange(MCl2, psycho_var)

col_df <- col_df%>%arrange(psycho_var, MCl2) 

mat2 <- mat2[, col_df$comb]



###
### setting color
max_val <- max(abs(max(mat2)), abs(min(mat2)))
max_val <- round(max_val, 2)
min_val <- -max_val
mycol <- colorRamp2(seq(min_val, max_val, length.out=20), colorRampPalette(c("blue", "white", "red"))(20))

###
### column annotation
col_var <- c("Social Support"="#29465B", "Perceived Stress"="#E34234")
col_MCls <- c("R0"="#E60A0A", "R1"="#F9A602", "R2"="#A65192", "R3"="#0570B0", "R4"="#00AF72", "R5"="#87CEFA")
anno_df <- data.frame(variable=gsub("_A.*", "", colnames(mat2)),
    MCls=gsub("_.*", "", gsub(".*_R", "R", colnames(mat2))) )

col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, MCls=col_MCls),
    show_annotation_name = F, show_legend = F)                  



###
### row annotation

## fn <- paste(outdir, "Figure2.3_row_cluster.tsv", sep="")
## df_cluster <- read_tsv(fn, show_col_types = F)
 
## anno_row <- data.frame(motif_name = rownames(mat))%>%
##     left_join(df_cluster, by="motif_name")%>%
##     dplyr::select(cluster)
## col_cluster <- c("cluster1"="#1b9e77", "cluster2"="#d95f02", "cluster3"="#e7298a", "cluster4"="#7570b3")
## row_ha <- rowAnnotation(df = anno_row, col=list(cluster=col_cluster),
##     show_annotation_name=F, show_legend=F)


 
## dendrograms from all zvalue
###library(dendsort)
row_dend <- dendsort(hclust(dist(mat2)))
 

###
### make raw heatmap and use all zval  
p0 <- Heatmap(mat2, col=mycol, cluster_rows=row_dend,  row_split = 4,
   cluster_columns=F,
   column_names_rot = -45, column_names_gp = gpar(fontsize=9),
   row_names_gp = gpar(fontsize=9),
   top_annotation = col_ha, ###right_annotation = row_ha,
   heatmap_legend_param = list(title = bquote(~italic(Z)~"-score"),
      at = round(seq(min_val, max_val, length.out=4), 1),
      grid_width = grid::unit(0.38, "cm"), legend_height = grid::unit(4, "cm")))

###
### save
figfn <- paste(outdir, "Figure2.4_zscore.heatmap.png", sep="")
png(figfn, width = 600, height=950, res = 120)
set.seed(0)
p0 <- draw(p0)
dev.off()
















###################################
### effects on gene expression 
###################################

## motif2 <- gsub("_.*", "", motif_sel)
## motif2 <- unique(unlist(str_split(motif2, "\\.\\.")))


## ###
## ### gene expression
## fn_DEG <- "./2_motif.outs/correct_excludeX/summary_Cluster_res0.07_psycho/Ali_RNA.diff.DESeq.txt.gz"
## res_DEG <- fread(fn_DEG, header=T, data.table=F)
## res_DEG <- res_DEG%>%filter(var%in%varSel, identifier%in%motif2)%>%
##     dplyr::rename(gene=identifier)   ### 37 genes


## ###
## ### plot data
## var_full <- c("ISEL_Mean" = "Social Support", "PSS_all_mean" = "Perceived Stress")
## MCls_name <- c("C0"="R0_CD4+", "C1"="R1_CD4+", "C2"="R2_CD8+", "C3"="R3_NK",
##                "C4"="R4_Monocyte", "C5"="R5_B", "C6"="R6_DC") 
## plotDF2 <- res_DEG%>%
##     mutate(is_sig=as.integer(padj<0.1), beta_eff=ifelse(padj<0.1, logFC, 0),
##            var2 = var_full[var], MCls=MCls_name[cluster], comb = paste(var2, MCls, sep = "_"),
##            comb2 = paste(gsub("_.*", "", MCls), var, sep="_"))

## mat <- plotDF2%>%filter(comb2%in%df_cl$comb)%>%
##     pivot_wider(id_cols=gene, names_from=comb, values_from=beta_eff, values_fill=0)%>%
##     column_to_rownames(var="gene")%>%
##     as.matrix()
 
## rnSel <- gsub("_.*", "", rownames(mat))
## rownames(mat) <- rnSel

## rn2 <- sort(rnSel)
## comb <- sort(colnames(mat))
## mat2 <- mat[rn2, comb]


## ###
## ### setting color
## max_val <- max(abs(round(max(mat), 2)), abs(round(min(mat), 2)))
## min_val <- -max_val
## mycol <- colorRamp2(seq(min_val, max_val, length.out=20), colorRampPalette(c("blue", "white", "red"))(20))

## ###
## ### column annotation
## col_var <- c("Social Support"="#29465B", "Perceived Stress"="#E34234")
## col_cl <- c("R0"="#E60A0A", "R1"="#F9A602", "R2"="#A65192", "R3"="#0570B0", "R4"="#00AF72", "R5"="#87CEFA")
## anno_df <- data.frame(variable=gsub("_R.*", "", colnames(mat2)),
##     cluster=gsub("_.*", "", gsub(".*_R", "R", colnames(mat2))) )

## col_ha <- HeatmapAnnotation(df = anno_df, col = list(variable=col_var, cluster=col_cl),
##     show_annotation_name = F, show_legend = F)                  

 
## p2 <- Heatmap(mat2, col=mycol, cluster_rows=T, cluster_columns=F,
##    column_names_rot = -45, column_names_gp = gpar(fontsize=9),
##    row_names_gp = gpar(fontsize=9),
##    top_annotation = col_ha, 
##    heatmap_legend_param = list(title = "Effects",
##       at = round(seq(min_val, max_val, length.out=4), 1),
##       grid_width = grid::unit(0.38, "cm"), legend_height = grid::unit(4, "cm")))

## ###
## ### save
## figfn <- paste(outdir, "Figure2.3_effects_gene.heatmap.png", sep="")
## png(figfn, width = 600, height=950, res = 120)
## p2 <- draw(p2)
## dev.off()





################################################
### glm enrichment results 
################################################


##################################
### combine glm results 
##################################

## outdir2 <- paste(outdir, "enrich_motif_glm/", sep="")
## fn_ls <- list.files(outdir2, "glm.txt")

## res_all <- map_dfr(fn_ls, function(ii){
##     ##
##     fn <- paste(outdir2, ii, sep="")
##     res <- read.table(fn, header = T, sep = "\t")
 
##     res <- res%>%arrange(pval)
##     n0 <- nrow(res)
##     res <- res%>%mutate(observed = -log10(pval), expected = -log10(ppoints(n0)))
##     ###  
##     res    
## })
 

## opfn <- paste(outdir2, "2_comb.txt.gz", sep="")
## fwrite(res_all, file=opfn, na=NA, sep="\t")


## ##################
## ### qq plots
## ##################
 
## fn_enrich <- paste(outdir2, "2_comb.txt.gz", sep="")
## plotDF <- fread(fn_enrich, header=T, sep="\t")

## p0 <- ggplot(plotDF, aes(x = expected, y = observed, color=cluster_rna))+
##     rasterise(geom_point(size=0.06), dpi=300)+
##     geom_abline(color="grey")+
##     scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
##         "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF"),
##          guide=guide_legend(override.aes=list(size=2)))+
##     facet_wrap(~comb, scales="free_y", ncol=6)+  ## ncol=7
##     xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
##     ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
##     theme_bw()+
##     theme(legend.title=element_blank(),
##           legend.text=element_text(size=9),
##           legend.key.size=grid::unit(0.4, "cm"),
##           axis.title=element_text(size=9),
##           axis.text=element_text(size=9),
##           strip.text=element_text(size=9))
## ###
## ###
## figfn <- paste(outdir2, "Figure2_enrich.glm.qq.png", sep="")
## ggsave(figfn, p0, width=1200, height=700, units="px", dpi=120)




###
### show 







## ## fn <- "../1.2_ArchR_process/4_motif_output/Motif_jaspar2022.motifinfor.ArchR.txt"
## ## motif_df <- read.table(fn, header=T, sep="\t")


## ###
## ### summary 
## fn <- paste(outdir, "2.1_motif.full.txt", sep="")
## df_full <- read.table(fn, header=T, sep="\t")

## ##
## fn_summ <- paste(outdir, "list_var_clust_nDEGs20_nTF1.txt", sep="")
## summ <- read.table(fn_summ, header=T, sep="\t")
## summ <- summ%>%mutate(comb=paste(Cluster_scRNA, variable, sep="_"))
 

## summ$olap <- map_dbl(1:nrow(summ), function(i){
##     ###
##     ii <- summ$comb[i]
##     n0 <- df_full%>%filter(comb==ii, padj<0.01, padj_activity<0.1)%>%pull(motif_name)%>%unique()%>%length()
##     n0
## })

## ###
## opfn <- paste(outdir, "2.1_motif.summ.tsv", sep="")
## write_tsv(summ, file=opfn)


###################################################
### LASSO select TFs are enriched 
####################################################





##################################################################
### cv.glmnet to determine which TF motifs are in DEGs  
#################################################################

 
## outdir2 <- "./3_regulatory.outs/enrich_motif_glmnet/"
## if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings = F, recursive = T)



## ###
## ### motif-gene match matrix 
## fn <- "../1.2_ArchR_process/4_motif_output/motif_anno/2_comb_motif.annot.txt.gz"
## x <- fread(fn, header=T, data.table=F, sep="\t")

## ###
## ###
## mat <- as.matrix(table(x$gene, x$motif_name))
## spMat <- as(mat, "sparseMatrix")

## opfn <- "./3_regulatory.outs/0_gene_motif.mat.rds"
## write_rds(spMat, file=opfn)


## ###
## ### cluster 
## fn_cl <- "./3_regulatory.outs/list_var_clust_nDEGs20_nTF1.txt"
## df_cl <- read.table(fn_cl, header=T, sep="\t")
## df_cl <- df_cl%>%
##     mutate(comb=paste(Cluster_scRNA, variable, sep="_"))%>%
##     filter(nTFs>1)%>%dplyr::select(comb, nTFs)%>%
##     arrange(comb)
    




## ## infn <- paste(outdir2, "C3_BPd_avg.glmnet.rds", sep="")
## ## res <- read_rds(infn)

## ## res <- res%>%mutate(comb=paste(psycho_variable, "_rna", cluster_rna, sep=""))%>%
## ##     group_by(comb)%>%mutate(padj = p.adjust(pval, "BH"))%>%ungroup()

## ###
## ###
## fn_ls <- list.files(outdir2, "glmnet.rds")

## ## fn_cl <- "./3_regulatory.outs/list_var_clust_nDEGs20_nTF1.txt"
## ## df_cl <- read.table(fn_cl, header=T, sep="\t")
## ## df_cl <- df_cl%>%mutate(comb = paste(Cluster_scRNA, variable, sep="_"))
## ## comb <- sort(df_cl$comb)

## res_all <- map_dfr(fn_ls, function(ii){
##    ###
##    ii0 <- gsub(".glmnet.rds", "", ii)
    
##    fn <- paste(outdir2, ii, sep="")
##    mat <- read_rds(fn)
    
##    ##
##    bse <- apply(mat[,-1], 1, sd)
##    b <- mat[,1]

##   ###  results  
##   res <- data.frame(motif_name=rownames(mat), estimate = b, stderror = bse)
##   res <- res%>%mutate(zval = ifelse(estimate==0, 0, estimate/stderror), pval = 2*pnorm(abs(zval), lower.tail = F))

##   ## res <- data.frame(motif_name = rownames(mat), estimate = b)  
##   res$comb <- ii0
##   res$cluster_rna <- str_replace(ii0, "_.*", "")
##   res$psycho_var <- str_replace(ii0, "C._", "")    
##   ##
##   res
## })


## opfn <- paste(outdir2, "2_comb.txt.gz", sep="")
## fwrite(res_all, file=opfn, na=NA, sep="\t")


## ###
## ###
## summ <- res_all%>%group_by(comb)%>%
##     summarize(nm_greater0=sum(estimate>0), nm_sig=sum(estimate>0&pval<0.05), .groups="drop")%>%
##     left_join(df_cl, by="comb")
## summ <- summ[, c("comb", "nTFs", "nm_greater0", "nm_sig")]

## opfn <- paste(outdir2, "2_glmnet.summ.xlsx", sep="")
## write.xlsx(summ, file=opfn)




## ###############################
## ### narrow down motifs
## ###############################


## fn <- "./3_regulatory.outs/list_var_clust_nDEGs20_nTF1.txt"
## summ <- read.table(fn, header=T, sep="\t")


## ##
## ### combine from RNA cluster and psycho_variable
## fn_enrich <- "./3_regulatory.outs/Celltype/enrich_motif/2_comb_enrich.motif.txt.gz"
## res_enrich <- fread(fn_enrich, header=T, sep="\t")


## ###
## ### diff motif, combine from ATAC cluster and psycho_variable
## fn_motif <- "./2_motif.outs/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz"
## res_motif <- fread(fn_motif, header=T, data.table=F)
## res_motif <- res_motif%>%
##     mutate(comb=paste(Cluster, psycho_variable, sep="_"))


## ###
## ### motif from glmnet from RNA cluster and psycho_variable
## fn_glmnet <- "./3_regulatory.outs/enrich_motif_glmnet/2_comb.txt.gz"
## res_glmnet <- fread(fn_glmnet, header=T, data.table=F)
## res_glmnet <- res_glmnet%>%filter(estimate > 0)%>%dplyr::select(comb, motif_name)





## ###
## ### motif
## df_full <- NULL
## for ( i in 1:nrow(summ)){
##     ##
##     var0 <- summ$variable[i]
##     var_full <- summ$varname[i]
    
##     MCl <- summ$celltype[i]    
##     cl_rna <- summ$Cluster_scRNA[i]
##     cl_atac <- summ$Cluster_scATAC[i]
##     ii <- paste(cl_rna, var0, sep="_")
##     ii2 <- paste(cl_atac, var0, sep="_")

##     n_DEG <- summ$nDEGs[i]
##     n_TF <- summ$nTFs[i]
##     ##
    
##     ###
##     ### motif activity results
##     res0 <- res_motif%>%filter(comb==ii2)%>%
##         dplyr::select(motif_name=gene, estimate, stderror, padj_activity=p.adjusted)%>%
##         filter(padj_activity<0.1)

##     ## if ( n_TF > 1){
##     ##     ##
##     ##     motifSel <- res_glmnet%>%filter(comb==ii)%>%pull(motif_name)
##     ##     res0 <- res0%>%filter(motif_name%in%motifSel)
##     ## }
##     ## ###

##     cat(cl_rna, cl_atac, MCl, var0, nrow(res0), "\n")
    

##     df0 <- res_enrich%>%filter(comb==ii, motif_name%in%res0$motif_name)%>%
##         select(motif_name, interest.in, interest.not, not.interest.in, not.interest.not, odds,
##                pval_enrich=pval, padj_enrich=padj, log_odds, log_lower, log_upper)
    
##     df0 <- df0%>%inner_join(res0, by="motif_name")  ##%>%filter(padj_activity<0.1, padj<0.01)
##     df0 <- df0%>%arrange(pval_enrich)
    
##     df0$cluster_RNA <- cl_rna
##     df0$comb <- ii    
##     df0$cluster_atac <- cl_atac
##     df0$celltype <- MCl

##     df0$psycho_variable <- var0
##     df0$variable_name <- var_full
    
##     df0$nDEG <- n_DEG
##     df0$nTFs <- n_TF
##     ## df0$nTFs_olap <- nrow(df0)

##     df_full <- rbind(df_full, df0)
## }
  
## ##
## opfn <- paste(outdir2, "2.1_motif.narrows.txt", sep="")
## write.table(df_full, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")


## ## fn <- paste(outdir2, "2_comb.txt.gz", sep="")
## ## enrich <- fread(fn, header=T, sep="\t")


## outdir2 <- "./3_regulatory.outs/enrich_motif_glmnet/"

## fn <- paste(outdir2, "2.1_motif.narrows.txt", sep="")
## x <- read.table(fn, header=T, sep="\t")

## summ <- x%>%
##     group_by(cluster_RNA, cluster_atac, celltype, psycho_variable, variable_name)%>%
##     summarize(nDEG = mean(nDEG), nTFs = mean(nTFs), nTFs_glmnet = length(unique(motif_name)),
##               nTFs_glmnet_enrich = sum(padj_enrich<0.1), .groups="drop")%>%ungroup()

## summ2 <- summ%>%arrange(psycho_variable, cluster_RNA, cluster_atac)

## ###
## ### save 
## opfn <- paste(outdir2, "2.2_summary.xlsx", sep = "")
## write.xlsx(summ2, file = opfn)

## ###
## ### significant motif

## ## df2 <- df_full%>%filter(padj_activity<0.1, padj<0.01)
## ## opfn <- paste(outdir, "2.1_motif.sig.txt", sep="")
## ## write.table(df2, file=opfn, row.names=F, col.names=T, quote=F, sep="\t")

## ## x2 <- x%>%filter(comb=="C0_PSS_all_mean")

## ###
## ### Forest plots for social support 


## ## col1 <- c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
## ##         "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7" = "#DB72FB")

## fn <- paste(outdir2, "2.1_motif.narrows.txt", sep="")
## plotDF <- fread(fn, header=T, sep="\t", data.table=F)

## plotDF2 <- plotDF%>%filter(comb=="C0_isel")%>%
##     mutate(motif2=gsub("_.*", "", motif_name), direct=ifelse(estimate>0, "gr1", "gr2"))

## p0 <- ggplot(plotDF2, aes(x = log_odds, y = fct_reorder(motif2, log_odds), color=factor(direct)  ))+
##     geom_errorbarh(aes(xmax = log_upper, xmin = log_lower), linewidth = 0.5, height = 0.2)+
##     geom_point(shape = 19, size = 1.5)+
##     scale_color_manual(values=c("gr1"="#F8766D", "gr2"="#00A9FF"),
##         labels=c("gr1"="Positive", "gr2"="Negative"),
##         guide=guide_legend(override.aes=list(size=2)))+    
##     geom_vline(aes(xintercept = 0), linewidth = 0.25, linetype = "dashed")+
##     xlab("log odds ratio")+xlim(-2, 4)+
##     theme_bw()+
##     theme(axis.title.y = element_blank(),
##           axis.title.x = element_text(size=10),
##           axis.text.x = element_text(size=10),
##           axis.text.y = element_text(size=10),
##           legend.title = element_blank(),
##           legend.text = element_text(size=9),
##           legend.key.size = unit(0.35, "cm"))

## ###
## ###
## figfn <- paste(outdir2, "Figure0_isel_C0_enrich.forest.png", sep="")
## ggsave(figfn, p0, width = 620, height = 600, units = "px", dpi=120)



## ###
## ###
## plotDF2 <- plotDF%>%filter(comb=="C4_isel")%>%
##     mutate(motif2=gsub("_.*", "", motif_name), direct=ifelse(estimate>0, "gr1", "gr2"))

## p0 <- ggplot(plotDF2, aes(x = log_odds, y = fct_reorder(motif2, log_odds), color=factor(direct)  ))+
##     geom_errorbarh(aes(xmax = log_upper, xmin = log_lower), linewidth = 0.5, height = 0.2)+
##     geom_point(shape = 19, size = 1.5)+
##     scale_color_manual(values=c("gr1"="#F8766D", "gr2"="#00A9FF"),
##         labels=c("gr1"="Positive", "gr2"="Negative"),
##         guide=guide_legend(override.aes=list(size=2)))+    
##     geom_vline(aes(xintercept = 0), linewidth = 0.25, linetype = "dashed")+
##     xlab("log odds ratio")+xlim(-2, 4)+
##     theme_bw()+
##     theme(axis.title.y = element_blank(),
##           axis.title.x = element_text(size=10),
##           axis.text.x = element_text(size=10),
##           axis.text.y = element_text(size=10),
##           legend.position="none",
##           legend.title = element_blank(),
##           legend.text = element_text(size=9),
##           legend.key.size = unit(0.35, "cm"))

## ###
## ###
## figfn <- paste(outdir2, "Figure0_isel_C4_enrich.forest.png", sep="")
## ggsave(figfn, p0, width = 500, height = 660, units = "px", dpi=120)






###
### END 







    

#####################################################
### filter union of DEGs and union of DARs
####################################################

## rm(list=ls())

## library(limma)
## library(edgeR)


## outdir2 <- "./3_regulatory.outs/TF_regulated/"
## if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings = F, recursive = T)



## ###
## ### 1-Normalize data

## rna_path <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/"
## fn_ls <- list.files(rna_path, "ALL.0.2.11.ComBat_seq")
## fn_ls <- fn_ls[grepl("generem", fn_ls)]

## df_file <- data.frame(file_name=fn_ls, cluster=gsub(".*ComBat_seq\\.|\\.SES.*", "", fn_ls))




## ###
## ### pc for no-related individuals
## path_wk <- "/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/"
## infn_pc <- paste(path_wk, "HOLD-CZI_geno_pc_nokin.txt", sep="")
## pc <- read.table(infn_pc, header=T)
## sampleIDsel <- unique(pc$Sample_ID)


## treat0 <- "CTRL"
## cl_rna_ls <- sort(unique(df_file$cluster))
## for ( cl_rna in cl_rna_ls){

##     cat(cl_rna, "\n")
    
##     ###
##     ### RNA counts data
##     fn0 <- df_file%>%filter(cluster==cl_rna)%>%pull(file_name)
##     fn0 <- paste(rna_path, fn0, sep="")
##     load(fn0)
##     count <- adjusted_counts
    
##     col_RNA <- str_split(colnames(count), "_", simplify=T)
##     col_RNA <- data.frame(comb = colnames(count), cluster = col_RNA[, 1],
##                           sampleID = col_RNA[,3], treat = gsub("RNA-", "", col_RNA[,4]))

##     col_RNA <- col_RNA%>%filter(treat==treat0, sampleID%in%sampleIDsel)

##     count2 <- count[, col_RNA$comb]
##     nsample <- ncol(count2)
      
##     dge <- DGEList(count2)
##     cpm <- cpm(dge)
##     keep.exprs <- rowSums(cpm>0.1) >= 0.2*nsample
##     dge <- dge[keep.exprs, ]
    
##     dge <- calcNormFactors(dge, method = "TMM")
##     ### cpm <- cpm(dge, log=T)
##     ### For plotting and descriptive purposes we recommend cpm(x, log=TRUE) from Gordon Smyth
##     v <- voom(dge, plot=F)
##     geneExpr <- v$E

##     if ( identical(col_RNA$comb, colnames(geneExpr)) ) colnames(geneExpr) <- paste(col_RNA$cluster, col_RNA$sampleID, sep = "_")

##     ###
##     ###
##     opfn <- paste(outdir2, "RNA_", cl_rna, "_normalized.rds", sep="")
##     write_rds(geneExpr, file=opfn)
##     ###
## }    
    



## ###
## ### TF




## ## path_wk <- "/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/"
## ## fn <- paste(path_wk, "HOLD_top9_variables_names.xlsx", sep="")
## ## df_var <- read.xlsx(fn)
## ## varSel <- df_var$Variable



## fn_rna <- "./2_motif.outs/summary_Cluster_res0.07_psycho/Ali_RNA.diff.DESeq.txt.gz"
## res_rna <- fread(fn_rna, header=T, data.table=F) 

## fn <- "./2_motif.outs/summary_Cluster_res0.07_psycho/Ali_RNA.DEGs.txt"
## DEGs <- read.table(fn, header=F)$V1


## fn_motif <- "./2_motif.outs/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz" 
## res_motif <- fread(fn_motif, header=T, data.table=F)


## ## top_motif<- res_motif%>%filter(p.adjusted<0.05)%>%pull(gene)%>%unique()

## ## top_motif <- res_motif%>%filter(gene%in%motif2)%>%
## ##     group_by(psycho_variable)%>%slice_max(order_by=abs(statistic), n=10)%>%pull(gene)%>%unique()

## ###








## ###
## ### covariates
## infn_cv <- paste(path_wk, "HOLD-CZI_covariates_Ali_updated.txt", sep="")
## cv <- read.table(infn_cv, header=T)

## ### Batch
## infn_batch <- paste(path_wk, "HOLD-CZI_covariates_Ali.dbgap.txt", sep="")
## batch_df <- read.table(infn_batch, header=T)%>%dplyr::select(sampleID=dbgap.ID, Batch)

## ###
## ### genotype pcs
## infn_pc <- paste(path_wk, "HOLD-CZI_geno_pc_nokin.txt", sep="")
## pcs <- read.table(infn_pc, header=T)
## pcs <- pcs%>%dplyr::select(sampleID=Sample_ID, PC1, PC2, PC3)



## ###
## ### 1

## fn_anno  <- "../1.2_ArchR_process/4_motif_output/motif_anno/2_comb_motif.annot.txt.gz"
## anno <- fread(fn_anno, header=T, data.table=F)



## ## cluster_df <- read_tsv("./4_find_TF_gene/TableS_summ.tsv", show_col_types = F)
## ## cluster_df <- cluster_df%>%
## ##     mutate(cluster_ATAC=ATAC_cluster, cluster_RNA=paste("C", RNA_cluster, sep=""),
## ##            MCls=gsub(".Cell", "", gsub(".T Cell", "", CellType)))%>%
## ##     dplyr::select(cluster_ATAC, cluster_RNA, MCls)

## ## summ_df <- data.frame(psycho_var=c(rep("isel", each=4), rep("pr_comp", each=3),
## ##                                    rep("PSS_all_mean", each=5), "chronic_sum"),
## ##     cluster_ATAC=c(c("C0", "C1", "C2", "C3"), c("C2", "C3", "C7"),
## ##       c("C0", "C1", "C2", "C3", "C6"), "C3"))

## summ_df <- read.table("./3_regulatory.outs/list_var_clust_nDEGs20_nTF1.txt", header=T, sep="\t")
## summ_df <- summ_df%>%filter(variable%in%c("isel", "pr_comp", "PSS_all_mean", "chronic_sum"), nTFs>1)
## summ_df <- summ_df%>%
##     dplyr::select(cluster_ATAC=Cluster_scATAC, cluster_RNA=Cluster_scRNA, MCls=celltype, variable) 


## motif2 <- sort(unique(res_motif$gene)) ## %>%filter(p.adjusted < 0.1)%>%pull(gene)%>%unique()
## res_all <- NULL
## for ( i in 1:nrow(summ_df)){
##     ###
##     ###
##     var0 <- summ_df$variable[i] 
##     cl_rna <- summ_df$cluster_RNA[i]
##     oneMCl <- summ_df$MCls[i]

##     fn <- paste(outdir2, "RNA_", cl_rna, "_normalized.rds", sep="") 
##     rna_norm <- read_rds(fn)
##     colnames(rna_norm) <- gsub("C._", "", colnames(rna_norm))
    
    
##     cv2 <- cv%>%dplyr::select(sampleID=dbgap.ID, age, sex_alph, z=all_of(var0))
##     cv2 <- cv2%>%inner_join(batch_df, by="sampleID")%>%
##         inner_join(pcs, by="sampleID")

##     cat(cl_rna, oneMCl, var0, "\n")

##     for ( motif0 in motif2){
##         ##
##         gene0 <- anno%>%filter(motif_name==motif0)%>%pull(gene)%>%unique()
##         gene2 <- intersect(intersect(DEGs, gene0), rownames(rna_norm))

##         y2 <- colMeans(rna_norm[gene2, ])
    
##         cv2i <- cv2%>%mutate(y=y2[sampleID])
##         cv2i <- cv2i%>%drop_na(all_of(colnames(cv2i)))

##         lm0 <- try(lm(y ~ 0 + factor(sex_alph) + factor(Batch) + PC1 + PC2 + age + z, data=cv2i), silent=T)

##         if ( class(lm0) == "try-error") next

##         bhat <- coef(lm0)["z"]
##         sdhat <- sqrt(diag(vcov(lm0))["z"])
##         zval <- bhat/sdhat
##         p <- 2*pnorm(-abs(zval))
##         ###

##         res0 <- data.frame(cluster_RNA=cl_rna, MCls=oneMCl, psycho_var=var0, motif_name=motif0,
##                estimate = bhat, stderror = sdhat, zscore = zval, p.value = p)            
##         res_all <- rbind(res_all, res0)
##      }   
##  }      
    
       
## opfn <- paste(outdir2, "1_RNA_TF.lm.txt", sep="")
## write.table(res_all, file=opfn, sep="\t", row.names=F, quote=F)

    
    
    
    

## for ( i in 1:nrow(
## cl_rna <- "C0"
## cl_atac <- "C0"
## oneMCl <- "CD4+"

## for ( var0 in c("PSS_all_mean", "isel")){
    
## cat(var0, "\n")


## res0 <- res_motif%>%filter(gene%in%top_motif, psycho_variable==var0, Cluster==cl_atac)%>%
##    dplyr::select(Cluster, psycho_variable, motif_name=gene, x_estimate_motif=estimate, x_zscore_motif=statistic) 

## res0_rna <- res_rna%>%filter(var==var0, cluster==cl_rna)%>%mutate(gene=identifier, zscore=logFC/SE)

## df2 <- NULL
## for ( ii in unique(res0$motif_name)){
##     ##
##     fn <- paste(prefix, "2_motif_", ii, ".annotation.txt", sep="")
##     x <- read.table(fn, header=T, sep="\t")
##     x <- x%>%filter(gene%in%DEGs)
##     geneSel <- unique(x$gene)

##     cat(ii, length(geneSel), "\n")
    
##     y_LFC <- res0_rna%>%filter(gene%in%geneSel)%>%pull(logFC)%>%median(na.rm=T)
##     y_zscore <- res0_rna%>%filter(gene%in%geneSel)%>%pull(zscore)%>%median(na.rm=T)

##     df0 <- data.frame(motif_name=ii, y_estimate_gene=y_LFC, y_zscore_gene=y_zscore)
##     df2 <- rbind(df2, df0)
## }    


## plotData <- res0%>%left_join(df2, by="motif_name")
fn <- paste(outdir2, "1_RNA_TF.lm.txt", sep="")
res_rna <- read.table(fn, header=T, sep="\t")

##
fn_motif <- "./2_motif.outs/summary_Cluster_res0.07_psycho/2_th20_plotData.comb.txt.gz" 
res_motif <- fread(fn_motif, header=T, data.table=F)
motif2 <- res_motif%>%dplyr::filter(p.adjusted<0.1)%>%pull(gene)%>%unique()

summ_df <- read.table("./3_regulatory.outs/list_var_clust_nDEGs20_nTF1.txt", header=T, sep="\t")
summ_df <- summ_df%>%filter(variable%in%c("isel", "pr_comp", "PSS_all_mean", "chronic_sum"), nTFs>1)
summ_df <- summ_df%>%
    dplyr::select(cluster_ATAC=Cluster_scATAC, cluster_RNA=Cluster_scRNA, MCls=celltype, variable)


###
### plot
varSel <- sort(unique(summ_df$variable))

var0 <- varSel[1]
summ0 <- summ_df%>%filter(variable==var0)

plotDF <- NULL
for ( i in 1:nrow(summ)){
   ##
   cl_atac <- summ0$cluster_ATAC[i]
   cl_rna <- summ0$cluster_RNA[i]
   res0 <- res_rna%>%filter(cluster_RNA==cl_rna, psycho_var==var0, motif_name%in%motif2)%>%
       select(cluster_RNA, MCls, psycho_var, motif_name, y_zscore_gene=zscore)

   res2 <- res_motif%>%filter(Cluster==cl_atac, psycho_variable==var0, gene%in%motif2)%>%
       select(cluster_ATAC=Cluster, motif_name=gene, x_zscore_motif=statistic)

   res0 <- res0%>%inner_join(res2, by="motif_name") 
   plotDF <- rbind(plotDF, res0)
}    


###
###
p0 <- ggplot(plotDF, aes(x=x_zscore_motif, y=y_zscore_gene))+
   geom_point(linewidth=1.5, color="#F87660")+
   scale_x_continuous(bquote(~italic(Z)~"score from TF activity"), expand=expansion(mult=0.1))+
   scale_y_continuous(bquote(~italic(Z)~"score from TF regulated genes"), expand=expansion(mult=0.1))+
   facet_wrap(~MCls, ncol=3, scales="free")+ 
   ggtitle(var0)+
   theme_bw()+
   theme(## axis.title = element_text(size=10),
         ## axis.text = element_text(size=10),
         plot.title = element_text(hjust=0.5),
         strip.text=element_text(size=12))

p0 <- p0 +
    geom_smooth(data=plotDF, aes(x=x_zscore_motif, y=y_zscore_gene),
    method="lm", formula=y~x, color="blue", linewidth=0.5, se=F)          



##
### correlation
corr <- cor.test(plotDF$x_zscore_motif, plotDF$y_zscore_gene)
r <- round(as.numeric(corr$estimate), digits=3)
pval <- corr$p.value

symb <- "NS"
if (pval<0.05) symb <- "*"
if (pval<0.01) symb <- "**"
if (pval<0.001) symb <- "***"
## annot_expression <- "paste(italic(R), \"=\", .(r),  .(symb))"
##
annot_express <-deparse( bquote(italic(R)==~.(r)~.(symb)))

###
xmin <- min(plotDF$x_zscore_motif)
xmax <- max(plotDF$x_zscore_motif)
xpos <- xmin+0.13*(xmax-xmin)
###
ymin <- min(plotDF$y_zscore_gene)
ymax <- max(plotDF$y_zscore_gene)
ypos <- ymax-0.03*(ymax-ymin)

## eq <- bquote(italic(R)==~"0.814,"~"**
p0 <- p0 + annotate("text", x = xpos, y = ypos, label = annot_express, parse = T, color = "blue")



##
figfn <- paste(outdir2, "Figure_", var0, ".scatter.png", sep="")
ggsave(figfn, p0, width=480, height=480, units="px", dpi=120)




###
###
feq <- function(x){
   ## 
   r <- round(as.numeric(x$estimate),digits=3)
   p <- x$p.value
   if(p<0.001) symb <- "***"
   if(p>=0.001 & p<0.01) symb <- "**"
   if (p>=0.01 & p<0.05) symb <- "*"
   if(p>0.05) symb <- "NS"

    eq <- bquote(italic(R)==.(r)~","~.(symb))
    eq
}

##
xFun <- function(dx,a=0.2){
   min1 <- min(dx$x_zscore_motif, na.rm=T)
   max1 <- max(dx$x_zscore_motif, na.rm=T)
   R <- max1-min1
   xpos <- min1+a*R
}
##
yFun <- function(dx,a=0.9){
  min1 <- min(dx$y_zscore_gene, na.rm=T)
  max1 <- max(dx$y_zscore_gene, na.rm=T)
  R <- max1-min1
  ypos <- min1+a*R
}





var0 <- varSel[2]
summ0 <- summ_df%>%filter(variable==var0)

motif2 <- res_motif%>%filter(psycho_variable==var0, p.adjusted<0.1)%>%pull(gene)%>%unique()
 
plotDF <- NULL
for ( i in 1:nrow(summ)){
   ##
   cl_atac <- summ0$cluster_ATAC[i]
   cl_rna <- summ0$cluster_RNA[i]
   res0 <- res_rna%>%filter(cluster_RNA==cl_rna, psycho_var==var0, motif_name%in%motif2)%>%
       select(cluster_RNA, MCls, psycho_var, motif_name, y_zscore_gene=zscore)

   res2 <- res_motif%>%filter(Cluster==cl_atac, psycho_variable==var0, gene%in%motif2)%>%
       select(cluster_ATAC=Cluster, motif_name=gene, x_zscore_motif=statistic)

   res0 <- res0%>%inner_join(res2, by="motif_name") 
   plotDF <- rbind(plotDF, res0)
}    


###
###
p0 <- ggplot(plotDF, aes(x=x_zscore_motif, y=y_zscore_gene))+
   geom_point(size=1.5, color="#F87660")+
   scale_x_continuous(bquote(~italic(Z)~"score from TF activity"), expand=expansion(mult=0.1))+
   scale_y_continuous(bquote(~italic(Z)~"score from TF regulated genes"), expand=expansion(mult=0.1))+
   facet_wrap(~MCls, ncol=2, scales="free")+  
   ggtitle(var0)+ 
   theme_bw()+
   theme(## axis.title = element_text(size=10),
         ## axis.text = element_text(size=10),
         plot.title = element_text(hjust=0.5))

p0 <- p0 +
    geom_smooth(data=plotDF, aes(x=x_zscore_motif, y=y_zscore_gene),
    method="lm", formula=y~x, color="blue", linewidth=0.5, se=F)          



##
### correlation



anno_df2 <- plotDF%>%group_by(MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$x_zscore_motif, (.x)$y_zscore_gene, method="pearson")),
           eq=map(corr, feq),
           r2=map_dbl(corr, ~(.x)$estimate),
           xpos=map_dbl(data, xFun),
           ypos=map_dbl(data, yFun))

           
p0 <- p0 + geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)


## 
figfn <- paste(outdir2, "Figure_", var0, ".scatter.png", sep="")
ggsave(figfn, p0, width=580, height=580, units="px", dpi=120)




####
####

var0 <- varSel[3]
summ0 <- summ_df%>%filter(variable==var0)

motif2 <- res_motif%>%filter(psycho_variable==var0, p.adjusted<0.1)%>%pull(gene)%>%unique()
 
plotDF <- NULL
for ( i in 1:nrow(summ)){
   ##
   cl_atac <- summ0$cluster_ATAC[i]
   cl_rna <- summ0$cluster_RNA[i]
   res0 <- res_rna%>%filter(cluster_RNA==cl_rna, psycho_var==var0, motif_name%in%motif2)%>%
       select(cluster_RNA, MCls, psycho_var, motif_name, y_zscore_gene=zscore)

   res2 <- res_motif%>%filter(Cluster==cl_atac, psycho_variable==var0, gene%in%motif2)%>%
       select(cluster_ATAC=Cluster, motif_name=gene, x_zscore_motif=statistic)

   res0 <- res0%>%inner_join(res2, by="motif_name") 
   plotDF <- rbind(plotDF, res0)
}    


###
###
p0 <- ggplot(plotDF, aes(x=x_zscore_motif, y=y_zscore_gene))+
   geom_point(size=1.5, color="#F87660")+
   scale_x_continuous(bquote(~italic(Z)~"score from TF activity"), expand=expansion(mult=0.1))+
   scale_y_continuous(bquote(~italic(Z)~"score from TF regulated genes"), expand=expansion(mult=0.1))+
   facet_wrap(~MCls, ncol=3, scales="free")+  
   ggtitle(var0)+ 
   theme_bw()+
   theme(## axis.title = element_text(size=10),
         ## axis.text = element_text(size=10),
         plot.title = element_text(hjust=0.5))

p0 <- p0 +
    geom_smooth(data=plotDF, aes(x=x_zscore_motif, y=y_zscore_gene),
    method="lm", formula=y~x, color="blue", linewidth=0.5, se=F)          



##
### correlation



anno_df2 <- plotDF%>%group_by(MCls)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$x_zscore_motif, (.x)$y_zscore_gene, method="pearson")),
           eq=map(corr, feq),
           r2=map_dbl(corr, ~(.x)$estimate),
           xpos=map_dbl(data, xFun),
           ypos=map_dbl(data, yFun))

           
p0 <- p0 + geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)


## 
figfn <- paste(outdir2, "Figure_", var0, ".scatter.png", sep="")
ggsave(figfn, p0, width=980, height=380, units="px", dpi=120)




###
###
var0 <- varSel[4]
summ0 <- summ_df%>%filter(variable==var0)

motif2 <- res_motif%>%filter(psycho_variable==var0, p.adjusted<0.1)%>%pull(gene)%>%unique()
 
plotDF <- NULL
for ( i in 1:nrow(summ)){
   ##
   cl_atac <- summ0$cluster_ATAC[i]
   cl_rna <- summ0$cluster_RNA[i]
   res0 <- res_rna%>%filter(cluster_RNA==cl_rna, psycho_var==var0, motif_name%in%motif2)%>%
       select(cluster_RNA, MCls, psycho_var, motif_name, y_zscore_gene=zscore)

   res2 <- res_motif%>%filter(Cluster==cl_atac, psycho_variable==var0, gene%in%motif2)%>%
       select(cluster_ATAC=Cluster, motif_name=gene, x_zscore_motif=statistic)

   res0 <- res0%>%inner_join(res2, by="motif_name") 
   plotDF <- rbind(plotDF, res0)
}    

plotDF <- plotDF%>%mutate(MCl2=paste(cluster_RNA, MCls, sep="-"))

###
###
p0 <- ggplot(plotDF, aes(x=x_zscore_motif, y=y_zscore_gene))+
   geom_point(size=1.5, color="#F87660")+
   scale_x_continuous(bquote(~italic(Z)~"score from TF activity"), expand=expansion(mult=0.1))+
   scale_y_continuous(bquote(~italic(Z)~"score from TF regulated genes"), expand=expansion(mult=0.1))+
   facet_wrap(~MCl2, ncol=3, scales="free")+  
   ggtitle(var0)+ 
   theme_bw()+
   theme(## axis.title = element_text(size=10),
         ## axis.text = element_text(size=10),
         plot.title = element_text(hjust=0.5))

p0 <- p0 +
    geom_smooth(data=plotDF, aes(x=x_zscore_motif, y=y_zscore_gene),
    method="lm", formula=y~x, color="blue", linewidth=0.5, se=F)          



##
### correlation



anno_df2 <- plotDF%>%group_by(MCl2)%>%
    nest()%>%
    mutate(corr=map(data, ~cor.test((.x)$x_zscore_motif, (.x)$y_zscore_gene, method="pearson")),
           eq=map(corr, feq),
           r2=map_dbl(corr, ~(.x)$estimate),
           xpos=map_dbl(data, xFun),
           ypos=map_dbl(data, yFun))

           
p0 <- p0 + geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), colour="blue", size=3, parse=T)


## 
figfn <- paste(outdir2, "Figure_", var0, ".scatter.png", sep="")
ggsave(figfn, p0, width=980, height=600, units="px", dpi=120)





