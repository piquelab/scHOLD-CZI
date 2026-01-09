####
####
library(tidyverse)
library(Matrix)
library(DESeq2)
library(biobroom)
library(data.table)
library("SummarizedExperiment")
##
library(cowplot)
library(RColorBrewer)
library(scales)
library(ComplexHeatmap)
library(openxlsx)
library(ggrastr)
 
###

rm(list=ls())

option <- "option_nFeature15K_cluster_res0.12"
type <- "psycho"
outdir2 <- paste("./1_DiffPeak.outs/summary_", option, "_", type,  "_CTRL/", sep="")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


###
###
outdir <- paste("./1_DiffPeak.outs/", option, "_th20_psycho/", sep="")
fn_ls <- list.files(outdir, "_DESeq.summ.tsv$")
fn_ls <- fn_ls[grepl(type, fn_ls)]
fn_ls <- fn_ls[grepl("CTRL", fn_ls)]
  
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
fn <- paste(outdir2, "1_", type, "_summary.sigs.tsv", sep="")
summ_df <- read_tsv(fn, show_col_types=F)%>%as.data.frame()

###
### variable
varSel <- sort(read.table("./psycho_var_dir/psychosocial_top10.txt")$V1)
varSel2 <- varSel[c(3:4, 7:10)]


###
### summarize sigs
summ2 <- summ_df%>%filter(variable%in%varSel2)%>%arrange(variable, .locale="en")%>%
    pivot_wider(id_cols=variable, names_from="Cluster", values_from=nsig_fdr0.1, values_fill=0)
###
opfn2 <- paste(outdir2, "1.2_", type, "_mat.sigs.xlsx", sep="")
write.xlsx(summ2, file=opfn2, overwrite=T)
    


###
### summarize #inds

## summ2 <- summ_df%>%filter(variable%in%varSel2)%>%arrange(variable, .locale="en")%>%
##     dplyr::select(variable, Cluster, nind, ngene=ngene_ind)
## summ2 <- summ2%>%    
##     pivot_wider(id_cols=variable, names_from="Cluster", values_from=c(nind, ngene), values_fill=0,
##                 names_glue="{Cluster}_{.value}")%>%as.data.frame()
## ###
## colSel <- c("variable", sort(colnames(summ2)[-1]))
## summ2 <- summ2[, colSel]
## ###
## opfn2 <- paste(outdir2, "1.2_th20_mat.xlsx", sep="")
## write.xlsx(summ2, file=opfn2, overwrite=T)





###########################
### summary cytokines
###############################

type <- "psycho"
fn <- paste(outdir2, "1_", type, "_summary.sigs.tsv", sep="")
summ_df <- read_tsv(fn, show_col_types=F)%>%as.data.frame()
 
## ###
## ### variable
## varSel <- sort(read.table("./psycho_var_dir/cytokines_var11.txt")$V1)


## ###
## ### summarize sigs
## summ2 <- summ_df%>%filter(variable%in%varSel)%>%arrange(variable, .locale="en")%>%
##     pivot_wider(id_cols=variable, names_from="Cluster", values_from=nsig_fdr0.1, values_fill=0)

## ###
## opfn2 <- paste(outdir2, "1.2_cytokines_mat.sigs.xlsx", sep="")
## write.xlsx(summ2, file=opfn2, overwrite=T)



###
#### full table

summ2 <- summ_df%>%filter(variable%in%varSel2)%>%arrange(variable, .locale="en")%>%
    dplyr::select(variable, Cluster, nind, ngene=ngene_ind, nsig=nsig_fdr0.1)
summ2 <- summ2%>%
    pivot_wider(id_cols=variable, names_from="Cluster", values_from=c(nind, ngene, nsig), values_fill=0,
                names_glue="{Cluster}_{.value}")%>%as.data.frame()
###
colSel <- c("variable", sort(colnames(summ2)[-1]))
summ2 <- summ2[, colSel]
###
opfn2 <- paste(outdir2, "1.2_", type, "_mat.full.xlsx", sep="")
write.xlsx(summ2, file=opfn2, overwrite=T)





###
### summarize #test genes
## summ2 <- summ_df%>%filter(variable%in%varSel2)%>%arrange(variable, .locale="en")%>%
##     pivot_wider(id_cols=variable, names_from="Cluster", values_from=ngene_test, values_fill=0)
## ###
## opfn2 <- paste(outdir2, "1.2_th", th0, "_mat.ngene.xlsx", sep="")
## write.xlsx(summ2, file=opfn2, overwrite=T)


    
## df_var <- read.table("./psycho_var_dir/psychosocial_variables", header=F)%>%
##     filter(!grepl("^PF|^Log_PF|^LPF|^Log_LPF", V1))
## write.table(df_var$V1, file="./psycho_var_dir/psychosocial_variables_interest.txt",
##             quote=F, row.names=F, col.names=F)

####################################################
### summarize psychosocial variables of interest
####################################################

## path_wk <- "/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/"
## th0 <- 20

## ###
## ### 10 psychosocial variable (top 9 from Ali and logil6) 
## fn <- paste(path_wk, "HOLD_top9_variables_names.xlsx", sep="")
## df_psycho <- read.xlsx(fn)
## names(df_psycho)[1:2] <- c("variable", "variable.name")

## d0 <- data.frame("variable"="Logil6", "variable.name"="Logil6", "category"="inflammation")
## df_psycho <- rbind(df_psycho, d0)

## ## opfn <- "./psycho_var_dir/psychosocial_top10.txt"
## ## write.table(sort(df_psycho$variable), file=opfn, row.names=F, quote=F, col.names=F)

## ###
## ### 
## fn <- paste(outdir2, "1_th", th0, "_summary.sigs.tsv", sep="")
## summ_df <- read_tsv(fn)%>%as.data.frame()

## mat_df <- summ_df%>%dplyr::rename(nsig=nsig_fdr0.1, ntest=ngene_test)%>%
##    pivot_wider(id_cols=variable, names_from=Cluster,
##    values_from=c(nsig, nind), values_fill=0,
##    names_glue="{Cluster}_{.value}")

 
## colSel <- paste(paste("C", rep(0:6, each=2), sep=""), "_", rep(c("nsig", "nind"), times=7), sep="")

## mat_val <- mat_df%>%dplyr::select(all_of(c("variable", colSel)))%>%as.data.frame()
 
## df2 <- mat_df%>%dplyr::select(variable)%>%inner_join(df_psycho, by="variable")%>%
##     left_join(mat_val, by="variable")

## ###
## opfn <- paste(outdir2, "1.2_th", th0, "_top10vars.mat.xlsx", sep="")
## write.xlsx(df2, file=opfn)

    




##########################
### QQ plots 
########################




## fn <- paste(outdir2, "1.2_th", th0, "_mat.sigs.xlsx", sep="")
## summ <- read.xlsx(fn)
## rnz <- rowSums(summ[,-1])
## names(rnz) <- summ$variable
## rnz <- sort(rnz, decreasing=T)

## ###
## varSel2 <- union(varSel, names(rnz[rnz>100]))


option <- "option_nFeature15K_cluster_res0.12"
type <- "psycho"
outdir2 <- paste("./1_DiffPeak.outs/summary_", option, "_", type,  "_CTRL/", sep="")
outdir <- paste("./1_DiffPeak.outs/", option, "_th20_psycho/", sep="")



###
### variable
varSel <- sort(read.table("./psycho_var_dir/psychosocial_top10.txt")$V1)
varSel2 <- varSel[c(3:4, 7:10)]

##
fn_ls <- list.files(outdir, "_DESeq.results.txt.gz$")
fn_ls <- fn_ls[grepl("CTRL", fn_ls)]
df_fn <- data.frame(fname=fn_ls,
     psycho_var = gsub("_DESeq.results.txt.gz", "", gsub(".*_CTRL_", "", fn_ls)),
     cluster = gsub("_.*", "", fn_ls))
df_fn <- df_fn%>%filter(psycho_var%in%varSel2)

res_all <- NULL
for ( ii in df_fn$fname){
    ###
    ### 
    fn0 <- paste(outdir, ii, sep="")
    res <- fread(fn0, header=T, data.table=F)
    cat(ii, res$nind[1], "\n") 

    ntest <- nrow(res)
    res <- res%>%
        arrange(p.value)%>%
        mutate(observed=-log10(p.value), expected=-log10(ppoints(ntest)))
    ###
    ### 
    res_all <- rbind(res_all, res)     
}

opfn <- paste(outdir2, "2_", type, "_plotData.comb.txt.gz", sep="")
fwrite(res_all, file=opfn, sep="\t", quote=F, na=NA)


###
### plots
 
col_cl <- hue_pal()(11)
names(col_cl) <- paste("C", 0:10, sep="")

infn <- paste(outdir2, "2_", type, "_plotData.comb.txt.gz", sep="")
plotDF <- fread(infn, header=T, data.table=F)


###
###
p0 <- ggplot(plotDF, aes(x=expected, y=observed, color=Cluster))+


    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=col_cl[1:8],
        guide=guide_legend(override.aes=list(size=3)))+
    facet_wrap(~psycho_variable, scales="free_y", nrow=2)+  ## ncol=7       
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=10),
          legend.key.size=grid::unit(0.4, "cm"),
          ## legend.position="inside",   
          ## legend.position.inside=c(0.9, 0.3),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          strip.text=element_text(size=12))

###
###
figfn <- paste(outdir2, "Figure2.1_", type, "_ctrl.qq.pdf", sep="")
ggsave(figfn, p0, width=9, height=5)
 





#####################################
### QQ plots cytokines, PFAS
#####################################

                                        
varSel <- sort(read.table("./psycho_var_dir/PFAS_vars.txt")$V1)
## varSel2 <- varSel[c(3:4, 6:10)]


outdir <- paste("./1_DiffPeak.outs/", option, "_th20_psycho_CTRL/", sep="")
fn_ls <- list.files(outdir, "_DESeq.results.txt.gz$")


###
### 
res_all <- NULL
for ( ii in fn_ls){
    ###
    psycho_var <- gsub(".*_CTRL_|_combat_seq.*", "", ii)

    if ( psycho_var%in%varSel){
       ### 
       fn0 <- paste(outdir, ii, sep="")
       res <- fread(fn0, header=T, data.table=F)
       cat(psycho_var, res$nind[1], "\n") 

       ntest <- nrow(res)
       res <- res%>%
           arrange(p.value)%>%
           mutate(observed=-log10(p.value), expected=-log10(ppoints(ntest)))
       ###
       ### 
       res_all <- rbind(res_all, res)
   }     
}

opfn <- paste(outdir2, "2_PFAS_plotData.comb.txt.gz", sep="")
fwrite(res_all, file=opfn, sep="\t", quote=F, na=NA)


###
### plots
## hue_pal()(9)
 
infn <- paste(outdir2, "2_PFAS_plotData.comb.txt.gz", sep="")
plotDF <- fread(infn, header=T, data.table=F)


###
###
var0 <- sort(unique(plotDF$psycho_variable))
i <- 2
i0 <- (i-1)*12+1
i1 <- min(i*12, length(var0))
var2 <- var0[i0:i1]
p0 <- ggplot(plotDF%>%filter(psycho_variable%in%var2), aes(x=expected, y=observed, color=Cluster))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
        "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB"),
        guide=guide_legend(override.aes=list(size=2)))+
    facet_wrap(~psycho_variable, scales="free_y", ncol=6)+  ## ncol=7       
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=10),
          legend.key.size=grid::unit(0.3, "cm"),
          legend.position="inside",   
          legend.position.inside=c(0.9, 0.12),
          ###legend.position="none",
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          strip.text=element_text(size=10))
 
###
###
figfn <- paste(outdir2, "Figure2.1_PFAS", i, "_ctrl.qq.png", sep="")
ggsave(figfn, p0, width=900, height=420, units="px", dpi=100)






#########################################
### scatter plots isel and ISEL_Mean
###########################################


infn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
df0 <- fread(infn, header=T, data.table=F)

###
plotDF <- map_dfr(sort(unique(df0$Cluster)), function(ii){
    ###
    df_old <- df0%>%filter(Cluster==ii, psycho_variable=="isel")%>%
        dplyr::select(Cluster, gene, beta_x=estimate, zscore_x=statistic,  padj_x=p.adjusted)
    df_new <- df0%>%filter(Cluster==ii, psycho_variable=="ISEL_Mean")%>%
        dplyr::select(gene, beta_y=estimate, zscore_y=statistic, padj_y=p.adjusted)

    df2 <- df_old%>%left_join(df_new, by="gene")
    df2
})    



###
###
plotDF2 <- plotDF%>%
    mutate(sig_gr=case_when(padj_x<0.1&padj_y>0.1 ~ "sig1",
            padj_x>0.1&padj_y<0.1 ~ "sig2",
            padj_x<0.1&padj_y<0.1 ~ "sig3",
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
   group_by(Cluster)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$beta_x, (.x)$beta_y, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=mypos(plotDF2$beta_x, 0.25),
          ypos=mypos(plotDF2$beta_x, 0.9))%>%
       dplyr::select(-data,-corr)

 
       

###
### plots 
pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=beta_x, y=beta_y))+ 
    geom_point(aes(color=sig_gr), size=0.8)+
    scale_color_manual(values=c("sig0"="grey", "sig1"="red", "sig2"="blue", "sig3"="green"), 
        labels=c("sig0" = "Not", "sig1"="DAM-old", "sig2"="DAM-new", "sig3"="Both"),
        breaks=c("sig0", "sig1", "sig2", "sig3"),
        guide=guide_legend(override.aes=list(size=2)))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    facet_wrap(~Cluster, scales="fixed", nrow=3)+
    scale_x_continuous(bquote(beta~"on DAR in isel"), expand=expansion(mult=0.1))+
    scale_y_continuous(bquote(beta~"on DAR in SEL_Mean"), expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_blank(),
          legend.key.size=grid::unit(0.4, "cm"))
 
figfn <- paste(outdir2, "Figure0.1_beta_isel.scatter.png", sep="")
ggsave(figfn, pcomb, width=900, height=800, units="px", dpi=120)


###
### compare z-score 
###


### annotation text 
anno_df2 <- plotDF2%>%
   group_by(Cluster)%>%
   nest()%>%
   mutate(corr=map(data, ~cor.test((.x)$zscore_x, (.x)$zscore_y, method="pearson")),
          eq=map(corr,feq),
          rr=map_dbl(corr, ~(.x)$estimate),
          xpos=mypos(plotDF2$zscore_x, 0.25),
          ypos=mypos(plotDF2$zscore_x, 0.9))%>%
       dplyr::select(-data,-corr)
        

###
### plots 
pcomb <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=zscore_x, y=zscore_y))+ 
    geom_point(aes(color=sig_gr), size=0.8)+
    scale_color_manual(values=c("sig0"="grey", "sig1"="red", "sig2"="blue", "sig3"="green"), 
        labels=c("sig0" = "Not", "sig1"="DAM-old", "sig2"="DAM-new", "sig3"="Both"),
        breaks=c("sig0", "sig1", "sig2", "sig3"),
        guide=guide_legend(override.aes=list(size=2)))+
    geom_text(data=anno_df2, aes(x=xpos, y=ypos, label=eq), size=3, parse=T)+
    geom_abline(slope=1, intercept=0, linewidth=0.3)+
    facet_wrap(~Cluster, scales="fixed", nrow=3)+
    scale_x_continuous(bquote(italic(Z)~"-score"~"on DAR in isel"), expand=expansion(mult=0.1))+
    scale_y_continuous(bquote(italic(Z)~"-score"~"on DAR in SEL_Mean"), expand=expansion(mult=0.1))+
    theme_bw()+
    theme(strip.text=element_text(size=12),
          axis.title=element_text(size=10),
          axis.text=element_text(size=10),
          legend.title=element_blank(),
          legend.key.size=grid::unit(0.4, "cm"))
 
figfn <- paste(outdir2, "Figure0.1_zscore_isel.scatter.png", sep="")
ggsave(figfn, pcomb, width=900, height=800, units="px", dpi=120)







#######################################################################################
### stratified qq plots and annotate peaks into two groups, in DAM and not DAMs
#######################################################################################

## ###
## outdir2 <- "./1_DiffPeak.outs/summary_Cluster_res0.07_psycho_CTRL/"

## ###
## ### variables 
## ## path_wk <- "/tier2/home/groups/pique20/schold/ATAC/sc-atac-cziHOLD/"
## ## th0 <- 20

## ## fn <- paste(path_wk, "HOLD_top9_variables_names.xlsx", sep="")
## ## df_psycho <- read.xlsx(fn)

## ## varSel2 <- c(df_psycho$Variable, "Logil6")


## #########################
## ### motif match Peak
## #########################

## fn_motif <- "../1.2_ArchR_process/4_motif_output/Motif_jaspar2022.match.mat.rds"
## motif_match <- read_rds(fn_motif)

## df_row <- rowRanges(motif_match)
## anno <- data.frame(chr=as.character(gsub("^chr", "", seqnames(df_row))),
##     peak=paste(seqnames(df_row), start(df_row), end(df_row), sep="_"))

## mat <- assay(motif_match)
## rownames(mat) <- anno$peak

## peak2 <- anno%>%dplyr::filter(chr%in%as.character(1:22))%>%pull(peak)

## mat2 <- mat[peak2,]
## ## rnz <- rowSums(mat2)
 

## ###################
## ### get DAMs
## ###################

## fn_DAM <- "./3_regulatory.outs/correct_excludeX/TableS_motif_all.infor.txt"
## res_DAM <- read.table(fn_DAM, header=T, sep="\t")
## res_DAM <- res_DAM%>%mutate(comb=paste(cluster_atac, psycho_variable, sep="_"))
 
## comb_sel <- sort(unique(res_DAM$comb))


## ###
## ### diff data
## infn <- paste(outdir2, "2_th20_plotData.comb.txt.gz", sep="")
## res <- fread(infn, header=T, data.table=F)
## res <- res%>%mutate(comb=paste(Cluster, psycho_variable, sep="_"))%>%
##     dplyr::filter(comb%in%comb_sel)


## res_new <- NULL
## summ_all <- NULL
## for ( ii in comb_sel){
##     ###
        
##     motif2 <- res_DAM%>%filter(comb==ii)%>%pull(motif_name)%>%unique()

##     mat0 <- mat2[, motif2, drop=F]
##     rnz <- rowSums(mat0)
##     peak2 <- names(rnz)[rnz>0]
##     ##
##     res0 <- res%>%filter(comb==ii)%>%mutate(is_motif=ifelse(gene%in%peak2, "gr1", "gr0"))   
##     ## res0 <- res0%>%group_by(is_motif)%>%mutate(padj2=p.adjust(p.value, "BH"))%>%ungroup()

##     gr <- unique(res0$is_motif)
##     res0 <- map_dfr(gr, function(ii){
##         ##
##         res00 <- res0%>%filter(is_motif==ii)
##         th0 <- res0%>%drop_na(p.adjusted)%>%pull(baseMean)%>%min()

##         ###
##         iwhich <- which(res00$baseMean>th0)
##         pval <- res00[iwhich, "p.value"]
##         padj_BH <- p.adjust(pval, method="BH")
        
##         res00$padj2 <- NA
##         res00$padj2[iwhich] <- padj_BH
        
##         ###
##         ###
##         ntest <- nrow(res00)
##         res00 <- res00%>%arrange(p.value)%>%
##             mutate(observed_gr=-log10(p.value), expected_gr=-log10(ppoints(ntest)))         
##         res00
##     })

##     summ0 <- data.frame(comb=ii, nsig_0=sum(res0$p.adjusted<0.1, na.rm=T), nsig_str=sum(res0$padj2<0.1, na.rm=T),
##        nmotif = length(motif2), npeak_gr1=sum(res0$is_motif=="gr1"), npeak_gr0=sum(res0$is_motif=="gr0"))         
##     summ_all <- rbind(summ_all, summ0)                    

##     ##
##     res_new <- rbind(res_new, res0)
##     cat(ii, sum(res0$is_motif=="gr1"), "\n")
## }
 

## ###
## ### save
## opfn <- paste(outdir2, "3_plotData_gr.comb.txt.gz", sep="")
## fwrite(res_new, file=opfn, sep="\t", na=NA)

## ###
## opfn <- paste(outdir2, "3_summary.tsv", sep="")
## write_tsv(summ_all, file=opfn)



## ###
## ###


    

## ######################
## ### qq plots
## ######################

## fn <- paste(outdir2, "3_th20_plotData_gr.comb.txt.gz", sep="")
## plotDF <- fread(fn, header=T, data.table=F)

## plotDF <- res_new
## p0 <- ggplot(plotDF, aes(x=expected_gr, y=observed_gr, color=factor(is_motif)))+
##     rasterise(geom_point(size=0.06), dpi=300)+
##     geom_abline(color="grey")+
##     scale_color_manual(values=c("gr0"="#f1b6da", "gr1"="#d01c8b"),
##         labels=c("gr0"="Not in DAMs", "gr1"="In DAMs"),
##         guide=guide_legend(override.aes=list(size=2)))+
##   facet_wrap(~comb, scales="free_x", ncol=3)+  ## ncol=7       
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
## figfn <- paste(outdir2, "Figure3.1_ctrl.stratified.qq.png", sep="")
## ggsave(figfn, p0, width=900, height=800, units="px", dpi=120)




###################################################
### peak annotation and stratified qq plots 
####################################################


option <- "option_nFeature15K_cluster_res0.12"
type <- "psycho"
outdir2 <- paste("./1_DiffPeak.outs/summary_", option, "_", type,  "_CTRL/", sep="")
if ( !file.exists(outdir2)) dir.create(outdir2, showWarnings=F, recursive=T)


anno_fn <- "../1.2_ArchR_process/4b_motif_output/motif_anno/1_peak.annotation.rds"
anno <- read_rds(anno_fn)
anno <- anno%>%filter(abs(distToGeneStart)<1e+05)

###
### gene expression
fn_rna <- paste("./2_motif.outs/summary_", option, "_psycho_CTRL/Ali_RNA_psycho.diff.DESeq.txt.gz", sep="")
res_rna <- fread(fn_rna, header=T, sep="\t", data.table=F)


###
### DAR results 
infn <- paste(outdir2, "2_psycho_plotData.comb.txt.gz", sep="")
res <- fread(infn, header=T, data.table=F)
 
var2 <- sort(unique(res$psycho_variable))
varSel <- var2[c(1, 3, 6)]

###
res <- res%>%filter(psycho_variable%in%varSel)%>%
    mutate(comb=paste(Cluster, psycho_variable, sep="_"))


###
### cluster 
## df_cl <- data.frame(cluster_atac=paste("C", 0:7, sep=""))

df_cl <- read_tsv("./Cluster/TableS_summ.tsv", show_col_type=FALSE)

cluster_rna <- gsub("R", "C", df_cl$RNA_cluster)
names(cluster_rna) <- gsub("A", "C", df_cl$ATAC_cluster)
## df_cl <- df_cl%>%mutate(cluster_rna=cl_rna[cluster_atac])

comb <- sort(unique(res$comb))
res_new <- NULL
summ_all <- NULL
for (ii in comb){
   ##
   res0 <- res%>%filter(comb==ii)
   var0 <- res0$psycho_variable[1]
   cl_atac <- res0$Cluster[1]
   cl_rna <- cluster_rna[cl_atac] 

   ####
   DEG <- res_rna%>%filter(padj<0.1, var==var0, cluster==cl_rna)%>%pull(identifier)%>%unique()
   peak2 <- anno%>%filter(gene%in%DEG)%>%pull(peak)%>%unique()          
   ###if ( length(peak2)>0){
        
       res0 <- res0%>%mutate(is_gr=ifelse(gene%in%peak2, "gr1", "gr0"))   

       ### stratified group qq plots
       gr <- unique(res0$is_gr)
       res0 <- map_dfr(gr, function(ii){
           ##
           res00 <- res0%>%filter(is_gr==ii)
           th0 <- res0%>%drop_na(p.adjusted)%>%pull(baseMean)%>%min()

           ###
           iwhich <- which(res00$baseMean>th0)
           pval <- res00[iwhich, "p.value"]
           padj_BH <- p.adjust(pval, method="BH")
        
           res00$padj2 <- NA
           res00$padj2[iwhich] <- padj_BH
        
           ###
           ###
           ntest <- nrow(res00)
           res00 <- res00%>%arrange(p.value)%>%
               mutate(observed_gr=-log10(p.value), expected_gr=-log10(ppoints(ntest)))         
           res00
       })

       summ0 <- data.frame(comb=ii, nsig_0=sum(res0$p.adjusted<0.1, na.rm=T), nsig_str=sum(res0$padj2<0.1, na.rm=T),
       nDEG = length(DEG), npeak_gr1=sum(res0$is_gr=="gr1"), npeak_gr0=sum(res0$is_gr=="gr0"))         
       summ_all <- rbind(summ_all, summ0)                    

       ##
       res_new <- rbind(res_new, res0)
       cat(ii, sum(res0$is_gr=="gr1"), "\n")
   ###}
   ### END if  
}

### END loop 
 

###
### save
opfn <- paste(outdir2, "3.2_plotData_grDEG.comb.txt.gz", sep="")
fwrite(res_new, file=opfn, sep="\t", na=NA)

###
opfn <- paste(outdir2, "3.2_summary_grDEG.tsv", sep="")
write_tsv(summ_all, file=opfn)



###
### make stratified qq 

summ2 <- summ_all%>%filter(nDEG>=20, npeak_gr1>0)
summ2 <- summ2%>%mutate(psycho_variable=gsub("^C._", "", comb))

opfn <- paste(outdir2, "3.2_summary_grDEG_short.tsv", sep="")
write_tsv(summ2%>%arrange(psycho_variable, .locale="en"), file=opfn)

### names
varSel <- sort(unique(res_new$psycho_variable))
var_val <- 1:length(varSel)
names(var_val) <- varSel

plotDF <- res_new%>%filter(comb%in%summ2$comb)%>%
    mutate(var_value=as.numeric(var_val[psycho_variable]),
           comb2=fct_reorder(comb, var_value))
nlen <- length(unique(plotDF$comb2))
p2 <- ggplot(plotDF, aes(x=expected_gr, y=observed_gr, color=factor(is_gr)))+
    rasterise(geom_point(size=0.06), dpi=300)+
    geom_abline(color="grey")+
    scale_color_manual(values=c("gr0"="#f1b6da", "gr1"="#d01c8b"),
        labels=c("gr0"="Not in DEG", "gr1"="In DEG"),
        guide=guide_legend(override.aes=list(size=2)))+
  facet_wrap(~comb2, scales="free", ncol=7)+  ## ncol=7       
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    theme_bw()+
    theme(legend.title=element_blank(),
          legend.text=element_text(size=9),
          legend.key.size=grid::unit(0.4, "cm"),
          ##legend.position="inside",
          ##legend.position.inside=c(0.95, 0.12),
          axis.title=element_text(size=9),
          axis.text=element_text(size=9),
          strip.text=element_text(size=9))

### 
###
figfn <- paste(outdir2, "Figure3.2_stratified_grDEG.qq.pdf", sep="")
ggsave(figfn, p2, width=12, height=6)



#############################
### Enrichment analysis 
###############################

outdir2 <-
    
fn <- paste(outdir2, "3.2_plotData_grDEG.comb.txt.gz", sep="")
res <- fread(fn, header=T, sep="\t", data.table=FALSE)

varSel <- c("ISEL_Mean", "PSS_all_mean")
clSel <- paste("C", 0:7, sep="")

res <- res%>%filter(psycho_variable%in%varSel, Cluster%in%clSel)


summ <- res%>%group_by(comb)%>%
    summarize(nsig=sum(p.adjusted<0.1, na.rm=TRUE), .groups="drop")%>%ungroup()%>%
    filter(nsig>0)

comb2 <- unique(summ$comb)
enriched <- map_dfr(comb2, function(ii){
    ##
    res0 <- res%>%filter(comb==ii)%>%mutate(is_sig=ifelse(p.adjusted<0.1, "sig", "sig_no"))
    mat <- as.matrix(table(res0$is_sig, res0$is_gr))
    mat2 <- mat[c("sig", "sig_no"), c("gr1", "gr0")]

    tmp <- fisher.test(mat2, alternative="greater")
    tmp2 <- fisher.test(mat2)
    ##
    df0 <- data.frame(cluster_var=ii,
                      interest.in=mat2[1,1], interest.not_in=mat2[1,2],
                      not_interest.in=mat2[2,1], not_interest.not_in=mat[2,2],
                      odds=tmp$estimate, pval=tmp$p.value,
                      lower=tmp2$conf.int[1], upper=tmp2$conf.int[2])
    df0
})

opfn <- paste(outdir2, "3.2_enriched.xlsx", sep="")
write.xlsx(enriched, file=opfn)
     
###
###
fn <- paste(outdir2, "3.2_enriched.xlsx", sep="")
x <- read.xlsx(fn)

opfn2 <- paste(outdir2, "3.2_enriched.tsv", sep="")
write_tsv(x, file=opfn2)

## summ_df <- res_new%>%group_by(comb, is_motif)%>%summarize(npeak=n(),.groups="drop")%>%ungroup()

## summ2 <- summ_df%>%pivot_wider(id_cols=comb, names_from="is_motif", values_from=npeak, values_fill=0)
## names(summ2)[2:3] <- c("Not in DAMs", "In DAMs")

## ### save
## opfn <- paste(outdir2, "3_summary_gr.tsv", sep="")
## write_tsv(summ2, file=opfn)


###
### Ali results


## x <- read.xlsx("./1_DiffPeak.outs/stats_DEG_Ali.xlsx")

## identical(summ2$variable, x$var)



###############################
#### MA plots for later use 
#################################

###
### for MA plots, we module load R/4.0.3
### module unload gnu9
### module load gnu7 R/4.0.3 
### In the old version of DESeq2, data.frame can be the input data of plotMA function. So We library old version DESeq2 rather than old version of R.     


## ###
## ###
## library(tidyverse)
## library(Matrix)
## library(data.table)
## library(DESeq2, lib.loc="/wsu/el7/groups/piquelab/R/4.0.3/lib64/R/library")



## outdir2 <- "./1_DiffPeak.outs/summary_Cluster_res0.07_psycho_CTRL/"

## fn <- paste(outdir2, "2_cytokines_plotData.comb.txt.gz", sep="")
## plotDF <- fread(fn, header=T, data.table=F)

## ###
## ###
## cluster <- sort(unique(plotDF$Cluster))
## varSel <- sort(unique(plotDF$psycho_variable))


## cl <- c("C6", "C7")
## var2 <- varSel[c(5,7,8)]



## ###
## ### MA plots from ggplot 
## plotDF2 <- plotDF%>%filter(baseMean>0, Cluster%in%cl, psycho_variable%in%var2)%>%
##     mutate(sig_gr=ifelse(p.adjusted<0.1, "gr1", "gr0"))

## th0 <- quantile(abs(plotDF2$estimate), probs=0.99, na.rm = T)

## plotDF2 <- plotDF2%>%
##     mutate(estimate2 = ifelse(abs(estimate)>th0, sign(estimate)*round(th0*1.5, 1), estimate))

## p0 <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=baseMean, y=estimate2, color=sig_gr))+
##     rasterise(geom_point(size=0.06), dpi=300)+
##     ## geom_point(size=0.1)+
##     scale_color_manual(values=c("gr0"="grey30", "gr1"="red"), guide = "none")+
##     facet_grid(psycho_variable~Cluster, scales="free")+
##     xlab("baseMean") + ylab("LFC")+
##     theme_bw()+
##     theme(axis.title = element_text(size=10),
##           axis.text = element_text(size=10),
##           strip.text = element_text(size=10))

## ###
## ### save 
## figfn <- paste(outdir2, "Figure4.2_cytokines.MA.png", sep="")
## ggsave(figfn, p0, width = 340, height = 450, units="px", dpi=100)




## ###
## ### MA plots from DESeq2 
## ## cl <- "C7"
## ## var2 <- varSel[c(5,7,8)]
## ## ##
## ## figfn <- paste(outdir2, "Figure4_cytokines_", cl, ".MA.png", sep="")
## ## png(figfn, width=180, height=400, res=100)
## ## par(oma=c(0, 0, 2, 0), mar=c(2, 2, 2, 2), mgp=c(2, 1, 0))
## ## x <- matrix(1:3, 3, 1)
## ## layout(x)

## ## ### 
## ## for ( ii in var2){
## ##    ##  
## ##    df0 <- plotDF%>%dplyr::filter(Cluster==cl, psycho_variable==ii)%>%
## ##       mutate(color=p.adjusted<0.1&!is.na(p.adjusted))%>% 
## ##       dplyr::select(baseMean, estimate, color, p.value, p.adjusted)
## ##    print(DESeq2::plotMA(df0[,1:3], colLine=NA, main=ii, cex.main=1, font.main=1, cex.axis=0.8, xlab="", ylab=""))
## ##    ###      
## ##    cat(cl, ii, "\n")
## ## }
## ## print(mtext(cl, side = 3, outer=T, cex = 1))
## ## dev.off()

## ###
## ### 

## outdir2 <- "./1_DiffPeak.outs/summary_Cluster_res0.07_psycho_CTRL/"
## var0 <- "PFAS"

## fn <- paste(outdir2, "2_", var0, "_plotData.comb.txt.gz", sep="")
## plotDF <- fread(fn, header=T, data.table=F)

## ###
## ###
## cluster <- sort(unique(plotDF$Cluster))
## varSel <- sort(unique(plotDF$psycho_variable))



## ###
## ### MA plots from ggplot
## cl <- c("C6", "C7")
## var2 <- varSel[c(1, 2, 9:12, 19, 20, 23)]

## plotDF2 <- plotDF%>%filter(baseMean>0, Cluster%in%cl, psycho_variable%in%var2)%>%
##     mutate(sig_gr=ifelse(p.adjusted<0.1, "gr1", "gr0"))

## th0 <- quantile(abs(plotDF2$estimate), probs=0.99, na.rm = T)

## plotDF2 <- plotDF2%>%
##     mutate(estimate2 = ifelse(abs(estimate)>th0, sign(estimate)*round(th0*1.5, 1), estimate))

## p0 <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=baseMean, y=estimate2, color=sig_gr))+
##     rasterise(geom_point(size=0.06), dpi=300)+
##     ## geom_point(size=0.1)+
##     scale_color_manual(values=c("gr0"="grey30", "gr1"="red"), guide = "none")+
##     facet_grid(Cluster~psycho_variable, scales="free")+
##     xlab("baseMean") + ylab("LFC")+
##     theme_bw()+
##     theme(axis.title = element_text(size=10),
##           axis.text = element_text(size=10),
##           strip.text = element_text(size=8))

## ###
## ### save  
## figfn <- paste(outdir2, "Figure4_", var0, "_1.MA.png", sep="")
## ggsave(figfn, p0, width = 950, height = 280, units="px", dpi=100)



## ###
## ### 
## cl <- c("C0", "C1", "C2", "C3")
## var2 <- c("PFHxS")

## plotDF2 <- plotDF%>%filter(baseMean>0, Cluster%in%cl, psycho_variable%in%var2)%>%
##     mutate(sig_gr=ifelse(p.adjusted<0.1, "gr1", "gr0"))

## th0 <- quantile(abs(plotDF2$estimate), probs=0.99, na.rm = T)

## th_mean <- quantile(plotDF2$baseMean, probs=0.99, na.rm=T)

## plotDF2 <- plotDF2%>%
##     mutate(estimate2 = ifelse(abs(estimate)>th0, sign(estimate)*round(th0*1.5, 1), estimate),
##            baseMean2 = ifelse(baseMean > th_mean, th_mean*1.5, baseMean))

 
## p2 <- ggplot(plotDF2%>%arrange(sig_gr), aes(x=baseMean2, y=estimate2, color=sig_gr))+
##     rasterise(geom_point(size=0.06), dpi=300)+
##     ## geom_point(size=0.1)+
##     scale_color_manual(values=c("gr0"="grey30", "gr1"="red"), guide = "none")+
##     facet_wrap(~Cluster, nrow=2, scales="free")+
##     xlab("baseMean") + ylab("LFC")+
##     theme_bw()+
##     theme(axis.title = element_text(size=10),
##           axis.text = element_text(size=10),
##           strip.text = element_text(size=10))


## figfn <- paste(outdir2, "Figure4_", var0, "_2.MA.png", sep="")
## ggsave(figfn, p2, width = 420, height = 420, units="px", dpi=100)



