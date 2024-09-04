library(tidyverse)
library(ggpubr)

args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_HOLD01-HOLD14_dbgap.ID_cziexp_fixed_05_03_2024.txt","alternative","CZ1_group.txt") #for testing
cov_file=args[2]

#read in samples file (just list of samples to run, each sample on newline)
if(!is.na(args[4])){
samples=read.table(paste0(args[1],args[4]),header=F)
samples$Batch <- sapply(strsplit(samples$V1,"-"),function(y) y[1])
cat("samplefile=",args[4],"\n")
project=sapply(strsplit(args[4],"_"),function(y)y[1])
} else{
   project="ALL"
}

cat("project=",project,"\n")

   outFolder=paste0(args[1],"1_demux_output/")
   opfn <- paste0(outFolder,"1_demux_New.SNG.rds")
   demux <- read_rds(opfn)
   demux <- demux %>% dplyr::filter(NUM.READS>100,NUM.SNPS>100)
   demux<- transform(demux,method="demux")
   cell.counts_demux <- demux %>% group_by(EXP,Sample_ID, BATCH, treats) %>% summarize(cell_count=n()) %>%as.data.frame()

   outFolder=paste0(args[1],"1_demux_alt_output/")
   opfn <- paste0(outFolder,project,".1_demux_alt_New.SNG.rds")
   fastdemux <- read_rds(opfn)
   fastdemux<- transform(fastdemux,method="fastdemux")
   cell.counts_fastdemux <- fastdemux %>% group_by(EXP,Sample_ID, BATCH, treats) %>% summarize(cell_count=n()) %>%as.data.frame()

figuredir=paste0(outFolder,"umatched_removed/")

df <- merge(demux,fastdemux,by=c("NEW_BARCODE","EXP","BATCH","treats","Sample_ID"))
df_cellcounts <- merge(cell.counts_demux,cell.counts_fastdemux,by=c("EXP","BATCH","treats","Sample_ID"))

fig0 <- ggplot(df,aes(x=NUM.READS.x,y=NUM.READS.y))+
geom_point()+
xlab("demux")+
ylab("fastdemux")+
stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
theme_bw()

figfn <- paste(figuredir, project,".demuxvsfastdemux_numreads.png", sep="")
png(figfn, width=4000, height=4000, res=380)
fig0
dev.off()

#number of cells
fig0 <- ggplot(df_cellcounts,aes(x=cell_count.x,y=cell_count.y))+
geom_point()+
xlab("demux")+
ylab("fastdemux")+
stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
theme_bw()

figfn <- paste(figuredir, project,".demuxvsfastdemux_numcells.png", sep="")
png(figfn, width=4000, height=4000, res=380)
fig0
dev.off()

remotes::install_github("js229/Vennerable")
library(Vennerable)

df <- unique(rbind(demux[,c("NEW_BARCODE","method")],fastdemux[,c("NEW_BARCODE","method")]))

d_list <- split(df[,"NEW_BARCODE"],df$method)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0(figuredir,project,".demuxvsfastdemux_BC_venn.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
print(p)
dev.off()

