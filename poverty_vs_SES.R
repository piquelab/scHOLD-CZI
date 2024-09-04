Poverty analysis: 
We have bulk RNA-seq data from a cohort similar to CZI/HOLD and we tested for DEG with poverty status. I am attaching here a table of the results. Could you look for correlation in the effect size (z-score) with the SES results you got for the first batch of CZI/HOLD? It would be for each cluster separately.

R
#in analysis folder
library(data.table)
library(plyr)
library(ggplot2)
library(ggpubr)

args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/HOLD-CZI_covariates_converted_dbgapID_key_n183_AR_18newvar_updated_08-05.txt","ALL","fastdemux") #for testing
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]

outFolder=paste0(base,method,"_pseudobulk_ctrl/")
figuredir=paste0(outFolder,"figures/")

cov_file=fread(args[2])
var="SES"
pov_version="nicole1"
if(pov_version=="firstrun"){
poverty <- fread(paste0(base,"SupplementalTable2-T1Both_DESeq_IndFiltT_BasicMod_3PCs_newrefeigen_newsamp204.txt"))
} else if(pov_version=="nicole1") {
  poverty <- fread(paste0(base,"handls_cole.txt"))
  names(poverty)[2] <- "GeneSymbol"
}
poverty <- transform(poverty, z=qnorm(pvalue))
poverty <- transform(poverty,z=ifelse(log2FoldChange<0, abs(z)*-1, abs(z)))

filenames <- c(paste0(project,".deseqres_SES-RNA-CTRL.txt"),paste0(project,".deseqres_SES-RNA-LPS.txt"))
data_names <- gsub(".txt", "", filenames) #remove file ending
shortnames <- gsub(".deseqres_", "", data_names)
shortnames <- gsub(project, "", shortnames)
for(i in 1:length(filenames)) assign(shortnames[i], fread(file.path(outFolder, filenames[i]),header = TRUE, sep='\t')[,analysis:=shortnames[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(shortnames, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(shortnames)
df <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) 
colnames(df) <- c('GeneSymbol', 'padj', 'pvalue', 'log2FoldChange','var','cluster','treats','analysis')
df <- transform(df, z=qnorm(pvalue))
df <- transform(df,z=ifelse(log2FoldChange<0, abs(z)*-1, abs(z)))
df[is.na(df$padj),] <- 1
poverty[is.na(poverty$padj),] <- 1

merged <- merge(df,poverty,by=c("GeneSymbol"))
threshold <- 0.1

for (treat in unique(df$treats)){
	merged_treat=subset(merged,treats==treat)
merged_treat <- transform(merged_treat, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, "3SES_sig", ifelse(padj.y<=threshold, "2poverty_sig", "1Not_Sig")))))

#didnt update this one for the colored by sig
p <- ggplot(merged_treat, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  facet_wrap(.~cluster, scales="free")+
  theme_bw()+
  geom_point(color="#3C3C3C",alpha = 0.1)+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
  ylab("SES") +
  xlab("poverty") +
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))

      png(width = 12, height = 12, file=paste0(figuredir,project,".",treat,".",pov_version,".poverty_vs_ses_log2FC.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

p <- ggplot(merged_treat, aes(x=z.x, y=z.y)) +
  facet_wrap(.~cluster, scales="free")+
  theme_bw()+
  geom_point(aes(color=sig))+ # color="#3C3C3C",alpha = 0.8
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
  ylab("SES") +
  xlab("poverty") +
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 12, height = 12, file=paste0(figuredir,project,".",treat,".poverty_vs_ses_Z.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
}

#ones significant in the poverty study, how many have nominal p<0.05 for SES in C4 and C5, separately
poverty_sig <- subset(poverty, padj<threshold)
merged_s <- merge(df,poverty_sig,by=c("GeneSymbol"))
ns <- ldply(lapply(unique(df$treats),function(treat){
  merged_treat=subset(merged_s,treats==treat)
  nomsig=subset(merged_treat, pvalue.x<0.05)
  return(nomsig)
}),data.frame)

df_LPS <- subset(df,treats=="RNA-LPS")
df_CTRL <- subset(df,treats=="RNA-CTRL")
ns <- ldply(lapply(unique(df$cluster),function(clus){
    df_clus=subset(df,cluster==clus)
  df_w <- reshape2::dcast(df_clus, GeneSymbol~treats,value.var="z")
  colnames(df_w) <- gsub("RNA-","",colnames(df_w))
  c <-cor.test(df_w$LPS,df_w$CTRL,method="spearman")
  cdf <- data.frame(rho=round(c$estimate,2),pvalue=c$p.value,cluster=clus)
  return(cdf)
  }),data.frame)

  dfn <- na.omit(df)
  df_w <- reshape2::dcast(df, GeneSymbol~treats,value.var="z")
  colnames(df_w) <- gsub("RNA-","",colnames(df_w))
  c <-cor.test(df_w$LPS,df_w$CTRL,method="spearman")
  cdf <- data.frame(rho=round(c$estimate,2),pvalue=c$p.value)


# plots also for income and education separately
#income=DSES_09,education=DSES_07
for (var in c("DSES_09","DSES_07")){
  cat("running ",var,"\n")
filenames <- c(paste0("deseqres_",var,"-RNA-CTRL.txt"),paste0("deseqres_",var,"-RNA-LPS.txt"))
data_names <- gsub(".txt", "", filenames) #remove file ending
shortnames <- gsub("deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(shortnames[i], fread(file.path(outFolder, filenames[i]),header = FALSE, sep='\t')[,analysis:=shortnames[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(shortnames, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(shortnames)
df <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) 
colnames(df) <- c('GeneSymbol', 'padj', 'pvalue', 'log2FoldChange','var','cluster','treats','analysis')
df <- transform(df, z=qnorm(pvalue))
df <- transform(df,z=ifelse(log2FoldChange<0, abs(z)*-1, abs(z)))
df[is.na(df$padj),] <- 1
poverty[is.na(poverty$padj),] <- 1

merged <- merge(df,poverty,by=c("GeneSymbol"))
threshold <- 0.1
for (treat in unique(df$treats)){
  merged_treat=subset(merged,treats==treat)
merged_treat <- transform(merged_treat, sig=as.factor(ifelse(padj.x<threshold & padj.y<threshold, "4both", ifelse(padj.x<=threshold, paste0("3",var,"_sig"), ifelse(padj.y<=threshold, "2poverty_sig", "1Not_Sig")))))

p <- ggplot(merged_treat, aes(x=z.x, y=z.y)) +
  facet_wrap(.~cluster, scales="free")+
  theme_bw()+
  geom_point(aes(color=sig))+ # color="#3C3C3C",alpha = 0.8
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", fill = NA)+
  ylab(paste0(var)) +
  xlab("poverty") +
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 12, height = 12, file=paste0(figuredir,treat,".poverty_vs_",var,"_Z.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
} 
}
