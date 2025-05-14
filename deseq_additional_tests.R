remotes::install_github("js229/Vennerable")
library(Vennerable)
library(data.table)
library(ggplot2)
library(ggpubr)
library(plyr);library(dplyr)

args <- c("CZI","/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_ISEL_Mean_03_26_2024.txt","ALL","fastdemux",0.2,11) #for testing
job <- args[1]
base <- args[2]
cov_file=fread(args[3]) #this is the psych cov file
project=args[4]
method=args[5]
resset <- args[6]
dimset=args[7]
sample_batch <- args[8]
zingeRrun=FALSE
outFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/")
figuredir=paste0(outFolder,"figures/")

######
noSESrun="SES_PCs_sex_age_and_treats_generem"
SESrun="SES_PCs_SES_sex_age_and_treats_generem"

combatrun="SES_PCs_sex_age_and_treats_adjusted_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
noSESoutFolder=paste0(baseoutFolder,noSESrun,"/")
SESoutFolder=paste0(baseoutFolder,SESrun,"/")
variablesL <- c("DSES_09","chronic_sum","SES","cytocomp","PSS_all_mean", "Chol_HDL", "BPd_avg", "isel", "pr_comp", 
    "SNI_NumPeople_r", "LogDED", "LogTrig", "nii_mean") 

for (var in c(psychvarstorun,"factor_HS_CRP")){
    #for (var in c(variablesL)){
    for (i in c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX")){
        cat("running ",var," ",i,"\n")

noSESres <- fread(paste0(noSESoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",noSESrun,".txt"))
SESres <- fread(paste0(SESoutFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",i,".",SESrun,".txt"))

noSESres <- transform(noSESres, model="noSES")
SESres <- transform(SESres, model="SES")

noSESressig <- subset(noSESres,padj<0.1)
SESressig <- subset(SESres,padj<0.1)

if(dim(noSESressig)[1]>0 & dim(SESressig)[1]>0){
df <- unique(rbind(noSESressig[,c("identifier","model")],SESressig[,c("identifier","model")]))

d_list <- split(df$identifier,df$model)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
png(width =6, height = 6, file=paste0(noSESoutFolder,"figures/",project,".",var,"-",i,".noSESvsSES_venn.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
print(p)
dev.off()
}
}}


#plotting LPS vs CONTROL
run="SES_PCs_sex_age_and_treats_generem"
#run="SES_PCs_sex_age_and_treats_generem_zingeR"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
threshold=0.10
variablesL <- c("chronic_sum","cytocomp", "SES","PSS_all_mean", "Chol_HDL", "BPd_avg", "isel", "pr_comp") 
variablesL <- variablesL[-c(3)]


treat="RNA-LPS"
control="RNA-CTRL"
subvars <- ldply(lapply(variablesL, function(var){
    cat("running",treat,"vs",control,var,"\n")
        ctrl_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",control,".",run,".txt"))
        treat_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",var,"-",treat,".",run,".txt"))
        merged_treat <- rbind(ctrl_deseqres,treat_deseqres)
        return(merged_treat)
}),data.frame)
intrun="SES_PCs_sex_age_and_treats_generem_treatint"
#intrun="SES_PCs_sex_age_and_treats_generem_zingeR_treatint"
intoutFolder=paste0(baseoutFolder,intrun,"/")
myDir <- paste0(intoutFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".treatinteraction.txt",filenames)]
data_names <- gsub(paste0(".",run,".treatinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".deseqres_"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
names(deseqres)[8] <- c("treats")
deseqres <- transform(deseqres, treats="interaction")
deseqresL <- subset(deseqres, var %in% variablesL)

all <- rbind(subvars,deseqresL)
all_sig <- subset(all, padj<0.1)
lapply(split(all_sig,all_sig$var),function(i){
    myvar=unique(i$var)
    cat("running ",myvar)
    df <- unique(i[,c("identifier","treats")])
    d_list <- split(df$identifier,df$treats)
    V_d_list <- Venn(d_list)
    Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
    #Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
    Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
    gp <- VennThemes(Vennlist)
    png(width =6, height = 6, file=paste0(outFolder,"figures/",resset,".",dimset,".",myvar,".ctrlvsplsvsint_venn.png"), pointsize=12, 
          bg = "transparent", units = "in", res = 1200)
    par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
    p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE)); grid.text(paste0(myvar), y=0.92, gp=gpar(col="black", cex=2))
    print(p)
    dev.off()
})

#plotting LPS vs CONTROL
run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
threshold=0.10
variablesL <- c("chronic_sum","cytocomp", "SES","PSS_all_mean", "Chol_HDL", "BPd_avg", "isel", "pr_comp") 
variablesL <- variablesL[-c(3)]

treat="RNA-LPS"
control="RNA-CTRL"
subvars <- ldply(lapply(variablesL, function(v){
    cat("running",treat,"vs",control,v,"\n")
    ctrl_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",control,".",run,".txt"))
    treat_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v,"-",treat,".",run,".txt"))
    merged_treat <- merge(ctrl_deseqres,treat_deseqres,by=c("identifier","cluster","var"),all=T)
    deseqressub <- subset(deseqres, var %in% v)
    merged_treat_all <- merge(merged_treat,deseqressub,by=c("identifier","cluster","var"),all.x=T)
    merged_treat_all <- transform(merged_treat_all, CTRL_zscore=logFC.x/SE.x, LPS_zscore=logFC.y/SE.y, interaction_zscore=logFC/SE, bothsig=ifelse(padj.x<threshold & padj.y<threshold,"both","notboth"))
    merged_treat_all <- merged_treat_all %>% mutate(sig = case_when(
    padj < threshold & bothsig=="both" ~ "6ALLsig",
    padj.x<threshold & padj.y<threshold ~ "5CTRLandLPS_sig",
    padj < threshold ~ "4interaction_sig",
    padj.x < threshold ~ "3CTRL_sig",
    padj.y < threshold ~ "2LPS_sig",    
    ))
    merged_treat_all$sig[is.na(merged_treat_all$sig)] <- "1Not_Sig"
    return(merged_treat_all)
}),data.frame)

cols <- c("4interaction_sig" = "orange", "6ALLsig" = "purple", "5CTRLandLPS_sig" = "green","1Not_Sig" = "grey","3CTRL_sig" = "red", "2LPS_sig" = "blue")
p <- ggplot(subvars, aes(x=CTRL_zscore, y=LPS_zscore)) +
  facet_wrap(.~var)+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
#ploting by var and cluster
p <- ggplot(subvars, aes(x=CTRL_zscore, y=LPS_zscore)) +
  facet_grid(cluster~var)+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline()+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 14, height = 14, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl_mat.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

#for logfc
subvars_sig1 <- subset(subvars, padj.x<threshold | padj.y<threshold | padj <threshold)
p <- ggplot(subvars_sig1, aes(x=logFC.x, y=logFC.y)) +
  facet_wrap(.~var, scales="free")+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  xlab("CTRL logFC") + ylab("LPS logFC")+
  geom_abline()+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
#ploting by var and cluster
p <- ggplot(subvars, aes(x=logFC.x, y=logFC.y)) +
  facet_grid(cluster~var, scales="free")+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  xlab("CTRL logFC") + ylab("LPS logFC")+
  geom_abline()+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 14, height = 14, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".lpsvsctrl_mat_logfc.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()

#comparing stress vs social support
#plotting LPS vs CONTROL
run="SES_PCs_sex_age_and_treats_generem"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/adjusted/",resset,".",dimset,"/")
outFolder=paste0(baseoutFolder,run,"/")
threshold=0.10
treat="RNA-CTRL"
v1="PSS_all_mean"
v2="ISEL_Mean"

stress_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v1,"-",treat,".",run,".txt"))
isel_deseqres <- fread(paste0(outFolder,"deseqres/",project,".",resset,".",dimset,".deseqres_",v2,"-",treat,".",run,".txt"))
merged_treat_all <- merge(stress_deseqres,isel_deseqres,by=c("identifier","cluster"),all=T)
merged_treat_all <- transform(merged_treat_all, stress_zscore=logFC.x/SE.x, socialsupport_zscore=logFC.y/SE.y)
merged_treat_all <- merged_treat_all %>% mutate(sig = case_when(
padj.x<threshold & padj.y<threshold ~ "4stress_and_socialsupport_sig",
padj.x < threshold ~ "3stress_sig",
padj.y < threshold ~ "2socialsupport_sig",    
))
merged_treat_all$sig[is.na(merged_treat_all$sig)] <- "1Not_Sig"

cols <- c("4stress_and_socialsupport_sig" = "green","1Not_Sig" = "grey","3stress_sig" = "red", "2socialsupport_sig" = "blue")
p <- ggplot(merged_treat_all, aes(x=stress_zscore, y=socialsupport_zscore)) +
  facet_wrap(.~cluster)+
  theme_bw()+
  geom_point(aes(color=sig))+ #aes(color=sig)
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  geom_abline(slope=-1)+
  scale_colour_manual(values = cols)+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
      png(width = 12, height = 12, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".",treat,".stressvssocialsupport.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()


stress_deseqres <- transform(stress_deseqres, model="Stress")
isel_deseqres <- transform(isel_deseqres, model="SocialSupport")

lapply(unique(merged_treat_all$cluster),function(c){
  stress_deseqressig <- subset(stress_deseqres,padj<0.1 & cluster==c)
  isel_deseqressig <- subset(isel_deseqres,padj<0.1 & cluster==c)

  if(dim(stress_deseqressig)[1]>0 & dim(isel_deseqressig)[1]>0){
  df <- unique(rbind(stress_deseqressig[,c("identifier","model")],isel_deseqressig[,c("identifier","model")]))

  d_list <- split(df$identifier,df$model)
  V_d_list <- Venn(d_list)
  Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
  #Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
  Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
  gp <- VennThemes(Vennlist)
  png(width =6, height = 6, file=paste0(outFolder,"figures/",project,".",resset,".",dimset,".",run,".",treat,".",c,".stressvssocialsupport_venn.png"), pointsize=12, 
        bg = "transparent", units = "in", res = 1200)
  par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
  p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
  print(p)
  dev.off()
  }
  })


###########################################################3
#################################################################
require(gridExtra)
png(width =6, height = 6, file=paste0(noSESoutFolder,"figures/",project,".",var,"-",i,".noSESvsSES_venn.png"), pointsize=12, 
      bg = "transparent", units = "in", res = 1200)
par(mar=c(5,5,4,2)+0.1) #,cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
grid.arrange(grid::grid.grabExpr(p), top=paste0(var,"-",i), bottom="subtitle")
dev.off()





vennplots <- lapply(c("RNA-CTRL","RNA-LPS","RNA-LPS-DEX"),function(i){
              cat("running ",var," ",i,"\n")

noSESres <- fread(paste0(noSESoutFolder,"deseqres/",project,".deseqres_",var,"-",i,".",noSESrun,".txt"))
SESres <- fread(paste0(SESoutFolder,"deseqres/",project,".deseqres_",var,"-",i,".",SESrun,".txt"))

noSESres <- transform(noSESres, model="noSES")
SESres <- transform(SESres, model="SES")

noSESressig <- subset(noSESres,padj<0.1)
SESressig <- subset(SESres,padj<0.1)

if(dim(noSESressig)[1]>0 | dim(SESressig)[1]>0){
df <- unique(rbind(noSESressig[,c("identifier","model")],SESressig[,c("identifier","model")]))

d_list <- split(df$identifier,df$model)
V_d_list <- Venn(d_list)
Vennlist <- compute.Venn(V_d_list, doWeights = TRUE)
#Vennlist@FaceLabels <- transform(Vennlist@FaceLabels, y=ifelse(Signature==100, 70, y))
Vennlist@FaceLabels <- distinct(Vennlist@FaceLabels,Signature, .keep_all= TRUE)
gp <- VennThemes(Vennlist)
p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))
return(Vennlist)
}
      })

plot_data_column = function (data, treat) {
      gp <- VennThemes(Vennlist)
      p <- plot(Vennlist, gpList=gp, show = list(FaceText = c("weight"), DarkMatter = TRUE))

    ggplot(data, aes_string(x = treat)) +
        geom_histogram(fill = "lightgreen") +
        xlab(column)
}

myplots <- lapply(Vennlist, plot_data_column, data = data2)