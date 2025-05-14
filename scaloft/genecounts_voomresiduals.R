
library(tidyverse)
library(edgeR)
library(limma)
library(annotables)
library(data.table)

#load in original counts info for clusters and treats used
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/")
opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.RData")
load(opfn)
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/"

# load annotation
u_eigenvec2 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/covariates/eigenvec2_u.txt")
cv <- unique(u_eigenvec2)
##check for individuals with more than one sample
cv %>% count(Sample_ID) %>% dplyr::filter(n>1)

library("AnnotationHub")
ah <- AnnotationHub()
if(length(ah["AH98047"]) == 0) {
  edb <- ah[["AH75011"]]
} else {
  edb <- ah[["AH98047"]]
}
geneIDs <- genes(edb) %>%
  as.data.frame() %>% 
  setDT(keep.rownames = "ensembl_gene_id") %>%
  .[, c("ensembl_gene_id","entrezid","symbol","seqnames","start","end","strand","gene_biotype", "description")]
names(geneIDs)[c(1,2,4,8)] <- c("ensgene","entrez","chr","biotype")
geneIDs <- subset(geneIDs, biotype=="protein_coding")
geneIDs <- transform(geneIDs, chr=as.character(chr),strand=as.character(strand))
#this 1bp positioning for bed file is based on examples from qtltools and tensorqtl
geneIDs <- transform(geneIDs, strand_start=ifelse(strand=="+",start,start+1),strand_end=ifelse(strand=="+",end,end-1))
geneIDs <- subset(geneIDs, chr %in% c(1:22))


cluster="C0"
treat="CTRL"
clusters <- names(counts_ls)
treatments <- unique(metadata_ls[[1]]$treats)
lapply(clusters,function(cluster){
#  lapply(treatments,function(treat){
    cat("running ",cluster,treat,"\n")

countFile <- fread(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/fastQTL/phenotypes_fastqtl.",cluster,".",treat,".bed"))
data1 <- countFile[,-c(1:3)]
data <- data1[data1$ensgene %in% geneIDs$ensgene,]
#data[is.na(data)] <- 0

## Normalization of data
# make edgeR object
dge <- DGEList(counts=data)
#sum(colnames(dge$counts)==cv$Sample_ID) #does the order need to match or is this just a check for how many samples have all info
sum(colnames(dge$counts) %in% cv$Sample_ID)

#Transform counts to counts per million
dge <- calcNormFactors(dge)
cpm <- cpm(dge)
rownames(cpm) <- rownames(data)
samples <- dim(data)[2]
#Remove genes that are lowly expressed
table(rowSums(dge$counts==0)==samples) #Shows how many transcripts have 0 count across all samples
# as did GTEx:
keep.exprs <- rowSums(cpm>=0.1)>=(0.2*samples) & rowSums(data>=6)>=(0.2*samples) #Only keep transcripts that have cpm>0.1 in at least 20% of the samples
#Julong does not use the 6 count min
dim(dge) #before filter
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge) #genes that pass filter
# subset samples in cv to match samples in express data
cv_d <- subset(cv, Sample_ID %in% colnames(data))

design <- model.matrix(~0+ as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))


v_e <- voom(dge, design, plot=FALSE, normalize.method="quantile")
genes_normed_baseline <- data.frame(v_e$E)
##
X <- design
H <- X %*% solve(t(X) %*% X) %*% t(X)
dim(H)
He <- (diag(rep(1,ncol(H)))-H)
Res <- as.matrix(genes_normed_baseline) %*% He
sum(abs(t(He)-He)) #Should be almost 0
#3.159864e-12
#save the residuals:
Res <- data.frame(Res)
colnames(Res) <- colnames(dge)
rownames(Res) <- dge$genes$ensgene


#compare to qnorm
mat_qnorm <- apply(dge,1,function(x){qqnorm(rank(x, ties.method = "random"), plot = F)$x})
#model 
# Fit the linear model
model <- lm(mat_qnorm ~ design, data = data.frame(mat_qnorm, cv_d))
# Extract residuals
residuals <- t(residuals(model))
Res1 <- data.frame(residuals)
colnames(Res1) <- colnames(dge)
rownames(Res1) <- dge$genes$ensgene

library(ggmatplot)
  png(width = 10, height = 10, file=paste0(figuredir,cluster,"_voomvsqnorm.png"), pointsize=12, 
        bg = "transparent", canvas = "white", units = "in", res = 1200)
   ggmatplot(Res, Res1, pch = 19, col = "black",xlab = "voom", ylab = "qnorm") # pch and col are ggplot2 arguments
dev.off()

Res$ensg <- rownames(Res)
Res1$ensg <- rownames(Res1)
Resm <- melt(Res)
Res1m <- melt(Res1)
names(Resm)[3] <- "voom"
names(Res1m)[3] <- "qnorm_lm"
Resb <- merge(Resm,Res1m)

png(width = 10, height = 10, file=paste0(figuredir,cluster,"_voomvsqnorm.png"), pointsize=12, 
    bg = "transparent", canvas = "white", units = "in", res = 1200)
p <- ggplot(Resb, aes(x=voom, y=qnorm)) +
theme_bw()+
geom_point()+ #aes(color=sig)
geom_smooth(method = "lm", se = FALSE)+
stat_cor(color="blue",method="pearson",cor.coef.name = "R", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
  theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
  axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
  axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
  legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
print(p)
dev.off()

#matrix residuals remove from qnorm
X <- design
genes_normed_baseline <- data.frame(t(mat_qnorm))

H <- X %*% solve(t(X) %*% X) %*% t(X)
dim(H)
He <- (diag(rep(1,ncol(H)))-H)
Res2 <- as.matrix(genes_normed_baseline) %*% He
sum(abs(t(He)-He)) #Should be almost 0
#3.159864e-12

#save the residuals:
Res2 <- data.frame(Res2)
colnames(Res2) <- colnames(dge)
rownames(Res2) <- dge$genes$ensgene

Res2$ensg <- rownames(Res2)
Res2m <- melt(Res2)

names(Res2m)[3] <- "qnorm_matrix"

Resb <- merge(Res1m,Res2m)

png(width = 10, height = 10, file=paste0(figuredir,cluster,"_qnormlmvsqnormmatrix.png"), pointsize=12, 
    bg = "transparent", canvas = "white", units = "in", res = 1200)
p <- ggplot(Resb, aes(x=qnorm_lm, y=qnorm_matrix)) +
theme_bw()+
geom_point()+ #aes(color=sig)
geom_smooth(method = "lm", se = FALSE)+
stat_cor(color="blue",method="pearson",cor.coef.name = "R", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
  theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
  axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
  axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
  legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
print(p)
dev.off()

identical(Res1, Res2) #false
all.equal(Res1, Res2) #true

Resb[!Resb$qnorm_lm==Resb$qnorm_matrix,]
row.values <- Resb[which(Resb$qnorm_lm != Resb$qnorm_matrix),]
Resb[(Resb$qnorm_lm != Resb$qnorm_matrix) | is.na(Resb$qnorm_lm) | is.na(Resb$qnorm_matrix),]