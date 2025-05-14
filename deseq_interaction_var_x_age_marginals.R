# Pull out variable number as contrast and age number as contrast (do we recover signal when correcting for agexvavr?)
marginal_outFolder <- "/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/adjusted/0.2.11/SES_PCs_sex_age_and_treats_generem_ageint_marginals/"
if (!file.exists(marginal_outFolder)) dir.create(marginal_outFolder, showWarnings=F)
if (!file.exists(paste0(marginal_outFolder,"stats/"))) dir.create(paste0(marginal_outFolder,"stats/"), showWarnings=F)
if (!file.exists(paste0(marginal_outFolder,"deseqres/"))) dir.create(paste0(marginal_outFolder,"deseqres/"), showWarnings=F)
if (!file.exists(paste0(marginal_outFolder,"figures/"))) dir.create(paste0(marginal_outFolder,"figures/"), showWarnings=F)

for (cluster in names(counts_ls)){
    #cluster="C0"
        cat("running ", cluster, "\n")
    fdr=0.05
    for(i in c("RNA-CTRL","RNA-LPS")){
        mclapply(secondrunvars,function(var){
            contrast=paste0(var,":age")
            opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_output-",contrast,".",cluster,".",i,".",run,".RDS")
            dds <- read_rds(opfn)
            contrast_marginal=paste0(var,"_marginal")
            res <- results(dds, name =var)
            sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
            names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
            sub.table <- sub.table[!is.na(sub.table$padj), ]
            cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
            sub.table$var=var
            sub.table$cluster=cluster
            sub.table$treat=i
            sub.table$contrast_marginal=contrast_marginal
            fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(marginal_outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",contrast_marginal,".",i,".",run,".txt"))
            sigDEGs <- subset(sub.table,padj<fdr)
            table <- data.frame(symb=var, description=variables_df[variables_df$variable==var,]$description, cluster=cluster,treat=i,contrast=contrast_marginal)
            table$number_samples <- paste(dim(dds@colData)[1])
            table$number_individuals <- paste(length(unique(dds@colData$Sample_ID)))
            table$gene_number <- paste(length(dds@rowRanges))
            table$tested_genes <- paste(nrow(sub.table))
            table$DEGs_FDR <- paste(nrow(sigDEGs))
            table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(marginal_outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast_marginal,".",i,".",run,".txt"))

            contrast_marginal=paste0(var,".age_marginal")
            res <- results(dds, name ="age")
            sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', res$'lfcSE',stringsAsFactors=FALSE)
            names(sub.table) <- c('identifier', 'padj', 'pvalue', 'logFC','SE')
            sub.table <- sub.table[!is.na(sub.table$padj), ]
            cat("BH diff. expressed ids.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
            sub.table$var=var
            sub.table$cluster=cluster
            sub.table$treat=i
            sub.table$contrast_marginal=contrast_marginal
            fwrite(sub.table, sep='\t', quote=F, row.names=F, col.names=T, paste0(marginal_outFolder,"deseqres/",project,".",resset,".",dimset,".",cluster,".deseqres_",contrast_marginal,".",i,".",run,".txt"))
            sigDEGs <- subset(sub.table,padj<fdr)
            table <- data.frame(symb=var, description=variables_df[variables_df$variable==var,]$description, cluster=cluster,treat=i,contrast=contrast_marginal)
            table$number_samples <- paste(dim(dds@colData)[1])
            table$number_individuals <- paste(length(unique(dds@colData$Sample_ID)))
            table$gene_number <- paste(length(dds@rowRanges))
            table$tested_genes <- paste(nrow(sub.table))
            table$DEGs_FDR <- paste(nrow(sigDEGs))
            table$DEGs_FDR_10 <- paste(nrow(subset(sub.table,padj<0.1)))
            fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(marginal_outFolder,"stats/",project,".",resset,".",dimset,".",cluster,".stats_all_cell_types-",contrast_marginal,".",i,".",run,".txt"))
        })
    }
}

for (var in secondrunvars){
    for (treat in c("RNA-CTRL","RNA-LPS")){
cat("running",var,treat,"\n")
myDir <- paste0(marginal_outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory 
filenames <- filenames[grep(paste0("-",var,"_marginal.",treat),filenames)]
if (length(filenames)>0){
data_names <- gsub(paste0(".",treat,".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("stats_all_cell_types-", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(stats, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".stats_all_cell_types-",var,"_marginal.",run,".",treat,".ageinteraction.txt"))

myDir <- paste0(marginal_outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0("_",var,"_marginal.",treat),filenames)]
data_names <- gsub(paste0(".",treat,".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".deseqres_",var,"_marginal.",run,".",treat,".ageinteraction.txt"))
}
}
}

for (var in secondrunvars){
    for (treat in c("RNA-CTRL","RNA-LPS")){
cat("running",var,treat,"\n")
myDir <- paste0(marginal_outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory 
filenames <- filenames[grep(paste0("-",var,".age_marginal.",treat),filenames)]
if (length(filenames)>0){
data_names <- gsub(paste0(".",treat,".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("stats_all_cell_types-", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(stats, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".stats_all_cell_types-",var,".age_marginal.",run,".",treat,".ageinteraction.txt"))

myDir <- paste0(marginal_outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(paste0("_",var,".age_marginal.",treat),filenames)]
data_names <- gsub(paste0(".",treat,".",run,".txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,"."), "", data_names)
data_names <- gsub("deseqres_", "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
fwrite(deseqres, sep='\t', quote=F, row.names=F, col.names=T, paste0(myDir,project,".",resset,".",dimset,".deseqres_",var,".age_marginal.",run,".",treat,".ageinteraction.txt"))
}
}
}

for (treat in c("RNA-CTRL","RNA-LPS")){
myDir <- paste0(marginal_outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".ageinteraction.txt",filenames)]
filenames <- filenames[grep(paste0(".age_marginal"),filenames,invert=T)]
filenames <- filenames[grep(paste0(treat),filenames)]
data_names <- gsub(paste0(".",run,".",treat,".ageinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".deseqres_"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#colnames(deseqres) <- c("identifier","padj","pvalue","logFC","SE","var","cluster","contrast")

myDir <- paste0(marginal_outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".ageinteraction.txt",filenames)]
filenames <- filenames[grep(paste0(".age_marginal"),filenames,invert=T)]
filenames <- filenames[grep(paste0(treat),filenames)]
data_names <- gsub(paste0(".",run,".",treat,".ageinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".stats_all_cell_types-"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does

var50 <- subset(stats, DEGs_FDR_10>50)
stats50 <- subset(stats, contrast %in% var50$contrast)
wstats50 <- reshape(stats50[,-c(6,8,10)], idvar = c("symb","description","contrast"), timevar = "cluster", v.names=c("number_individuals","tested_genes","DEGs_FDR_10"), direction = "wide",sep=":")

 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

degcols <- grep("DEGs_FDR_10",colnames(wstats50))
stats50degcols <- degcols[ apply(wstats50[,degcols],MARGIN=2,FUN=my.max)<50]
#subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit( colnames(wstats50)[stats50degcols],"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(wstats50))] #old version
subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit( colnames(wstats50)[stats50degcols],"[.]"),function(y) y[1]),"NA"),collapse="|"),colnames(wstats50))]
dft <- subsubvars %>% flextable() %>% span_header(sep=":")
bodycol=ncol(subsubvars)-3
cols <- seq(4,bodycol,by=3)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1:2, j = c(4:length(subsubvars)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(cols), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(marginal_outFolder,"figures/",project,".",resset,".",dimset,".ageinteraction_variablemarginal.",run,".",treat,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("number_individuals|tested_genes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() %>% split_header(sep=":")
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = c(4), border = border, part = "all")
dft <- align(dft, i = 2, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars_degonly)), align = "center", part = "body")

flextable::save_as_image(
  dft,
  path = paste0(marginal_outFolder,"figures/",project,".",resset,".",dimset,".ageinteraction_variablemarginal.",run,".",treat,".ndegonly.png"))
}

for (treat in c("RNA-CTRL","RNA-LPS")){
myDir <- paste0(marginal_outFolder,"deseqres/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".ageinteraction.txt",filenames)]
filenames <- filenames[grep(paste0(".age_marginal"),filenames)]
filenames <- filenames[grep(paste0(treat),filenames)]
data_names <- gsub(paste0(".",run,".",treat,".ageinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".deseqres_"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
deseqres <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does
#colnames(deseqres) <- c("identifier","padj","pvalue","logFC","SE","var","cluster","contrast")

myDir <- paste0(marginal_outFolder,"stats/")
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep(".ageinteraction.txt",filenames)]
filenames <- filenames[grep(paste0(".age_marginal"),filenames)]
filenames <- filenames[grep(paste0(treat),filenames)]
data_names <- gsub(paste0(".",run,".",treat,".ageinteraction.txt"), "", filenames) #remove file ending
data_names <- gsub(paste0(project,".",resset,".",dimset,".stats_all_cell_types-"), "", data_names)
for(i in 1:length(filenames)) assign(data_names[i], fread(file.path(myDir, filenames[i]),header = TRUE, sep='\t')) #read in specific files and set the df object names. can dro0p unwanted columns
ddf <- lapply(data_names, function(x) get(x)) #grab data from list of df names
names(ddf) <- c(data_names)
stats <- as.data.frame(Reduce(function(x, y) rbind(x, y),ddf,accumulate=F)) #recursively merge df from list of df. still not sure what accumulate does

var50 <- subset(stats, DEGs_FDR_10>50)
stats50 <- subset(stats, contrast %in% var50$contrast)
wstats50 <- reshape(stats50[,-c(6,8,10)], idvar = c("symb","description","contrast"), timevar = "cluster", v.names=c("number_individuals","tested_genes","DEGs_FDR_10"), direction = "wide",sep=":")

 my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)

degcols <- grep("DEGs_FDR_10",colnames(wstats50))
stats50degcols <- degcols[ apply(wstats50[,degcols],MARGIN=2,FUN=my.max)<50]
#subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit( colnames(wstats50)[stats50degcols],"[.]"),function(y) y[1]),"NA"),collapse="[.]|"),colnames(wstats50))] #old version
subsubvars <- wstats50[,!grepl(paste(c(sapply(strsplit( colnames(wstats50)[stats50degcols],"[.]"),function(y) y[1]),"NA"),collapse="|"),colnames(wstats50))]
dft <- subsubvars %>% flextable() %>% span_header(sep=":")
bodycol=ncol(subsubvars)-3
cols <- seq(4,bodycol,by=3)
border <- fp_border()
big_border <- fp_border(color = "black", width = 2)
dft <- align(dft, i = 1:2, j = c(4:length(subsubvars)), align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars)), align = "center", part = "body")
dft <- vline(dft, i = NULL, j = c(cols), border = border, part = "all")
dft <- border_outer(dft, part = "all", border = big_border)

flextable::save_as_image(
  dft,
  path = paste0(marginal_outFolder,"figures/",project,".",resset,".",dimset,".ageinteraction_agemarginal.",run,".",treat,".50degs.png"))

subsubvars_degonly <- subsubvars[,!grepl("number_individuals|tested_genes",colnames(subsubvars))]
dft <- subsubvars_degonly %>% flextable() %>% split_header(sep=":")
dft <- border_outer(dft, part = "all", border = big_border)
dft <- vline(dft, i = NULL, j = c(4), border = border, part = "all")
dft <- align(dft, i = 2, j = NULL, align = "center", part = "header")
dft <- align(dft, i = NULL, j = c(4:length(subsubvars_degonly)), align = "center", part = "body")

flextable::save_as_image(
  dft,
  path = paste0(marginal_outFolder,"figures/",project,".",resset,".",dimset,".ageinteraction_agemarginal.",run,".",treat,".ndegonly.png"))
}


#compare marginals to baseline
distinct_colors <-c("#db6e00",
"#4329cb",
"#5fab00",
"#c331df",
"#019b32",
"#b95aff",
"#c2a700",
"#0036b7",
"#f9bc28",
"#6d0096",
"#79db87",
"#ff7dfa",
"#bcd05c",
"#014ab7",
"#9a9200",
"#688bff",
"#ff3123",
"#01b281",
"#a80044",
"#01aea8",
"#c14600",
"#57aeff",
"#6d0004",
"#0086ce",
"#ff9b6c",
"#0155a6",
"#6b6600",
"#1f2865",
"#d7b282",
"#81004c",
"#006a34",
"#d5bbf9",
"#273a00",
"#ffabc9",
"#724500",
"#714771",
"#ff928d",
"#55133b")
baselinerun="SES_PCs_sex_age_and_treats_generem"
marginalrun=run
baselineoutFolder=paste0(baseoutFolder,baselinerun,"/")
for (i in c("RNA-CTRL","RNA-LPS")){
    varsrun<-na.omit(ldply(lapply(variables_df$variable, function(var){
    if(isTRUE(file.size(paste0(marginal_outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"_marginal.",marginalrun,".",i,".ageinteraction.txt")) > 0)){
    df <- data.frame(variable=var)
    } else {
    df <- data.frame(variable=NA)
    }
    return(df)
        }),data.frame))
    subvars <- ldply(lapply(varsrun$variable, function(var){
    cat("running",i,var,"\n")
    if(file.exists(paste0(baselineoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",baselinerun,".txt"))){
        stats <- fread(paste0(baselineoutFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"-",i,".",baselinerun,".txt"))
        corrected_stats <- fread(paste0(marginal_outFolder,"stats/",project,".",resset,".",dimset,".stats_all_cell_types-",var,"_marginal.",marginalrun,".",i,".ageinteraction.txt"))
        merged_treat <- merge(corrected_stats,stats,by=c("cluster","symb"))
        names(merged_treat)[c(11,18)] <- c("age_var_corrected_DEGs","DEGs")
        return(merged_treat)
    }
}),data.frame)

p <- ggplot(subvars, aes(x=DEGs, y=age_var_corrected_DEGs)) +
  theme_bw()+
  geom_point(size=7,aes(shape=cluster,color=variable,fill=variable))+ #aes(color=sig)
  scale_shape_manual(values = c(7,8,15:18,25))+
  scale_color_manual(values = distinct_colors) +
  #scale_fill_discrete(guide="none")+
  geom_vline(xintercept = 50,linetype="dotted",colour="red")+
  geom_hline(yintercept = 50,linetype="dotted",colour="red")+
  geom_abline()+
#  geom_text(aes(x=-Inf, y=Inf, hjust=-0.2, vjust=1.2, label = r2_eqn(lm(as.numeric(AF) ~ as.numeric(value)))), data=bbinom_ind_pop[bbinom_ind_pop$cell %in% c(unique(snp1$cell))], parse = TRUE, size=6,colour="black") 
  stat_cor(color="blue",method="spearman",cor.coef.name = "rho", size=6, label.sep="\n", r.digits=2,na.rm=T)+ #label.x = -6,label.y = 5
    theme(panel.background = element_rect(fill="white",colour = "black",size=1.3),
    axis.text.x = element_text(colour = "black",size = rel(1.3)),axis.text.y = element_text(colour = "black",size = rel(1.3)),
    axis.title.y = element_text(colour = "black",size = rel(1.5)),axis.title.x = element_text(colour = "black",size = rel(1.5)),
    legend.text=element_text(size = rel(1.3)),legend.title=element_text(size = rel(1.5)),strip.text.x = element_text(size = rel(1.3))) #+ coord_cartesian(ylim = c(-8,8), xlim = c(-8,8))
png(width = 10, height = 10, file=paste0(marginal_outFolder,"figures/",project,".",resset,".",dimset,".",i,".agexvarcorrectedvsbaseline_numDEGs.png"), pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 1200)
print(p)
dev.off()
}

