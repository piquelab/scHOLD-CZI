library("limma")

## Run limma
fit <- lmFit(transformedCov, mod)
fit0 <- lmFit(transformedCov, mod0)

## Determine DE status for the regions
## Also in https://github.com/LieberInstitute/jaffelab with help and examples
getF <- function(fit, fit0, theData) {
    rss1 <- rowSums((fitted(fit) - theData)^2)
    df1 <- ncol(fit$coef)
    rss0 <- rowSums((fitted(fit0) - theData)^2)
    df0 <- ncol(fit0$coef)
    fstat <- ((rss0 - rss1) / (df1 - df0)) / (rss1 / (ncol(theData) - df1))
    f_pval <- pf(fstat, df1 - df0, ncol(theData) - df1, lower.tail = FALSE)
    fout <- cbind(fstat, df1 - 1, ncol(theData) - df1, f_pval)
    colnames(fout)[2:3] <- c("df1", "df0")
    fout <- data.frame(fout)
    return(fout)
}

ff <- getF(fit, fit0, transformedCov)

## Get the p-value and assign it to the regions
limma <- regionMat$chr21$regions
limma$fstat <- ff$fstat
limma$pvalue <- ff$f_pval
limma$padj <- p.adjust(ff$f_pval, "BH")


lmcoefs <- ldply(lapply(names(counts_ls),function(cluster){
    #cluster="C6"
    cluster_metadata_sce <- metadata_ls[[cluster]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata <- transform(cluster_metadata, treats=as.factor(treats))
    opfn <- paste0(baseoutFolder,project,".",resset,".",dimset,".ComBat_seq.",cluster,".",combatrun,".RData")
    load(opfn)
    vars <- ldply(lapply(secondrunvars,function(var){
        cat("running ", cluster, var, "\n")
    df <- ldply(lapply(1:length(contrastdf$control),function(x){
        con=contrastdf[x,]
        contrast=paste0(con$treatment,"_vs_",con$control)
        cluster_metadata <- within(cluster_metadata, treats <- relevel(treats, ref = con$control))
        cluster_metadata_var <- subset(cluster_metadata, treats %in% c(con$control,con$treatment)) #removing due to low ind counts
        cluster_metadata_var <- cluster_metadata_var[,c("Sample_ID","Wave","genPC1","genPC2","genPC3","Sex","cage1","treats",var)]
        cluster_metadata_var <- cluster_metadata_var[complete.cases(cluster_metadata_var), ] #if there are missing covariates, this removes those individuals as deseq can't handle NAs
        #cdf <- subset(plyr::count(cluster_metadata_var$Sample_ID),freq<2)
        #cluster_metadata_var <- subset(cluster_metadata_var, !Sample_ID %in% cdf$x)
        cluster_counts_t <- adjusted[,which(colnames(adjusted) %in% rownames(cluster_metadata_var))]
        cv_d <- cluster_metadata_var[which(rownames(cluster_metadata_var) %in% colnames(cluster_counts_t)),]
        cv_d <- transform(cv_d, treats=factor(treats))
        ## Normalization of data
        # make edgeR object
        dge <- DGEList(counts=cluster_counts_t)
        #Transform counts to counts per million
        dge <- calcNormFactors(dge)
        cpm <- log2(cpm(dge)+1)
        rownames(cpm) <- rownames(cluster_counts_t)
        samples <- dim(cluster_counts_t)[2]
        varcol=cv_d[,var]
        design <- model.matrix(~ as.numeric(varcol):cv_d$treats + as.numeric(varcol) + cv_d$treats + as.factor(cv_d$Sex) + as.numeric(cv_d$cage1) + factor(cv_d$Wave) + as.numeric(cv_d$genPC1) + as.numeric(cv_d$genPC2) +as.numeric(cv_d$genPC3))
        transformedCov <- t(cpm)
        fit <- lmFit(cpm, design)
fit2 <- eBayes(fit)
rslt <- decideTests(fit2)

        fit <- eBayes(fit, trend = TRUE, robust = TRUE)
        tt <- topTable(fit, n = Inf, adjust.method = "BH")

fit$coefficients ##extract coefficients
fit$p.value ## p-values
data.frame(cluster=cluster, contrast=contrast, variable=var, fit$coefficients, fit$p.value)
coefstreats <- grepl("treats", colnames(fit$coefficients))
coefsvarcol <- grepl("varcol", colnames(fit$coefficients))
pvaltreats <- grepl("treats", colnames(fit$p.value))
pvalvarcol <- grepl("varcol", colnames(fit$p.value))

int <- data.frame(coef=tt[,coefstreats&coefsvarcol],pval=fit$p.value[,pvaltreats&pvalvarcol])
        limmamodel <- data.frame(cluster=cluster, contrast=contrast, variable=var, term=ifelse(grepl("treats", term) & grepl("varcol", term), paste0(contrast,":",var), ifelse(grepl("treats", term),contrast, ifelse(grepl("varcol", term), var, term))))

        model <- lm(t(cpm) ~ design, data = data.frame(t(cpm), cv_d))
        tmodel <- as.data.frame(tidy(model))
        tmodel <- transform(tmodel, cluster=cluster, contrast=contrast, variable=var, term=ifelse(grepl("treats", term) & grepl(var, term), paste0(contrast,":",var), ifelse(grepl("treats", term),contrast, ifelse(grepl(var, term), var, term))))
        tmodeltreat <- subset(tmodel, term==paste0(contrast,":",var) | term==contrast | term==var)
        tmodeltreata <- ldply(lapply(unique(tmodeltreat$term),function(t){
            dft <- subset(tmodeltreat, term==t)
            names(dft)[1] <- "identifier"
            dft <- transform(dft, padj=p.adjust(p.value,method="BH"))
            fwrite(dft, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",cluster,".",var,".lm_",t,".treatvarint.txt"))
            table <- data.frame(variable= var, cluster=cluster,term=t)
        table$number_samples <- paste(nrow(cluster_metadata_var))
        table$number_individuals <- paste(length(unique(cluster_metadata_var$Sample_ID)))
    table$gene_number <- paste(nrow(cluster_counts_t))
    table$tested_genes <- paste(nrow(dft))
    table$DEGs_FDR <- paste(length(which(dft$padj<0.05)))
    table$DEGs_FDR_10 <- paste(length(which(dft$padj<0.1)))
    fwrite(table, sep='\t', quote=F, row.names=F, col.names=T, paste0(outFolder,project,".",resset,".",dimset,".",cluster,".",var,".lmstats_",t,".treatvarint.txt"))
        return(dft)
        }),data.frame)
        return(tmodeltreata)
    }), data.frame)
    return(df)
    }), data.frame)
    return(vars)
}), data.frame)