
library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)
library(plyr)

timestamp()


###################################################################
########## make histogram score distributions #####################
########### for 165 HOLD samples #################################

rm(list=ls())


setwd("/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/distributions/")
outFolder="/rs/rs_grp_schold/CZI/RNA/analysis/fastdemux_pseudobulk_ctrl/nodex/noCombat_DESeq/figures/distributions/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)


####### Load the correacter master covariate table: #############
cov.file <- "/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt"
#cov.file <- "/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_SNI_NoP_LogDED-Trig_chronic-sum-categ_nii-mean_11_25_2024.txt"
covariates <- read.table(cov.file, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')


#variables <- c("cytocomp", "z_ifny_0_log_w", "z_il6_0_log_w",  "z_tnfa_0_log_w", "ifny", "il6", "tnfa", "ISEL_Mean", "PSS_all_mean", "HS_CRP", "LogCRP")
#variable_names <- c("cytocomp", "z_ifny_0_log_w", "z_il6_0_log_w", "z_tnfa_0_log_w", "ifny", "il6", "tnfa", "Social Support", "Psychological Stress", "High Sensitivity CRP", "LogCRP")

## grab PFAS variable names form deseq_nocombat.R script
variables <- PFASvars[-21]
variable_names <- PFASvars[-21]

#####

	for(i in 1:length(variables)){
	    myvar <- variables[i]
	    var_name <- variable_names[i]

	    cat("##########################################################\n")
	    cat("Now testing distribution of", var_name, "\n")
	    cat("#################\n")

	    #fname <- paste(iddir, platePrefix, "C0", ".sampleid_", myvar, plateSufix, ".txt", sep="")
	    #samplesV <- read.table(fname, sep="\t", header=T, quote='"', comment="")

	    # Read in the data

	    #cov <- covariates %>% filter(dbgap.ID %in% samplesV$Sample_ID)

	    cov <- covariates

	    cov <- cov[!is.na(cov[, myvar]),]
	    cov$variable <- myvar

		# Function to insert newlines into long titles
		wrap_title <- function(title, width = 55) {
		  paste(strwrap(title, width = width), collapse = "\n")
		}

	    #fname=paste(myvar,"_ALOFT_SCAIP_batch_2_3_removed_","histogram",".pdf",sep="")
	    #pdf(fname)
	    p1 <- ggplot(cov, aes(x=get(myvar)))+#, y=get(myvar))) +
	        geom_histogram(color = "black", fill = "dodgerblue4", alpha=0.5) +
	        #stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="firebrick", fill="firebrick") +
  			#ggtitle(wrap_title(paste0(myvar, ": ", var_name), width = 55)) +  # Apply the wrapping function here
  			ggtitle(wrap_title(paste0(var_name), width = 55)) +  # Apply the wrapping function here
	        theme(plot.title = element_text(hjust = 0.5, size=20), 
	              axis.text.y = element_text(size = 20), 
	              axis.text.x = element_text(size = 20), 
	              axis.title = element_text(size = 20),
	              plot.subtitle = element_text(size = 18, hjust = 0.5)) +
	        #theme(legend.position="none", axis.title.x = element_blank()) +
	        labs(x=paste0(myvar), subtitle = paste0("n=",nrow(cov))) 
	        
	    #print(p1)
	    #dev.off()


	    figfn <- paste(getwd(),"/", myvar,"_histogram_distribution",".png",sep="")
	    png(figfn, width=1000, height=1000, res=120)
	    #print(plot_grid(fig1, fig3, ncol=2, labels="AUTO", label_fontface="plain"))
	    print(p1)
	    dev.off()

	}




