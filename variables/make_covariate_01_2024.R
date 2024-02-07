
library(foreign)
library(tidyr)
library(dplyr)
library(haven) 

#####################################################################
###### script for making covariate file (restructurring)     ########
######  from the spss file                                   ########
###### uses a supp file witih all the variables of interest  ########
###### /rs/rs_grp_scaloft/variables/make_cov/across_waves_variables_AR_01_2024.txt
###### script modified from Justyna's: /nfs/rprdata/ALOFT/AL1-6_ln/covariates/merge_with_Slatcher_cv.R
#####################################################################

rm(list=ls())

setwd("/rs/rs_grp_scaloft/variables/make_cov/")

# load new covariate file from Sam's group
#slatch <- read.spss("/rs/rs_grp_scaloft/variables/Luca _lab_data 9.14.23.sav",to.data.frame = TRUE)
slatch <- read_sav("/rs/rs_grp_scaloft/variables/Luca _lab_data 9.14.23.sav")

# load the samples and waves to make the covariate files for. (can also load an old covariate file)
cv <- read.table("/rs/rs_grp_scaloft/variables/ALOFT_Sampels_all_per_project_cohort_corrected_01_24.txt", header=T)
cv <- cv %>% select(dbgap.ID, Participant.ID, Wave, ChildWaveId)

# load the rest of the variables we are interested in. This table should be modified if additional variabels need to be included and analyzed. 
acrossWaves <- read.table("/rs/rs_grp_scaloft/variables/make_cov/across_waves_variables_AR_01_2024.txt", sep="\t", stringsAsFactors=F, header=T,comment="", quote='"')

#load W 0 only variables. Out of the following list, we are only interested in csex and ceth from w 0 only
WAonly_vars <- c("ceth1_0", "csex1_0")
slatchWA <- select(slatch, c("childID",WAonly_vars))

# add the W1 only variables to covariate:
cv <- left_join(cv, slatchWA, by=c("Participant.ID"="childID"))

# add appropriate wave to the sample for across-wave variables:
# split by wave:
cvA <- filter(cv, Wave=="a" | Wave=="A")
cvB <- filter(cv, Wave=="b" | Wave=="B")
cvC <- filter(cv, Wave=="c" | Wave=="C")
cvD <- filter(cv, Wave=="D")

#merge appropriate columns:
cvAm <- left_join(cvA, slatch[,c("childID",acrossWaves$W1.var.name)], by=c("Participant.ID"="childID"))
cvBm <- left_join(cvB, slatch[,c("childID",acrossWaves$W3.var)], by=c("Participant.ID"="childID"))
cvCm <- left_join(cvC, slatch[,c("childID",acrossWaves$w4.var)], by=c("Participant.ID"="childID"))

# rename the added columns to universal names:
colnames(cvAm)[(ncol(cvA)+1):ncol(cvAm)] <- acrossWaves$W3.var.name
colnames(cvBm)[(ncol(cvB)+1):ncol(cvBm)] <- acrossWaves$W3.var.name
colnames(cvCm)[(ncol(cvC)+1):ncol(cvCm)] <- acrossWaves$W3.var.name
#add those columns to D:
cvDm <- cbind(cvD, matrix(nrow=nrow(cvD), ncol=ncol(cvAm)-ncol(cvA)))
colnames(cvDm)[(ncol(cvD)+1):ncol(cvDm)] <- acrossWaves$W3.var.name

# merge everything back together:
identical(colnames(cvA),colnames(cvB))
identical(colnames(cvC),colnames(cvD))
cvm <- rbind(cvAm,cvBm,cvCm,cvDm)

# add a column indicating missing basic covariate (age, wght, :
cvm[,"missing_basic_covariate"] <- "no"
cvm[!((!is.na(cvm$cwght1))*(!is.na(cvm$chght1))),"missing_basic_covariate"] <- "yes"
cvm[is.na(cvm$cage1),"missing_basic_covariate"] <- "yes"
# just checking that none miss other basic covariate if not missing wght,chght:
cvnonmiss <- filter(cvm,missing_basic_covariate=="no")
(sum(is.na(cvnonmiss$csex1_0))==0)

############## with using the heaven package and the read_sav function, i don't need to do this conversion from factor level anymore. The numeric value now matches the master data description questionare. 
# some columns are factor levels and can be converted to numeric
#factor_columns <- sapply(cvm, function(x) is.factor(x))
#num_factor_columns <- sum(factor_columns)
#cat("Number of factor columns:", num_factor_columns, "\n") 
# Number of factor columns: 6
#cat("Factor columns:", names(factor_columns[factor_columns]), "\n")
# Factor columns: ceth1_0 csex1_0 pincme pedu psesl cgpd5


# Assuming 'your_data_frame' is the name of your data frame
# Assuming 'factor_columns' is a logical vector indicating factor columns

# Convert factor columns to numeric
#cvm[, factor_columns] <- lapply(cvm[, factor_columns], as.numeric)

# save the covariate file:
#write.table(cvm, "/rs/rs_grp_scaloft/variables/make_cov/ALOFT_covariate_resuttructured_not_fixed_01_2024.txt", sep="\t", row.names=F, quote=F)
write.table(cvm, "/rs/rs_grp_scaloft/variables/make_cov/ALOFT_covariate_resuttructured_not_fixed_01_09_2024.txt", sep="\t", row.names=F, quote=F)



