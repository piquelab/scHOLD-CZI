
library(foreign)
library(tidyr)
library(plyr)
library(dplyr)
library(haven) 

#####################################################################
###### script for making covariate file (restructurring)     ########
######  from the spss file                                   ########
###### ommited from using other files to select samples      
###### /rs/rs_grp_scaloft/variables/make_cov/across_waves_variables_AR_01_2024.txt
###### script modified from Justyna's: /nfs/rprdata/ALOFT/AL1-6_ln/covariates/merge_with_Slatcher_cv.R
#####################################################################

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/rs/rs_grp_schold/CZI/","HOLD MASTERFILE 3.27.24_FL_final.sav","6055") #for testing
basefolder=args[1]
savfile=args[2]

# load new covariate file from Sam's group
#slatch <- read.spss("/rs/rs_grp_scaloft/variables/Luca _lab_data 9.14.23.sav",to.data.frame = TRUE)
slatch <- read_sav(paste0(basefolder,savfile)
data <- as.data.frame(slatch)

if(!is.na(args[3])){
samplestoremove=args[3]
data <-subset(data, !pID %in% samplestoremove)
}

# all -99 to NA
sum(data == -99, na.rm = TRUE)
cnts <- sapply(data, function(x) sum(x == -99, na.rm = TRUE)) #
cnts <- cnts[cnts > 0] # 3 PCs and cage for wave b2 have -99. covert all to NA

data[] <- lapply(data, function(x) replace(x, x == -99, NA))

# select columns to decide which participantIDWave should be exluded. infer they did not complete the study. 
# selected 9 variables to decide which sampleID waves to drop. 
columns <- c("cage1", "csex1", "ceth1", "cwght1", "chght1", "FEVPP", "ctasfq", "lymp_av", "pedu")
data$valid_basic <- rowSums(!is.na(data[, columns]))
summary(data$valid_basic)
dffilt <- data %>% filter(valid_basic > 0)


# get dbgap.ID and check if all samples that we want are included: 
# load the samples and waves to make the covariate files for. (can also load an old covariate file)
cv <- read.table("/rs/rs_grp_scaloft/variables/ALOFT_Sampels_all_per_project_cohort_corrected_01_24.txt", header=T)
#cv <- cv %>% select(dbgap.ID, Participant.ID, Wave, ChildWaveId)

## re-assign wave, and then sample.ID, ChildWaveID
cv$wave_old <- cv$Wave
cv$Wave[cv$Wave == "A"] <- "A1"
cv$Wave[cv$Wave == "B"] <- "B1"
cv$Wave[cv$Wave == "C"] <- "C1"
cv$Wave[cv$Wave == "D"] <- "D1"
cv$Wave[cv$Wave == "a"] <- "a2"
cv$Wave[cv$Wave == "b"] <- "b2"
cv$Wave[cv$Wave == "c"] <- "c2"
table(cv$Wave)
cv$ParticipantIDwave <- paste0(cv$Participant.ID, "-", cv$Wave)

cv %>% filter(!(ParticipantIDwave %in% dffilt$ParticipantIDwave))
# 7 wave D samples. 2 wave c2 samples. 1 wave C1 samples. don't include them. 


cv$childID <- cv$Participant.ID

colcv <-  cv[,c("dbgap.ID", "childID")]#, "wave_old")]#, "puberty", "RNA_puberty", "eLife", "SCAIP1_6", "SCAIP7_18")]
cvunique <- colcv[!duplicated(colcv$dbgap.ID),]

dffilt <- dffilt %>% left_join(cvunique, by = "childID")
dffilt$sample.ID <- paste0(dffilt$dbgap.ID, "-", dffilt$Wave)

dffilt <- dffilt %>% left_join(cv[c("ParticipantIDwave", "wave_old")], by = "ParticipantIDwave")

#dffilt %>% filter(is.na(dbgap.ID)) %>% select(childID, Wave, dbgap.ID, valid_basic)
missing <- dffilt %>% filter(is.na(dbgap.ID)) %>% select(childID, Wave)
write.table(missing, "/rs/rs_grp_scaloft/scALOFT_2024/covariates/make_cov/Samples_valid_ALOFT_variables.txt", sep="\t", row.names=F, quote=F)
# 33 samples have missing dbgap.ID. no waves from any of these samples has been used in the luca lab before
# these two have PBMC and DNA, check if enough variables are valid for these:
#dffilt %>% filter(childID == "3317") %>% select(valid_basic, ParticipantIDwave)
#dffilt %>% filter(childID == "3292") %>% select(valid_basic, ParticipantIDwave)

data <- dffilt %>% filter(!is.na(dbgap.ID))
length(unique(data$childID))
length(unique(data$ParticipantIDwave))


# add which project the data is included in:
cv$sample.ID <-  paste0(cv$dbgap.ID, "-", cv$Wave)
cvproj <-  cv[,c("sample.ID", "puberty", "RNA_puberty", "eLife", "SCAIP1_6", "SCAIP7_18")]
data1 <- data %>% left_join(cvproj, by = "sample.ID")


write.table(data1, "/rs/rs_grp_scaloft/scALOFT_2024/covariates/make_cov/ALOFT_covariate_resuttructured_not_fixed_unq265_n629_03_30_2024.txt", sep="\t", row.names=F, quote=F)















################################################### delete #############################################
#length(unique(dffilt$childID))
#length(unique(cv$dbgap.ID))


#### check if the 33 are present in old covariate files: 
cv <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_puberty_4-1-2021.txt",sep="\t", stringsAsFactors=F, header=T, comment="", quote='"')
cv <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_LSI.txt",sep="\t", stringsAsFactors=F, header=T, comment="", quote='"')
ids <- dffilt %>% filter(is.na(dbgap.ID)) %>% pull(childID) #, Wave, dbgap.ID)
ids <- as.vector(ids)

cv %>% filter(Participant.ID %in% ids)

data %>% filter(childID %in% ids)


dffilt <- merge(dffilt, colcv, by="childID", all.x=T)
dffilt <- dffilt[order(dffilt$ParticipantIDwave), ]


ids <- dffilt %>% filter(is.na(dbgap.ID)) %>% pull(childID) #, Wave, dbgap.ID)
idwave <- dffilt %>% filter(is.na(dbgap.ID)) %>% pull(ParticipantIDwave) #, Wave, dbgap.ID)

cv %>% filter(ParticipantIDwave %in% idwave)

cv %>% filter(Participant.ID %in% ids)
cov %>% filter(Participant.ID %in% ids)



dffilt %>% filter(childID %in% ids) %>% select(childID, Wave, dbgap.ID, ParticipantIDwave)


sum(is.na(dffilt$dbgap.ID))

dffilt %>% filter(childID == "3178")

df %>% filter(childID == "3178")
df %>% filter(childID %in% ids) %>% select(childID, Wave)#, dbgap.ID) #, ParticipantIDwave)




dffilt %>% filter(is.na(cage1)) %>% pull(childID, ParticipantIDwave, cage1, csex1, cwght1, FEVPP)
ids <- dffilt %>% filter(is.na(cage1)) %>% pull(ParticipantIDwave)
cv %>% filter(ParticipantIDwave %in% ids)

sum(is.na(dffilt$csex1))


# figure out which childwaveIDs should be kept
df$valid_count <- rowSums(!is.na(df))
dffilt <- df %>% filter(valid_count > 10)
# 6 of the columns are family and child IDs. so this filters for > 4 valid columns



# load covariate file: 
cov.file <- "/rs/rs_grp_scaloft/scALOFT_2024/covariates/covariate_tables/ALOFT_covariate_fixed_age_PBMC_IFNGswap_IL13_dbgapW_height-weight_n521_uniq-n265_03_18_2024.txt"
cov <- read.table(cov.file, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')


# load the bulk covariate file:
cov <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_puberty_4-1-2021.txt",sep="\t", stringsAsFactors=F, header=T, comment="", quote='"')
cov %>% filter(Participant.ID %in% ids) %>% select(Participant.ID, Wave, dbgap.ID, ChildWaveId)



cv %>% filter(ParticipantIDwave == "3299-a2")

df %>% filter(ParticipantIDwave == "3219-c2")
df %>% filter(ParticipantIDwave == "3222-c2")



df %>% filter(valid_basic == 1)


df$valid_count <- rowSums(!is.na(df))


dim(df)
sum(!is.na(df$cage1))
sum(!is.na(df$csex1))
sum(!is.na(df$ceth1))
sum(!is.na(df$ceth1))
sum(!is.na(df$cwght1))
sum(!is.na(df$chght1))
sum(!is.na(df$FEVPP))


df$valid_count <- rowSums(!is.na(df))
dfsub <- df %>% filter(valid_count < 10)
dim(dfsub)


df %>% filter(valid_count == 14)

dfsub <- df %>% filter(valid_count > 10)
sum(is.na(dfsub$cage1))

dfsub %>% filter(is.na(cage1))


df %>% select(childID, Wave, csslat...68)

length(intersect(colnames(df2), colnames(df1)))
length(intersect(colnames(df0), colnames(df1)))
setdiff(colnames(df1), colnames(df0))
setdiff(colnames(df2), colnames(df0))
setdiff(colnames(df2), colnames(df1))


