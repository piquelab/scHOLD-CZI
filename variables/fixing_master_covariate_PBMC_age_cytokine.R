library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)

rm(list=ls())

##############################################################################
######## manually fix the master covariate for identified issues  ############
######## age -99, PBMC wave b (and any other), cytokine 99/-99    ############
########
########
##############################################################################

setwd("/rs/rs_grp_scaloft/variables/make_cov/")

#covariates <- read.csv("/rs/rs_grp_scaloft/variables/covariates_new_all_scALOFT_SCAIP1_18_and_other_aloft_samples_bulk_n508.csv")
#cov1 <- read.table("/rs/rs_grp_scaloft/variables/make_cov/ALOFT_covariate_resuttructured_not_fixed_01_2024.txt", header=T, sep="\t")

cov <- read.table("/rs/rs_grp_scaloft/variables/make_cov/ALOFT_covariate_resuttructured_not_fixed_01_09_2024.txt", header=T, sep="\t")

#setdiff(colnames(cov1), colnames(cov2))
#setdiff(cov1$dbgap.ID, cov2$dbgap.ID)
#identical(cov1$ChildWaveId, cov2$ChildWaveId)
#identical(cov1$csex1_0, cov2$csex1_0)
#identical(cov1$chght1, cov2$chght1)
#setdiff(cov1$ceth1_0, cov2$ceth1_0)


############################## add genotype PCs #############################################################
## add genotype PCs. Both old (for ALOFT) and new genotype PCs (Ones I generated for SCAIP7-18)

# newly generated genotype pcs
allgpc <- read.table("/rs/rs_grp_scaloft/genotypes/ALOFT_n258L_genotype_PCs1-3_filtered.txt", header=TRUE)

merger <- data.frame(t(allgpc)[-1,])
merger$dbgap.ID <- rownames(merger)
merger$dbgap.ID <- gsub('\\.', '-', merger$dbgap.ID)
colnames(merger)[1:3] <- c("genPC1","genPC2","genPC3")

length(intersect(unique(cov$dbgap.ID), merger$dbgap.ID))
merged_cov <- merge(cov, merger, by="dbgap.ID", all.x=T)
merged_cov$genPC1 <- as.numeric(as.character(merged_cov$genPC1))
merged_cov$genPC2 <- as.numeric(as.character(merged_cov$genPC2))
merged_cov$genPC3 <- as.numeric(as.character(merged_cov$genPC3))

##### this should be used to get SCAIP1-6 old genotype PCs: 

#Covariates with LSI variables added:
cov.file <- "/nfs/rprdata/SCAIP-ALOFT-Shreya/SCAIP-ALOFT/covariates/SCAIP_ALOFT_covariate_all_treatments_LSI_added_8-11-22.txt"
covscold <- read.table(cov.file, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')
length(unique(covscold$dbgap.ID))
covscold <- covscold %>% filter(!(duplicated(dbgap.ID)))

oldvcg <- covscold %>% select(dbgap.ID, SCAIP1_6_genPC1, SCAIP1_6_genPC2, SCAIP1_6_genPC3) %>% distinct(dbgap.ID, .keep_all=T)

mm_cov <- merge(merged_cov, oldvcg, by="dbgap.ID", all.x=T)
cov <- mm_cov

###################################### delete this. the genePCs in this file dont match with other SCIAP1-6 gene pcs. 
# add SCAIP7-18 genotypes
#oldcv <- read.table("/nfs/rprdata/SCAIP-ALOFT-Shreya/SCAIP-ALOFT/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_Wave-a-added_4-8-2021.txt", sep="\t", header=T, quote='"', comment="")
#oldcv <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/covariates/ALOFT_RNA-seq_sample_masterfile_corrected_merged_8-23-19_puberty_4-1-2021_updated_astma-added.txt",sep="\t", stringsAsFactors=F, header=T, comment="", quote='"')
#length(unique(oldcv$dbgap.ID))
#oldvcg <- oldcv %>% select(dbgap.ID, genPC1, genPC2, genPC3) %>% distinct(dbgap.ID, .keep_all=T)
#colnames(oldvcg) <- c("dbgap.ID", "SCAIP1_genPC1", "SCAIP1_genPC2", "SCAIP1_genPC3")

#mm_cov <- merge(merged_cov, oldvcg, by="dbgap.ID", all.x=T)
#cov <- mm_cov

########################################## Fix age issue #############################################################
cov %>% filter(cage1<0) %>% select(dbgap.ID, ChildWaveId, cage1) # 5 samples in wave b with -99
cov$cage1[cov$cage1 == "-99"] <- NA


# assign age of missing samples using the provided visit dates. Samples from scALOFT with missing age that will need calculation: 3227-b, 357-C, and 3249-b
csna <- c("3227-b", "357-C", "3249-b")
csna_prior <- c("3227-a", "357-B", "357-A", "3249-a")
#slatch %>% filter(childID == "357") %>% select(childID, cohort, cage1_0, cage1_1, cage1_2)

cov %>% filter(ChildWaveId %in% csna) %>% select(dbgap.ID, ChildWaveId, cage1) # 5 samples in wave b with -99
cov %>% filter(ChildWaveId %in% csna_prior) %>% select(dbgap.ID, ChildWaveId, cage1) # 5 samples in wave b with -99

#ID     w0         w1               w2
#3227   1/25/2017  3/7/2018         NA*(7/19/2019)
#357    12/5/2011  NA*(12/10/2012)  3/17/2014
#3249   4/12/2017  5/2/2018         NA*(5/23/2019)

#  childID cohort  cage1_0    cage1_1    cage1_2
#  357     1       13.37988   14.45      NA

# dbgap.ID  ChildWaveId    cage1
#   AL-181  3227-a         12.83231
#   AL-199  3249-a         11.78097

#calculate age for 357-C
initial_age <- 13.37988
date_data1 <- as.Date("2011-12-05")
date_target <- as.Date("2014-03-17")
age_on_target_date <- initial_age + as.numeric(difftime(date_target, date_data1, units = "days")) / 365.25
print(age_on_target_date)
cov$cage1[cov$ChildWaveId == "357-C"] <- age_on_target_date

#calculate age for 3227-b
initial_age <- 12.83231
date_data1 <- as.Date("2017-01-25")
date_target <- as.Date("2018-03-07")
age_on_target_date <- initial_age + as.numeric(difftime(date_target, date_data1, units = "days")) / 365.25
print(age_on_target_date)
cov$cage1[cov$ChildWaveId == "3227-b"] <- age_on_target_date

#calculate age for 3249-b
initial_age <- 11.78097
date_data1 <- as.Date("2017-04-12")
date_target <- as.Date("2018-05-02")
age_on_target_date <- initial_age + as.numeric(difftime(date_target, date_data1, units = "days")) / 365.25
print(age_on_target_date)
cov$cage1[cov$ChildWaveId == "3249-b"] <- age_on_target_date

########################################## Fix cytokine -99/99 issue #############################################################

cov %>% filter(IL5_co==99) %>% select(dbgap.ID, ChildWaveId, IL5_co)

cols <- c("IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc")

#cov <- cov %>% filter_at(vars(cols), any_vars(. == 99)) %>% 
 # 		mutate_at(vars(cols), ~ifelse(. == 99, NA, .))  #%>%
  		#select(dbgap.ID, ChildWaveId, all_of(cols))

cov <- cov %>%
  		mutate(across(all_of(cols), ~ifelse(. == 99, NA, .)))

#two samples at wave B had 99 for multiple cytokine variable. they were changed to NA.
#AL-059       384-B  
#AL-098      3126-B  

########################################## Fix PBMC fraction issue #############################################################

cov %>% filter(lymp_av>1) %>% select(dbgap.ID, ChildWaveId, lymp_av)
cov %>% filter(baso_av>1) %>% select(dbgap.ID, ChildWaveId, baso_av)

cov %>% filter(Wave=="b") %>% select(dbgap.ID, ChildWaveId, lymp_av) #%>% summary() ### all of samles in b wave are fraction not percentage. devide them by 10
cov %>% filter(Wave=="b") %>% select(dbgap.ID, ChildWaveId, mono_av) ### all of samles in b wave are fraction not percentage. devide them by 10
cov %>% filter(Wave=="b") %>% select(dbgap.ID, ChildWaveId, eosi_av) ### all of samles in b wave are fraction not percentage. devide them by 10
cov %>% filter(Wave=="b") %>% select(dbgap.ID, ChildWaveId, neut_av) ### all of samles in b wave are fraction not percentage. devide them by 10
cov %>% filter(Wave=="b") %>% select(dbgap.ID, ChildWaveId, baso_av) ### all of samles in b wave are fraction not percentage. devide them by 10

cov <- cov %>% mutate(lymp_av= ifelse(Wave == "b", lymp_av/100, lymp_av))
cov <- cov %>% mutate(mono_av= ifelse(Wave == "b", mono_av/100, mono_av))
cov <- cov %>% mutate(eosi_av= ifelse(Wave == "b", eosi_av/100, eosi_av))
cov <- cov %>% mutate(neut_av= ifelse(Wave == "b", neut_av/100, neut_av))
cov <- cov %>% mutate(baso_av= ifelse(Wave == "b", baso_av/100, baso_av))


### found out some samples were assigned to wrong cohorot. they were fixed when generating the new cov. 

########################################## SWAP IFNG_co with IFNG_hc for a (b and c dont have cytokine levels) ########################################
#cov.file <- "/rs/rs_grp_scaloft/variables/make_cov/ALOFT_covariate_resutructured_fixed_age_PBMC-b_cytkine_n519_uniq-n265_01_2024.txt"
#cov <- read.table(cov.file, sep="\t", stringsAsFactors=F, header=T, comment="", quote = '"')

cov$IFNG_co_tmp <- cov$IFNG_co
cov$IFNG_hc_tmp <- cov$IFNG_hc

tail(cov %>% filter(Wave=="a") %>% select(dbgap.ID, ChildWaveId, IFNG_hc_tmp, IFNG_co_tmp)) #%>% summary() ### all of samles in b wave are fraction not percentage. devide them by 10

cov <- cov %>%  
    mutate(
        IFNG_co = if_else(Wave == "a", IFNG_hc_tmp, IFNG_co_tmp),
        IFNG_hc = if_else(Wave == "a", IFNG_co_tmp, IFNG_hc_tmp)
        )

tail(cov %>% filter(Wave=="a") %>% select(dbgap.ID, ChildWaveId, IFNG_hc, IFNG_co))


###### add additional columns: Sex

cov$Sex <- cov$csex1_0

# save the covariate file:
write.table(cov, "/rs/rs_grp_scaloft/variables/make_cov/ALOFT_covariate_resutructured_fixed_age_PBMC-b_cytkine_IFNGswap_age-added_n521_uniq-n265_01_09_2024.txt", sep="\t", row.names=F, quote=F)


