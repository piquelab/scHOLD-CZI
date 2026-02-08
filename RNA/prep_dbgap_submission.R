library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(data.table)
library(ggrepel)
library(qvalue)



fname = "/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt"
cov <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")


fname = "/rs/rs_grp_schold/covariates/dbgap/HOLD_library_metadata_prelim_n165_dbgapIDs_batch_09_17_2024.txt"
batches <- read.table(fname, sep = "\t", header = TRUE, quote = '"', comment.char = "")


df2a <- cov %>% select(dbgap.ID, sex_alph) %>% left_join(batches, by="dbgap.ID")
dbgap_2a <- df2a %>% mutate(SUBJECT_ID=dbgap.ID, CONSENT="1", SEX=sex_alph, POOLED_SUBJECT_ID=Batch) %>% select(SUBJECT_ID, CONSENT, SEX, POOLED_SUBJECT_ID)
dbgap_2a$POOLED_SUBJECT_ID <- sub("^HOLD0([1-9])$", "HOLD\\1", dbgap_2a$POOLED_SUBJECT_ID)

fname="/rs/rs_grp_schold/covariates/dbgap/submiussion/2a_dbGaP_SubjectDS.txt"
write.table(dbgap_2a, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



### 4a scRNA library and scATAC library
# rna libraries: 
fname="/.rs/volumes/rs_grp_schold/rs_vol_schold/9cf75849-f9b5-4d4a-923d-cd5651afac0a/CZI/RNA/counts_cellranger_2024-04-19/libList.txt"
rnalib <- read.table(fname, sep = "\t", header = FALSE, quote = '"', comment.char = "")
colnames(rnalib) <- "SAMPLE_ID"
df4a1 <- rnalib %>% mutate(BODY_SITE="peripheral blood", ANALYTE_TYPE="RNA", IS_TUMOR="NO") %>% filter(!grepl("DEX", SAMPLE_ID))

# atac libraries
fname="/.rs/volumes/rs_grp_pique20/rs_vol_pique20/122ad1b4-ec4d-4d9f-9154-d5805e1e1913/schold/ATAC/counts_cellranger_atac/libList.txt"
ataclib <- read.table(fname, sep = "\t", header = FALSE, quote = '"', comment.char = "")
colnames(ataclib) <- "SAMPLE_ID" 
df4a2 <- ataclib %>% mutate(BODY_SITE="peripheral blood", ANALYTE_TYPE="DNA", IS_TUMOR="NO") %>% filter(!grepl("DEX", SAMPLE_ID))

#genotype
df4a3 <- dbgap_2a %>% mutate(SAMPLE_ID=SUBJECT_ID, BODY_SITE="peripheral blood", ANALYTE_TYPE="DNA", IS_TUMOR="NO") %>% select(SAMPLE_ID, BODY_SITE, ANALYTE_TYPE, IS_TUMOR)


fname="/rs/rs_grp_schold/covariates/dbgap/submiussion/4a_SampleAttributes_DS_scRNAseq.txt"
write.table(df4a1, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

fname="/rs/rs_grp_schold/covariates/dbgap/submiussion/4a_SampleAttributes_DS_scATACseq.txt"
write.table(df4a2, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)


fname="/rs/rs_grp_schold/covariates/dbgap/submiussion/4a_SampleAttributes_DS_Genotype.txt"
write.table(df4a3, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

##### 5a covariates metadata
df5a <- cov %>% select(dbgap.ID, age, sex_alph, ISEL_Mean, PSS_all_mean, cytocomp, SES, PHip, PWaist, WHR)  #%>% mutate(SUBJECT_ID="dbgap.ID")

# load PCs
args <- c("/rs/rs_grp_schold/CZI/RNA/analysis/","/rs/rs_grp_schold/covariates/dbgap/HOLD_covariates_n165_dbgapIDs_updated_WHR_05_28_2025.txt","ALL","fastdemux",13,0.1) #for testing

base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
dimset=args[5]
resset=args[6]

eigenvec2_o <- fread(file=paste0(base,"genotypePCnokin/",project,".eigenvec_pc.txt")) #will use col PC1
genepc <- eigenvec2_o %>% as.data.frame() %>% select(Sample_ID, PC1, PC2) %>% mutate(dbgap.ID=Sample_ID)

df5a <- left_join(df5a, genepc, by="dbgap.ID")

dbgap_5a <- df5a %>% mutate(SUBJECT_ID=dbgap.ID, AGE=age, SEX=sex_alph) %>% 
			select(-Sample_ID, -dbgap.ID, -age, -sex_alph, -PHip, -PWaist) %>% 
			select(SUBJECT_ID, AGE, SEX, everything())


fname="/rs/rs_grp_schold/covariates/dbgap/submiussion/5a_SubjectPhenotypes_DS.txt"
write.table(dbgap_5a, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)



### 3a SUbJECT ID to SAMPLE ID

# RNA SUBJECT ID to sample ID
bitdfrna <- df4a1 %>% separate(SAMPLE_ID, into = c("POOLED_SUBJECT_ID", "Bio", "treat"),sep = "-",remove = FALSE)
df3a <- dbgap_2a %>% inner_join(bitdfrna, by = "POOLED_SUBJECT_ID") %>% select(SUBJECT_ID, SAMPLE_ID)


#ATAC
bitdfatac <- df4a2 %>% separate(SAMPLE_ID, into = c("POOLED_SUBJECT_ID", "Bio", "treat"),sep = "-",remove = FALSE)
df3a2 <- dbgap_2a %>% inner_join(bitdfatac, by = "POOLED_SUBJECT_ID") %>% select(SUBJECT_ID, SAMPLE_ID)

#Genoyotype
df3a3 <- dbgap_2a %>% mutate(SAMPLE_ID=SUBJECT_ID) %>% select(SUBJECT_ID, SAMPLE_ID)

dbgap3a <- rbind(df3a, df3a2, df3a3)

fname="/rs/rs_grp_schold/covariates/dbgap/submiussion/3a_dbGaP_SubjectSampleMappingDS.txt"
write.table(dbgap3a, file = fname, sep = "\t", append = F, quote = F, col.names = T, row.names = F)

