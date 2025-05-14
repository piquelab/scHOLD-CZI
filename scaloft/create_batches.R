library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
set.seed(905)

#genotype issues:
cd /rs/rs_grp_scaloft/locoglimpse2/ALOFT2only/
cat plink2.kin0 | awk '$6>0.05'

dat_e <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_samples_batchassignment_with5new.txt")
kin <- fread("/rs/rs_grp_scaloft/locoglimpse2/ALOFT2only/plink2.kin0")
names(kin)[1]<-"ID1"
idmatch <- fread("/rs/rs_grp_scaloft/locoglimpse2/aloft2b1.v2.ids.match.txt") #/rs/rs_grp_scaloft/locoglimpse2/aloft2b1.ids.txt

almatch <- merge(idmatch,dat_e,by.x="Participant.ID",by.y="c1pid")
famc <- plyr::count(almatch,"fam_id")
fams <- subset(famc,freq>1)
almatch_fam <- subset(almatch, fam_id %in% fams$fam_id)
famcombo <- ddply(almatch_fam,"fam_id",plyr::summarize,
  combo1=paste0(dbgap.ID[1],"_",dbgap.ID[2]),
  combo2=paste0(dbgap.ID[2],"_",dbgap.ID[1]))
almatch_notfam <- subset(almatch, !fam_id %in% fams$fam_id)

kin$combo <- paste0(kin$ID1,"_",kin$ID2)
kin_famm1 <- merge(kin, famcombo[,-3], by.x="combo",by.y="combo1")
kin_famm2 <- merge(kin, famcombo[,-2], by.x="combo",by.y="combo2")
kinfam <- rbind(kin_famm1,kin_famm2)
kinfam_kin05 <- subset(kinfam,KINSHIP>0.05)
kinfam_notkin05 <- subset(kinfam,KINSHIP<=0.05)

kinfam_kin05[order(kinfam_kin05$fam_id),]
kinfam_kin05[order(kinfam_kin05$KINSHIP),]

kinnotfam <- subset(kin, !combo %in% kinfam$combo)
kinnotfamm1 <- merge(almatch,kinnotfam[,-2],by.x="dbgap.ID",by.y="ID1")
kinnotfamm2 <- merge(almatch,kinnotfam[,-1],by.x="dbgap.ID",by.y="ID2")
kinnotfamm <- rbind(kinnotfamm1,kinnotfamm2)

kinnotfam_kin05 <- subset(kinnotfamm,KINSHIP>0.05)
kinnotfam_kin05[order(kinnotfam_kin05$KINSHIP),]

idmatch[idmatch$V2=="AL-380",]

halfsib <- subset(kinfam_kin05, KINSHIP<0.2)

#two samples in pathway B (Kinship reflects siblings, but they have different family IDs) are 3378 and 3371, we may end up removing one so treat like the low viability samples, also assign a later batch in case we make a different decision
#move 3351 3352 , as well as 4352 to batch 8 as we are waiting on resequencing 
#move half siblings (k>0.05 <0.2) to later batches

#now remaking batches based on info we got from genotyping
###########################################################
data_o <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_samples.txt")
#REMOVE samples missing vials, low viability etc
data_o <- transform(data_o, fam_id = gsub("^[34]","",c1pid))
data <- subset(data_o, useable=="Y")
data_nodna <- subset(data_o, useable=="no_dna" | useable=="no_sample")
data_half <- subset(data_o, fam_id %in% halfsib$fam_id)
data_unknownsibs <- subset(data_o, c1pid %in% c("3378","3371"))
data_problemsibs <- subset(data_o, c1pid %in% c("3351","3352","4352"))
#subsetting individuals to those we know are fine to proceed with 
earlybatch <- subset(data, !c1pid %in% c(data_half$c1pid,data_unknownsibs$c1pid,data_problemsibs$c1pid))
#these are separate so they can be spread across these early batches
data_lowviability <- subset(data_o, useable=="low_viability" & !c1pid %in% c(data_half$c1pid,data_unknownsibs$c1pid,data_problemsibs$c1pid))
latebatch_low <- subset(latebatch, useable=="low_viability")

datsib <- subset(earlybatch, asthmatic==0)
dat <- earlybatch
num_batches <- ceiling(nrow(dat)/12)
dat$group <- rep_len(1:num_batches, length.out = nrow(dat)) #if you wanted to make 12 groups
dat$group <- as.factor(dat$group)
a <- 0.1; b <- 0.1; c <- 0.1; d <- 0.1
thresh <- 0.84 #Minimum threshold for p value. original 0.85
z <- 1
#while (a < thresh | b < thresh |c < thresh |d < thresh) {
while (a < thresh | b < thresh |c < thresh) {
  dat <- transform(dat, group = sample(group)) #shuffles the groups
  x <- summary(aov(age ~ group, dat)) #ANOVA for continuous variables
  a <- x[[1]]['group','Pr(>F)']
  x <- summary(table(dat$group, dat$sex)) #Chi Sq for categorical variables
  b <- x[['p.value']]
  x <- summary(table(dat$group, dat$race))
  c <- x[['p.value']]
  #x <- summary(table(dat$group, dat$asthmatic))
  #d <- x[['p.value']]
  z <- z + 1
  if (z > 10000) {
    print('10,000 tries, no solution, reduce threshold')
    break
  }
}
dat_s <- subset(dat, asthmatic==1)
#match sibling into group
datsib_m <- merge(datsib,dat_s[,c("fam_id","group")],by="fam_id")
dat_b <- rbind(dat_s,datsib_m)
#the previous will not be completely even groups so need to even them out
table(dat_b$group)
 1  2  3  4  5  6
10  8 11 10 11  9

g3 <-dat_b[dat_b$group=="6",]
g3[order(g3$fam_id),]
dat_c <-transform(dat_b, group=ifelse(fam_id=="348" | fam_id=="391","2",ifelse(fam_id=="323","1",ifelse(fam_id=="355" | fam_id=="396","4",group)))) #make sure 1-5 are full as only 4 low viability to add in
table(dat_c$group)
 1  2  3  4  5  6
12 10 11 12 11  3

#add in low viability samples, spread out
dat_c[dat_c$fam_id %in% data_lowviability$fam_id,] # no matching fams to make sure get put together
data_lowviability[order(data_lowviability$fam_id),]
data_lowviability <- transform(data_lowviability, group=ifelse(fam_id=="327","3",ifelse(fam_id=="350","2",ifelse(fam_id=="338","5",NA))))
dat_d <- rbind(data_lowviability,dat_c)
table(dat_d$group)
 1  2  3  4  5  6
12 12 12 12 12  3
#dat_d %>% group_by(group) %>% summarize(age = mean(age, na.rm = TRUE),sex=mean(sex, na.rm = TRUE),race=mean(race, na.rm = TRUE),asthmatic=mean(asthmatic, na.rm = TRUE))

#move 1-2 14yr from group 6 to group2
#dat_d[dat_d$group=="6",]  #347 to g2, 344 to g6
#dat_d <- transform(dat_d, group=ifelse(fam_id=="388"| fam_id=="396","2",ifelse(fam_id=="325","6",group)))

fwrite(dat_d, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_samples_batchassignment.txt")

#dat_d <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_samples_batchassignment.txt")
#dat_d <- transform(dat_d, group=as.factor(group))

#now need to add new samples
data_half <- transform(data_half, note="halfsib")
data_unknownsibs <- transform(data_unknownsibs, note="unknown_sibs")
data_problemsibs <- transform(data_problemsibs, note="problem_sibs")
v25 <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_samples_5new.txt")
v25 <- transform(v25, fam_id = gsub("^[34]","",c1pid))
latebatch <- rbind(data_half,data_unknownsibs,data_problemsibs,v25,dat_d[dat_d$group==6,-9])

datsib <- subset(latebatch, asthmatic==0)
dat <- latebatch
num_batches <- ceiling(nrow(dat)/12)
dat$group <- rep_len(1:num_batches, length.out = nrow(dat)) #if you wanted to make 12 groups
dat$group <- as.factor(dat$group)
a <- 0.1; b <- 0.1; c <- 0.1; d <- 0.1
thresh <- 0.57 #Minimum threshold for p value. original 0.85
z <- 1
#while (a < thresh | b < thresh |c < thresh |d < thresh) {
while (a < thresh | b < thresh |c < thresh) {
  dat <- transform(dat, group = sample(group)) #shuffles the groups
  x <- summary(aov(age ~ group, dat)) #ANOVA for continuous variables
  a <- x[[1]]['group','Pr(>F)']
  x <- summary(table(dat$group, dat$sex)) #Chi Sq for categorical variables
  b <- x[['p.value']]
  x <- summary(table(dat$group, dat$race))
  c <- x[['p.value']]
  #x <- summary(table(dat$group, dat$asthmatic))
  #d <- x[['p.value']]
  z <- z + 1
  if (z > 10000) {
    print('10,000 tries, no solution, reduce threshold')
    break
  }
}
dat_s <- subset(dat, asthmatic==1)
#match sibling into group
datsib_m <- merge(datsib,dat_s[,c("fam_id","group")],by="fam_id")
dat_b <- rbind(dat_s,datsib_m)
#the previous will not be completely even groups so need to even them out
table(dat_b$group)
dat_b[order(dat_b$group),] #unknown sibs correctly in diff groups, problemsibs in same group. low viability not all in 1 group 

g3 <-dat_b[dat_b$group=="3",]
sample_n(g3,5)
dat_c <-transform(dat_b, group=ifelse(c1pid=="3351" ,"1",ifelse(c1pid %in% c("3403","3400","4400","3402","3360"),"2",group))) 
dat_late <-transform(dat_c, group=as.numeric(group)+5)
table(dat_late$group)

dat_e <- rbind(dat_d[!dat_d$group=="6",],dat_late)
table(dat_e$group)

png(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_","age","_hist.png"), width=1000, height=1000, res=120)
ggplot(dat_e, aes(age, fill = group)) + geom_histogram()+ facet_wrap(~group)
dev.off()
png(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_","sex","_hist.png"), width=1000, height=1000, res=120)
ggplot(dat_e, aes(sex, fill = group)) + geom_histogram()+ facet_wrap(~group)
dev.off()
png(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_","race","_hist.png"), width=1000, height=1000, res=120)
ggplot(dat_e, aes(race, fill = group)) + geom_histogram()+ facet_wrap(~group)
dev.off()
png(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_","asthmatic","_hist.png"), width=1000, height=1000, res=120)
ggplot(dat_e, aes(asthmatic, fill = group)) + geom_histogram()+ facet_wrap(~group)
dev.off()

png(paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_","fam_id","_bar.png"), width=1000, height=1000, res=120)
ggplot(dat_e, aes(fam_id, fill = group)) + geom_bar(position = 'stack')
dev.off()

almatch_new <- merge(idmatch,dat_e,by.x="Participant.ID",by.y="c1pid")
almatch_new[almatch_new$Participant.ID %in% c(data_lowviability$c1pid, data_unknownsibs$c1pid, data_problemsibs$c1pid,latebatch_low$c1pid),]

fwrite(almatch_new, sep='\t', quote=F, row.names=F, col.names=T, file="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_samples_batchassignment_040325.txt")

##
library(glmnet)
almatch_new <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/v2_samples_batchassignment_040325.txt")
almatch_new <- transform(almatch_new, group=as.factor(group),sex=as.factor(sex),race=as.factor(race))
x = model.matrix( as.formula(~group+sex+age+race) , almatch_new)
qr(x)$rank
ncol(x) 
t(x) %*% x