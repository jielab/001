setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, dplyr, tidyverse, crosswalkr, lubridate)
indir = "D:/data/ukb/phe/"
dates_bad=as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1999-01-01", "2037-07-07"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main PHE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vip <- read.table(paste0(indir,"common/ukb.vip.dat"), header=F, flush=T)
	vip$V1[duplicated(vip$V1)]; vip$V2[duplicated(vip$V2)] #check duplication
	vip$V2 <- paste0("f.", vip$V2, ".0.0")
source(paste0(indir,"raw-670287/vip.r")) # first use the latest version
	phe1 <- bd %>% subset(select=grep("f.eid|\\.0\\.0", names(bd))) %>% rename(eid=f.eid)
	names1 <- names(phe1)
source(paste0(indir,"raw-50136/vip.r")) # then add previous version
	phe2 <- bd %>% subset(select=!names(bd) %in% names1) 
	phe2 <- phe2 %>% subset(select=grep("f.eid|\\.0\\.0", names(phe2))) %>% rename(eid=f.eid)
phe0 <- merge(phe1, phe2, by="eid") %>% 
	crosswalkr::renamefrom(vip, V2, V1, drop_extra=F) %>% filter(eid>0)
	phe0[phe0 <0] = NA # BE CAREFUL
	dates <- names(phe0)[sapply(phe0, is.Date)]; summary(phe0[,dates]) # check invalid dates
	for (date1 in dates) { phe0[[date1]][phe0[[date1]] %in% dates_bad] <- NA }
phe <- phe0 %>% 
	mutate(
	ethnic_cat = ifelse(grepl("^1",ethnicity),1, ifelse(grepl("^2",ethnicity),2, ifelse(grepl("^3",ethnicity),3, ifelse(grepl("^4",ethnicity),4, ethnicity)))),
	ethnic_cat = factor(ethnic_cat, levels=c(1,2,3,4,5,6), labels=c("White","Mixed","Asian","Black","Chinese","Other")),	
	age_cat = cut(age, breaks=seq(35,75,5)),
	sex_f = ifelse(sex==0, "female", "male"),
	bmi_cat = cut(bmi, breaks=c(10,18.5,25,30,100), labels=c("lean","healthy","overweight","obese")),
	bmi_cat = factor(bmi_cat, levels=c("healthy","lean","overweight","obese")),
	smoke_status = factor(smoke_status, levels=0:2, labels=c("never","previous","current")), 
	alcohol_status = factor(alcohol_status, levels=0:2, labels=c("never","previous","current"))
	)
saveRDS(phe, file="D:/data/ukb/Rdata/ukb.phe.rds")
	# covid0 <- read.table("D:/data/ukb/hes/covid19_result_england.txt", header=T)[,c(1,2,5,6)]
	# covid <- aggregate(result ~ eid, data=covid0, sum) %>% rename(inf=result) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ICD 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(indir,"raw-670287/icd.r")); icd <- bd
source(paste0(indir,"raw-670287/icdDate.r")); icdDate <- bd
#icd <- icd[1:10, 1:15]; icdDate <- icdDate[1:10, 1:15]
dates <- names(icdDate)[sapply(icdDate, is.Date)]; summary(icdDate[,-1])
for (date1 in dates) { icdDate[[date1]][icdDate[[date1]] %in% as.Date(c("1900-01-01", "1999-01-01"))] <- NA }
vip <- read.table(paste0(indir,"common/ukb.vip.icd"), header=F, flush=T) %>%
	rename(trait=V1, code=V2)
dat <- subset(icd, select=f.eid) %>% rename(eid=f.eid)
for (i in 1:nrow(vip)) {
	dat1 <- icd[,-1]
	dat2 <- icdDate[,-1]
	exclude <- as.matrix(apply(dat1, 1:2, function(x){!grepl(vip$code[i],x)}))
	dat1[exclude] <- NA
	dat2[exclude] <- NA
	dat[[paste0("icd_",vip$trait[i])]] <- apply(dat1, 1, function(x){sum(!is.na(x))})
	dat[[paste0("icdDate_",vip$trait[i])]] <- as.Date(apply(dat2, 1, FUN=min, na.rm=T)) 
}
saveRDS(dat, file="D:/data/ukb/Rdata/ukb.icd.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FOD (first occurrence date)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(indir,"raw-50136/fod.r")) 
dates <- names(bd)[sapply(bd, is.Date)]; summary(bd[,dates])
for (date1 in dates) { bd[[date1]][bd[[date1]] %in% dates_bad] <- NA }
code12 <- read.table(paste0(indir,"common/ukb.fod.code12"), header=F) %>%
	rename(icd=V1, field=V2)
vip <- read.table(paste0(indir,"common/ukb.vip.fod"), header=F, flush=T) %>% 
	rename(trait=V1, type=V2, first=V3, last=V4) %>%
	merge(code12, by.x="first", by.y="icd", all.x=T) %>%
	merge(code12, by.x="last", by.y="icd", all.x=T) %>%
	mutate(
		first=ifelse(!is.na(field.x), field.x, first),
		last =ifelse(!is.na(field.y), field.y, last)
	)
	subset(vip, is.na(field.x) & type=="icd-10") #check
fields_all <- read.table(paste0(indir,"raw-50136/fields.ukb"), header=F)
for (i in 1:nrow(vip)) {
	fields <- paste0("f.", subset(fields_all, V1 >=vip[i,"first"] & V1 <=vip[i,"last"] & V1%%2==0)$V1, ".0.0")
	subdata <- subset(bd, select=fields)
	bd[[paste0("fod_",vip$trait[i])]] <- as.Date(apply(subdata, 1, FUN=min, na.rm=T)) 
}
fod <- subset(bd, select=grep("eid|fod_", names(bd), value=T)) %>% rename(eid=f.eid)
saveRDS(fod, file="D:/data/ukb/Rdata/ukb.fod.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge all together
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phe <- readRDS("D:/data/ukb/Rdata/ukb.phe.rds")
icd <- readRDS("D:/data/ukb/Rdata/ukb.icd.rds")
fod <- readRDS("D:/data/ukb/Rdata/ukb.fod.rds")
gen <- readRDS("D:/data/ukb/Rdata/ukb.gen.rds") 
hla <- readRDS("D:/data/ukb/Rdata/ukb.hla.rds") 
#hla <- hla %>% select(which(colMeans(.,na.rm=T) >=0.02))
prs <- read.table("D:/data/ukb/prs/all.prs.txt.gz", header=T, as.is=T)
prs <- prs %>% subset(select=!grepl("eid\\.", names(prs)))
dat0 = Reduce(function(x,y) merge(x,y,by="eid",all=T), list(phe, icd, fod, gen, prs)) 
saveRDS(dat0, file="D:/data/ukb/Rdata/all.Rdata")