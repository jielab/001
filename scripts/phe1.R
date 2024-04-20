setwd("D:/")
pacman::p_load(data.table, dplyr, tidyverse, crosswalkr, lubridate)
indir="D:/data/ukb/phe/"
dates_bad=as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1999-01-01", "2037-07-07"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main PHE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(indir,"raw-670287/vip.r")); phe1 <- bd %>% rename(eid=f.eid) # use the latest version first
source(paste0(indir,"raw-50136/vip.r")); phe2 <- subset(bd, select=!names(bd) %in% names(phe1)) 
phe0 <- merge(phe1, phe2, by.x="eid", by.y="f.eid", all=TRUE) %>% 
	rename(f.25019.0.0=f.25019.2.0, f.25020.0.0=f.25020.2.0) 
phe0 <- phe0 %>% 
	dplyr::select(grep("eid|\\.0.0", names(phe0), value=T), f.48.2.0, f.49.2.0, f.21001.2.0) # 目前只保留 baseline 数据
	dates <- names(phe0)[sapply(phe0, is.Date)]; summary(phe0[,dates]) # check invalid dates
	for (date1 in dates) { phe0[[date1]][phe0[[date1]] %in% dates_bad] <- NA }
vip <- read.table(paste0(indir,"common/ukb.vip.dat"), header=FALSE, flush=TRUE)
	vip$V1[duplicated(vip$V1)]; vip$V2[duplicated(vip$V2)] #check duplication
	fid <- as.data.frame(names(phe0)[-1]); names(fid)="id1"
	fid[c('f0','f1','f2')] <- str_split_fixed(fid$id1, "\\.", 3)
	fid <- merge(fid, vip, by.x="f1", by.y="V1", all.x=TRUE) %>% mutate(id2 = ifelse(f2=="0.0", V2, paste0(V2,".",f2))) 
phe0 <- phe0 %>% crosswalkr::renamefrom(fid, id1, id2, drop_extra=F) %>% 
	filter(eid>0) %>% rename(chunk=height_sitting) 
phe <- phe0 %>% 
	mutate(
		edu = ifelse(edu==-7, 7, edu),
		across(grep("deprivation|date", names(phe0), invert=T, value=T), ~ifelse(.x<0, NA, .x)), # 该命令也不适用于date变量
		ethnic_cat = ifelse(ethnicity==2002,12, ifelse(ethnicity==2003,13, ifelse(ethnicity %in% c(2004,3004),6,  ifelse(grepl("^1",ethnicity),1, ifelse(grepl("^2",ethnicity),2, ifelse(grepl("^3",ethnicity),3, ifelse(grepl("^4",ethnicity),4, ethnicity))))))),
		ethnic_cat = factor(ethnic_cat, levels=c(1,2,3,4,5,6,12,13), labels=c("White","Mixed","Asian","Black","Chinese","Other","White-Black", "White-Asian")),	
		age_cat = cut(age, breaks=seq(35,75,5)),
		sex_cat = ifelse(sex==0, "female", "male"),
		bmi_cat = cut(bmi, breaks=c(10,18.5,25,30,100), labels=c("lean","healthy","overweight","obese")),
		bmi_cat = factor(bmi_cat, levels=c("healthy","lean","overweight","obese")),
		whr = waist / hip, whr.2.0 = waist.2.0 / hip.2.0, whr_diff = whr.2.0 - whr, bmi_diff = bmi.2.0 - bmi, 
		smoke_quite_year = age_smoke_quit - age_visit,
		smoke_score = ifelse(smoke_status==0,1, 0),
		smoke_status = factor(smoke_status, levels=0:2, labels=c("never","previous","current")), 
		cannabis_score = ifelse(cannabis==0,1, 0),
		alcohol_score = ifelse(alcohol_status==0,1, 0),
		alcohol_status = factor(alcohol_status, levels=0:2, labels=c("never","previous","current")),
		edu_score = ifelse(edu %in% 1:2,3, ifelse(edu %in% 3:6,2, 1)),  # 1-2: College or above; 3-6: High school or equivalent; -7: Less than high school
		emp_score = ifelse(emp %in% c(1,2,6,7), 2, 1), # 1 paid employment or self-employed; 2 retired; 6 doing unpaid or voluntary work; 7 full or part time students; 3 Looking after home and/or family; 4 Unable to work because of sickness or disability; 5 Unemployed 
		income_cat = ifelse(income==1,1, ifelse(income %in% 2:3,2, 3)),
		leg = height - chunk, leg_ratio = leg / chunk, chunk_ratio = chunk / leg,
		leg_ratio_adj = leg_ratio / height, chunk_ratio_adj = chunk_ratio / height,
		walk_brisk = ifelse(walk_pace==3, 1, 0),
		walk_freq = ifelse(walk_freq %in% 1:2, "low", ifelse(walk_freq %in% 3:4, "average", "high")),
		walk_time = ifelse(walk_time %in% 1:2, "short", ifelse(walk_time %in% 3:4, "average", "long")),
		walk_freq = factor(walk_freq, levels=c("low", "average", "high")),
		walk_time = factor(walk_time, levels=c("short", "average", "long")),
		birth_year = ifelse((birth_year<1936 | birth_year>1970), NA, birth_year), # 1934年1个，1971年2个，去掉
		birth_5year = cut(birth_year, breaks=seq(1935,1970,5)),
		age_2021 = 2021-birth_year, age_2021 = ifelse(is.na(date_lost) & is.na(date_death), age_2021, NA),
		lifespan = year(date_death) - birth_year
	)
saveRDS(phe, file="D:/data/ukb/Rdata/ukb.phe.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PC 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
umap <- read.table(paste0(indir,"raw-50136/pc-umap.txt"), header=F)
names(umap) <- c("eid", "umap1", "umap2")
source(paste0(indir,"raw-50136/pc.r"))
names(bd) <- gsub("f.22009.0.", "PC", names(bd))
pc <- bd %>% rename(eid=f.eid) %>% merge(umap, by="eid", all=TRUE)
saveRDS(pc, file="D:/data/ukb/Rdata/ukb.pc.rds")
write.table(pc, "PC.txt", append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ICD 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(indir,"raw-670287/icd10.r")); icd <- bd
source(paste0(indir,"raw-670287/icd10Date.r")); icdDate <- bd; bd <- NULL
dates <- names(icdDate)[sapply(icdDate, is.Date)]; summary(icdDate[,-1])
for (date1 in dates) { icdDate[[date1]][icdDate[[date1]] %in% as.Date(c("1900-01-01", "1999-01-01"))] <- NA }
vip <- read.table(paste0(indir,"common/ukb.vip.icd"), header=F, flush=T) %>%
	rename(trait=V1, code=V2)
dat <- subset(icd, select=f.eid) %>% rename(eid=f.eid)
for (i in 1:nrow(vip)) {
	dat1 <- icd[,-1]
	dat2 <- icdDate[,-1]
	exclude <- apply(dat1, 2, function(x){!grepl(vip$code[i],x)}) # 2 或者 1:2 得到 same dimension
	dat1[exclude] <- NA
	dat2[exclude] <- NA
	dat[[paste0("icd_",vip$trait[i])]] <- apply(dat1, 1, function(x){sum(!is.na(x))}) # apply 1
	dat[[paste0("icdDate_",vip$trait[i])]] <- as.Date(apply(dat2, 1, FUN=min, na.rm=T)) 
}
dat$icdDate_vte <- as.Date(apply(subset(dat, select=c(icdDate_dvt, icdDate_pe)), 1, FUN=min, na.rm=T))
hist(dat$icdDate_covid, breaks="weeks", freq=TRUE)
saveRDS(dat, file="D:/data/ukb/Rdata/ukb.icd.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FOD (first occurrence date)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(indir,"raw-670287/fod.r")) 
dates <- names(bd)[sapply(bd, is.Date)]; summary(bd[,dates])
for (date1 in dates) { bd[[date1]][bd[[date1]] %in% dates_bad] <- NA }
icdField <- read.table(paste0(indir,"common/ukb.icd-dataField.txt"), header=F) %>%
	rename(icd=V1, field=V2)
vip <- read.table(paste0(indir,"common/ukb.vip.fod"), header=F, flush=T) %>% 
	rename(trait=V1, code=V2) %>% 
	mutate(first=gsub("-.*", "", code), last=gsub(".*-", "", code)) %>%
	merge(icdField, by.x="first", by.y="icd", all.x=T) %>%
	merge(icdField, by.x="last", by.y="icd", all.x=T) %>%
	mutate(first=ifelse(!is.na(field.x), field.x, first), last =ifelse(!is.na(field.y), field.y, last))
fields_all <- read.table(paste0(indir,"raw-50136/fields.ukb"), header=F)
for (i in 1:nrow(vip)) {
	fields <- paste0("f.", subset(fields_all, V1 >=vip[i,"first"] & V1 <=vip[i,"last"] & V1%%2==0)$V1, ".0.0")
	subdata <- subset(bd, select=fields)
	bd[[paste0("fod_",vip$trait[i])]] <- as.Date(apply(subdata, 1, FUN=min, na.rm=T)) 
}
fod <- subset(bd, select=grep("eid|fod_", names(bd), value=T)) %>% rename(eid=f.eid)
saveRDS(fod, file="D:/data/ukb/Rdata/ukb.fod.rds")