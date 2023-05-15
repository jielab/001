setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, dplyr, tidyverse, crosswalkr, lubridate)
indir = "D:/data/ukb/phe/"
dates_bad=as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1999-01-01", "2037-07-07"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main PHE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PHElink <- read.table(paste0(indir,"common/ukb.vip.dat"), header=F, flush=T)
	PHElink$V1[duplicated(PHElink$V1)]; PHElink$V2[duplicated(PHElink$V2)] #check duplication
	PHElink$V2 <- paste0("f.", PHElink$V2, ".0.0")
source(paste0(indir,"raw-670287/vip.r")) # first use the latest version
	phe1 <- bd %>% subset(select=grep("f.eid|\\.0\\.0", names(bd))) %>% rename(eid=f.eid)
	names1 <- names(phe1)
source(paste0(indir,"raw-50136/vip.r")) # then add previous version
	phe2 <- bd %>% subset(select=!names(bd) %in% names1) 
	phe2 <- phe2 %>% subset(select=grep("f.eid|\\.0\\.0", names(phe2))) %>% rename(eid=f.eid)
phe0 <- merge(phe1, phe2, by="eid") %>% 
	crosswalkr::renamefrom(PHElink, V2, V1, drop_extra=F) %>% filter(eid>0)
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