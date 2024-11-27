setwd("D:/")
pacman::p_load(data.table, tidyverse, stringr, crosswalkr, lubridate)
indir="D:/data/ukb/phe"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main PHE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vip.dat <- read.table(paste0(indir,"/common/ukb.vip.dat"), header=FALSE, flush=TRUE)
	vip.dat$V1[duplicated(vip.dat$V1)]; vip.dat$V2[duplicated(vip.dat$V2)] #check duplication
phe0 <- read.delim(paste0(indir,"/rap/vip.tab"), sep="\t", header=T) # 不能用 read.table 🈲
	names(phe0) <- gsub("_i0_a0$|_i0$", "", names(phe0))
	phe0 <- phe0[, !grepl("_", names(phe0))]; names(phe0) 
	setdiff(vip.dat$V1, names(phe0)) # 查找不overlap的column
	phe0 <- phe0 %>% crosswalkr::renamefrom(cw_file=vip.dat, raw=V1, clean=V2, drop_extra=F)
phe <- phe0 %>% 
	mutate(
		across(grep("date", names(phe0), value=T), ~as.Date(.x)),
		across(grep("deprivation|date", names(phe0), invert=T, value=T), ~ifelse(.x<0, NA, .x)), # 该命令也不适用于date变量
		ethnic_cat = ifelse(ethnicity==2002,12, ifelse(ethnicity==2003,13, ifelse(ethnicity %in% c(2004,3004),6,  ifelse(grepl("^1",ethnicity),1, ifelse(grepl("^2",ethnicity),2, ifelse(grepl("^3",ethnicity),3, ifelse(grepl("^4",ethnicity),4, ethnicity))))))),
		ethnic_cat = factor(ethnic_cat, levels=c(1,2,3,4,5,6,12,13), labels=c("White","Mixed","Asian","Black","Chinese","Other","White-Black", "White-Asian")),	
		age_cat = cut(age, breaks=seq(35,75,5)),
		sex_cat = ifelse(sex==0, "female", "male"),
		bmi_cat = cut(bmi, breaks=c(10,18.5,25,30,100), labels=c("lean","healthy","overweight","obese")),
		bmi_cat = factor(bmi_cat, levels=c("healthy","lean","overweight","obese")),
		smoke_quite_year = age_smoke_quit - age_visit,
		smoke_score = ifelse(smoke_status==0,1, 0),
		smoke_status = factor(smoke_status, levels=0:2, labels=c("never","previous","current")), 
		cannabis_score = ifelse(cannabis==0,1, 0),
		alcohol_score = ifelse(alcohol_status==0,1, 0),
		alcohol_status = factor(alcohol_status, levels=0:2, labels=c("never","previous","current")),
		sexp_cat =ifelse(sex_partner > 10, "high", ifelse(sex_partner>3, "middle", "low")),
		sexp_cat = factor(sexp_cat, levels=c("low", "middle", "high")),
		leg = height - chunk, leg_ratio = leg / chunk, chunk_ratio = chunk / leg,
		leg_ratio_adj = leg_ratio / height, chunk_ratio_adj = chunk_ratio / height,
		birth_year = ifelse(birth_year>1970, NA, birth_year), # 1936-1971, 去掉1个1971年的
		birth_5year = cut(birth_year, breaks=seq(1935,1970,5)),
		age_facial = factor(age_facial, labels=c("young", "old", "normal")),
		lifespan = year(date_death) - birth_year,
		sexp25 =ifelse(sex_partner > 5, 1, ifelse(sex_partner>2, NA,0)),
		sexp35 =ifelse(sex_partner > 5, 1, ifelse(sex_partner>3, NA,0)),
		age_2021 = 2021-birth_year, age_2021 = ifelse(is.na(date_lost) & is.na(date_death), age_2021, NA),
		death60 = ifelse(!is.na(lifespan) & lifespan <=60, 1, ifelse(age_2021>=70, 0, NA))

	)
saveRDS(phe, paste0(indir,"/Rdata/ukb.phe.rds"))

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ICD10 (41270 和 41280)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
lines <- readLines(paste0(indir,"/rap/icd10.tab")) # 这个文件
	split_lines <- strsplit(lines, split = "\t") 
	max_cols <- max(sapply(split_lines, length))
	icd <- do.call(rbind, lapply(split_lines, function(row) {
		length(row) <- max_cols
		return(row)
	}))
	icd <- as.data.frame(icd, stringsAsFactors=FALSE)
	colnames(icd) <- c("eid", paste0("I", seq_len(max_cols-1)))
icdDate <- read.table(paste0(indir,"/rap/icd10Date.tab"), header=T)
	colnames(icdDate) <- c("eid", paste0("D", seq_len(max_cols-1)))
	icdDate[-1] <- lapply(icdDate[-1], as.Date)
	min(icdDate$D1, na.rm=TRUE)
vip.icd <- read.table(paste0(indir,"/common/ukb.vip.icd"), header=F, flush=T) %>% rename(trait=V1, code=V2)
dat <- subset(icd, select=eid)
	for (i in 1:nrow(vip.icd)) {
		datCode <- icd[,-1]
		datDate <- icdDate[,-1]
		exclude <- apply(datCode, 2, function(x){!grepl(vip.icd$code[i],x)}) # 2 或者 1:2 得到 same dimension
		datCode[exclude] <- NA
		datDate[exclude] <- NA
		dat[[paste0("icdCt_",vip.icd$trait[i])]] <- apply(datCode, 1, function(x){sum(!is.na(x))}) # apply 1
		dat[[paste0("icdDate_",vip.icd$trait[i])]] <- as.Date(apply(datDate, 1, FUN=min, na.rm=T)) 
	}
hist(dat$icdDate_cad, breaks="weeks", freq=TRUE)
saveRDS(dat, paste0(indir,"/Rdata/ukb.icd.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SRD (self-report disease) (20002 和 87) 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srd <- read.table(paste0(indir,"/rap/srd.tab"), header=T)
srdTime <- read.table(paste0(indir,"/rap/srdTime.tab"), header=T)
vip.srd <- read.table(paste0(indir,"/common/ukb.vip.srd"), header=F, flush=T) %>% rename(trait=V1, code=V2)
dat <- subset(srd, select=eid)
for (i in 1:nrow(vip.srd)) {
	datCode <- srd[,-1]
	datAge <- sapply(srdTime[,-1], as.numeric)
	datYear <- sapply(srdTime[,-1], as.numeric)
	exclude <- apply(datCode, 2, function(x){!grepl(vip.srd$code[i],x)}) 
	datCode[exclude] <- NA
	datAge[exclude] <- NA; datAge[datAge >100] <- NA
	datYear[exclude] <- NA; datYear[datYear <1930] <- NA
	dat[[paste0("srdCt_",vip.srd$trait[i])]] <- apply(datCode, 1, function(x){sum(!is.na(x))})
	dat[[paste0("srdAge_",vip.srd$trait[i])]] <- apply(datAge, 1, FUN=min, na.rm=TRUE)
	dat[[paste0("srdYear_",vip.srd$trait[i])]] <- apply(datYear, 1, FUN=min, na.rm=TRUE)
}
dat[sapply(dat, is.infinite)] <- NA #!!
saveRDS(dat, paste0(indir,"/Rdata/ukb.srd.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# prot 和 met 数据🏮🏮，暂不合并 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prot <- read.table(paste0(indir,"/rap/raw/prot.tab"), sep="\t", header=T)
	saveRDS(prot, paste0(indir,"/Rdata/ukb.prot.rds"))
vip.met <- read.table(paste0(indir,"/common/ukb.vip.met"), header=FALSE, flush=TRUE)
	vip.met$V1[duplicated(vip.met$V1)]; vip.met$V2[duplicated(vip.met$V2)]
	met <- read.delim(paste0(indir,"/rap/met.tab"), sep="\t", header=T)
	names(met) <- gsub("_i0$", "", names(met))	
	setdiff(vip.met$V1, names(met)) # 查找不overlap的column
	met <- met %>% crosswalkr::renamefrom(cw_file=vip.met, raw=V1, clean=V2, drop_extra=F)
	saveRDS(met, paste0(indir,"/Rdata/ukb.met.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HLA 和 PCA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hla <- read.table(paste0(indir,"/rap/hla.tab"), sep="\t", header=T)
	alleles <- read.table(paste0(indir,"/50136/ukb_hla_v2.txt"), header=T, as.is=T)
	hla <- hla %>% separate(p22182, into=names(alleles), remove=T, convert=T, sep=","); str(hla)
	names(hla)=paste0("hla_", names(hla))
	saveRDS(phe, paste0(indir,"/Rdata/ukb.hla.rds"))
pca <- read.table(paste0(indir,"/rap/pca.tab"), header=T)
	colnames(pca) <- c("eid", paste0("PC", 1:40))
	saveRDS(pca, paste0(indir,"/Rdata/ukb.pca.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main GEN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# abo.O2.rs41302905_T (G/T, T for O2), abo.A|B.rs8176746_T, abo.O1.rs8176719_TC (TC/T, for ins/del)
# On = ifelse( (abo.O1.rs8176719_T==2 | abo.O2.rs41302905_C==0 | (abo.O1.rs8176719_T==1 & abo.O2.rs41302905_C<2)), 2, ifelse( (abo.O1.rs8176719_T==1 | abo.O2.rs41302905_C==1), 1, 0)),
# abo = ifelse(On==2, "OO", ifelse(On==1 & abo.AB.rs8176746_G==2, "OA", ifelse(On==1, "OB", ifelse(abo.AB.rs8176746_G==0, "BB", ifelse(abo.AB.rs8176746_G==1, "AB", "AA"))))),
imp <- read.table('D:/data/ukb/gen/imp/vip.raw.gz', header=T, as.is=T) %>% rename(eid=IID) 
abo <- read.table('D:/data/ukb/phe/hes/covid19_misc.txt', header=T)
apoe <- read.table("D:/data/ukb/gen/hap/apoe.hap", header=T, as.is=T) %>% rename(eid=IID) 
dat0 = Reduce(function(x,y) merge(x,y,by="eid",all=T), list(imp, abo, apoe))
dat <- dat0 %>%
	mutate(
	abo=recode_factor(blood_group, "AA"="A", "BB"="B", "OO"="O", "AO"="A", "BO"="B"),
		abo.a.add = ifelse(blood_group =="OO", 0, ifelse(blood_group =="AO", 1, ifelse(blood_group =="AA", 2, NA))),
		abo.a.dom = ifelse(abo =="O", 0, ifelse(abo =="A", 1, NA)),
		se=ifelse(fut2.rs601338_A==2, 0, 1), # [PMID: 30345375]
	apoe = factor(apoe, levels=c("e3e3", "e2e2", "e2e3", "e2e4", "e3e4", "e4e4")),
		apoe.e4 = ifelse(grepl("e4", apoe), 1, 0),
	sp1 = ifelse(sp1.Z.rs28929474_C==0, "ZZ", 
		ifelse (sp1.S.rs17580_T ==0, "SS", 
		ifelse( (sp1.Z.rs28929474_C==1 & sp1.S.rs17580_T==1), "SZ", 
		ifelse( (sp1.Z.rs28929474_C==1 & sp1.S.rs17580_T==2), "MZ", 
		ifelse( (sp1.Z.rs28929474_C==2 & sp1.S.rs17580_T==1), "MS", 
		ifelse( (sp1.Z.rs28929474_C==2 & sp1.S.rs17580_T==2), "MM", 
		NA)))))),
	sp1.s.add = ifelse(grepl("SS", sp1), 2, ifelse(grepl("S", sp1), 1, 0)),
	sp1.z.add = ifelse(grepl("ZZ", sp1), 2, ifelse(grepl("Z", sp1), 1, 0)),
	sp1.s = ifelse(sp1.s.add==0, 0, 1), 
	sp1.z = ifelse(sp1.z.add==0, 0, 1)
	)
table(dat$sp1)
saveRDS(dat, paste0(indir,"/Rdata/ukb.gen.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge all together
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phe <- readRDS(paste0(indir,"/Rdata/ukb.phe.rds"))
icd <- readRDS(paste0(indir,"/Rdata/ukb.icd.rds"))
srd <- readRDS(paste0(indir,"/Rdata/ukb.srd.rds"))
pca <- readRDS(paste0(indir,"/Rdata/ukb.pca.rds"))
gen <- readRDS(paste0(indir,"/Rdata/ukb.gen.rds")) 
prs0 <- read.table("D:/data/ukb/prs/all.prs.txt.gz", header=T, as.is=T)
	prs <- subset(prs0, select=grepl("^eid$|score_sum$|allele_cnt$", names(prs0)))
dat0 <- Reduce(function(x,y) merge(x,y,by="eid",all=T), list(phe, icd, srd, pca, gen, prs))
dat0$icdDate_vte <- as.Date(apply(subset(dat0, select=c(icdDate_dvt, icdDate_pe)), 1, FUN=min, na.rm=T))
vip.srd <- read.table(paste0(indir,"/common/ukb.vip.srd"), header=F, flush=T) %>% rename(trait=V1, code=V2)
for (t in vip.srd$trait) {
	dat0[[paste0("icdDate_",t,".2")]] <- as.Date(
		ifelse(!is.na(dat0[[paste0("icdDate_",t)]]), as.character(dat0[[paste0("icdDate_",t)]]), 
		ifelse(!is.na(dat0[[paste0("srdYear_",t)]]), (paste0(dat0[[paste0("srdYear_",t)]],"-07-01")), 
		ifelse(!is.na(dat0[[paste0("srdAge_",t)]]) & dat0[[paste0("srdAge_",t)]] >0, paste0(dat0$birth_year+dat0[[paste0("srdAge_",t)]],"-07-01"), NA)))
	)
}
sum(!is.na(dat0$icdDate_t2dm2))	
saveRDS(dat0, paste0(indir,"/Rdata/all.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 做一点 sanity check
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(phe$ethnicity, phe$ethnic_cat, useNA="always")
	naniar::gg_miss_var(subset(phe, select=grep("sex|bb_", names(phe), value=TRUE)), facet=sex)
	hist(phe$bmi_diff, nclass=100)
group_by(dat0, abo, se) %>% summarise(count=n(), mean=mean(bb_TES, na.rm=TRUE))
	aggregate(bb_TES ~ abo*se, dat0, FUN=function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2)} )
	bp <- boxplot(bb_TES ~ abo*se, dat0, las=2, col=rainbow(6), font=2); bp$stats
	dat0 %>% drop_na(abo, se) %>% ggplot(aes(x=abo, y=bb_TES, fill=se)) + geom_boxplot() + stat_summary(fun.y=mean, color="darkred", position=position_dodge(0.75), geom="point", shape=18, size=3)
dat <- dat0 %>% select(grep("bb_ALB|bb_APOB|bb_ALP|bb_CYS|bb_HDL|bb_LDL", names(dat0), value=TRUE)) %>% na.omit() %>% dplyr::sample_n(10000)
	car::scatterplotMatrix(dat, spread=FALSE, smoother.args=list(lty=0.1))
	psych::pairs.panels(dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 生成GWAS所需的表型数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnic_cat=="White") %>% rename(IID=eid) 
	dat <- dat %>% 
	mutate(
		FID = IID,
		bald12 = ifelse(bald==1,0, ifelse(bald==2,1,NA)),
		bald13 = ifelse(bald==1,0, ifelse(bald==3,1,NA)),
		bald14 = ifelse(bald==1,0, ifelse(bald==4,1,NA)),
		bald134 = ifelse(bald==1,0, ifelse(bald==3 | bald==4,1,NA)),
		bald1234 = ifelse(bald==1,0, 1),
		bald1ab = ifelse(bald!=1, NA, ifelse(seq_len(nrow(dat)) %% 2 ==1, "a", "b")),
		# lifespan里面的"绝大多数"是NA，第一个ifelse()如果不用 !is.na(lifespan)，那"绝大多数"就直接变成NA了，后面的">=80"就没用了
	) %>%
	dplyr::select(FID, IID, age, sex, PC1, PC2, PC3, PC4, bmi, height, leg, chunk, leg_ratio, leg_ratio_adj, chunk_ratio, chunk_ratio_adj, walk_pace, walk_brisk, a01,b01,a12,b12,a012,b012, bald12,bald13,bald14,bald1a2,bald1b2,bald1a3,bald1b3, bald1a4,bald1b4, bald134, bald1234, sexp25, sexp35, death60)
	write.table(dat, "ukb.pheno", append=FALSE, quote=FALSE, row.names=FALSE)
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🈲 main PHE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(indir,"/raw-670287/vip.r")); phe1 <- bd %>% rename(eid=f.eid) # use the latest version first
source(paste0(indir,"/raw-50136/vip.r")); phe2 <- subset(bd, select=!names(bd) %in% names(phe1)) 
phe0 <- merge(phe1, phe2, by.x="eid", by.y="f.eid", all=TRUE) 
	dates <- names(phe0)[sapply(phe0, is.Date)]; summary(phe0[,dates]) # check invalid dates
	dates_bad=as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1999-01-01", "2037-07-07"))
	for (date1 in dates) { phe0[[date1]][phe0[[date1]] %in% dates_bad] <- NA }
vip.dat <- read.table(paste0(indir,"/common/ukb.vip.dat"), header=FALSE, flush=TRUE)
	vip.dat$V1[duplicated(vip.dat$V1)]; vip.dat$V2[duplicated(vip.dat$V2)] #check duplication
	fid <- as.data.frame(names(phe0)[-1]); names(fid)="id1"
	fid[c('f0','f1','f2')] <- stringr::str_split_fixed(fid$id1, "\\.", 3)
	fid <- merge(fid, vip.dat, by.x="f1", by.y="V1", all.x=TRUE) %>% mutate(id2 = ifelse(f2=="0.0", V2, paste0(V2,".",f2))) 
phe0 <- phe0 %>% crosswalkr::renamefrom(fid, id1, id2, drop_extra=F) %>% 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🈲 FOD (first occurrence date)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(indir,"/raw-670287/fod.r")) 
dates <- names(bd)[sapply(bd, is.Date)]; summary(bd[,dates])
for (date1 in dates) { bd[[date1]][bd[[date1]] %in% dates_bad] <- NA }
icdField <- read.table(paste0(indir,"/common/ukb.icd-dataField.txt"), header=F) %>%
	rename(icd=V1, field=V2)
vip.fod <- read.table(paste0(indir,"/common/ukb.vip.fod"), header=F, flush=T) %>% 
	rename(trait=V1, code=V2) %>% 
	mutate(first=gsub("-.*", "", code), last=gsub(".*-", "", code)) %>%
	merge(icdField, by.x="first", by.y="icd", all.x=T) %>%
	merge(icdField, by.x="last", by.y="icd", all.x=T) %>%
	mutate(first=ifelse(!is.na(field.x), field.x, first), last =ifelse(!is.na(field.y), field.y, last))
fields_all <- read.table(paste0(indir,"/raw-50136/fields.ukb"), header=F)
for (i in 1:nrow(vip.fod)) {
	fields <- paste0("f.", subset(fields_all, V1 >=vip.fod[i,"first"] & V1 <=vip.fod[i,"last"] & V1%%2==0)$V1, ".0.0")
	subdata <- subset(bd, select=fields)
	bd[[paste0("fod_",vip.fod$trait[i])]] <- as.Date(apply(subdata, 1, FUN=min, na.rm=T)) 
}
fod <- subset(bd, select=grep("eid|fod_", names(bd), value=T)) %>% rename(eid=f.eid)
saveRDS(fod, file="D:/data/ukb/Rdata/ukb.fod.rds")