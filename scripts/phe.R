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
# LE8
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat %>% 
mutate( 运动 🏃‍ MET minutes per week for moderate [22038] or vigorous [22039] activity
	physical_mins = met_mod + met_vig * 2, PA_pts = case_when(physical_mins >= 150 ~ 100, 
	physical_mins >= 120 & physical_mins < 150 ~ 90, physical_mins >= 90 & physical_mins < 120 ~ 80, physical_mins >= 60 & physical_mins < 90 ~ 60, physical_mins >= 30 & physical_mins < 60 ~ 40, physical_mins >= 1 & physical_mins < 30 ~ 20, physical_mins == 0 ~ 0)
	# 吸烟 🚬 smoke_quit_age [p2897], smoke_history [p1249], smoke_in_house [p1259]
    smoke_quit_age = ifelse(smoke_quit_age %in% c(-3, -1), NA, smoke_quit_age), smoke_history = ifelse(smoke_history == -3, NA, smoke_history), smoke_in_house = ifelse(smoke_in_house == -3, NA, smoke_in_house),
    smkquit_y = ifelse(smoke_quit_age > 0 & age >= smoke_quit_age, age - smoke_quit_age, NA), smoke_cat = case_when(smoking == 0 ~ 100, smoking == 1 & smkquit_y >= 5 ~ 75, smoke_history == 3 ~ 75, smoking == 1 & smkquit_y >= 1 & smkquit_y < 5 ~ 50, smoke_history == 2 ~ 50, smoking == 1 & smkquit_y < 1 ~ 25, smoking == 2 ~ 0),
    smoke_pts = ifelse((smoke_in_house == 1 | smoke_in_house == 2) & smoke_cat > 0, smoke_cat - 20, smoke_cat)
	# 睡眠🛏 sleep_dura [p1160]
    sleep_dura = ifelse(sleep_dura %in% c(-1, -3), NA, sleep_dura),
    sleep_pts = case_when(sleep_dura >= 7 & sleep_dura < 9 ~ 100, sleep_dura >= 9 & sleep_dura < 10 ~ 90, sleep_dura >= 6 & sleep_dura < 7 ~ 70, (sleep_dura >= 5 & sleep_dura < 6) | sleep_dura >= 10 ~ 40, sleep_dura >= 4 & sleep_dura < 5 ~ 20, sleep_dura >= 0 & sleep_dura < 4 ~ 0)
	# BMI 🎈
	bmi_pts = case_when(bmi_i > 0 & bmi_i < 25 ~ 100, bmi_i >= 25 & bmi_i < 30 ~ 70, bmi_i >= 30 & bmi_i < 35 ~ 30, bmi_i >= 35 & bmi_i < 40 ~ 15, bmi_i >= 40 ~ 0))
	# 血脂 
    nonhdl = pmax((chol_i - hdl_chol_i) * 38.67, 0),
    nonhdl_cat = case_when(nonhdl >= 0 & nonhdl < 130 ~ 100, nonhdl >= 130 & nonhdl < 160 ~ 60, nonhdl >= 160 & nonhdl < 190 ~ 40, nonhdl >= 190 & nonhdl < 220 ~ 20, nonhdl >= 220 ~ 0),
    nonhdl_pts = ifelse(lipid_drug == 1 & nonhdl_cat > 0, nonhdl_cat - 20, nonhdl_cat)
	# 血糖
    hba1c_n = hba1c_i * 0.0915 + 2.15,
    hba1c_pts = case_when(diabetes_his == 0 & hba1c_n > 0 & hba1c_n < 5.7 ~ 100, diabetes_his == 0 & hba1c_n >= 5.7 & hba1c_n < 6.5 ~ 60, hba1c_n > 0 & hba1c_n < 7 ~ 40, hba1c_n >= 7 & hba1c_n < 8 ~ 30, hba1c_n >= 8 & hba1c_n < 9 ~ 20, hba1c_n >= 9 & hba1c_n < 10 ~ 10, hba1c_n >= 10 ~ 0)
	# 血压
	BP_cat = case_when(sbp_i >= 160 | dbp_i >= 100 ~ 0, (sbp_i >= 140 & sbp_i < 160) | (dbp_i >= 90 & dbp_i < 100) ~ 25, (sbp_i >= 130 & sbp_i < 140) | (dbp_i >= 80 & dbp_i < 90) ~ 50, (sbp_i >= 120 & sbp_i < 130) & (dbp_i >= 0 & dbp_i < 80) ~ 75, sbp_i > 0 & sbp_i < 120 & dbp_i > 0 & dbp_i < 80 ~ 100),)
)
dat <- dat %>% rowwise() %>% # 🎇 汇总
	mutate( LE8score = mean(c(dashpts, PA_pts, smoke_pts, sleep_pts, bmi_pts, nonhdl_pts, hba1c_pts, BP_pts), na.rm = FALSE)) %>% 
	ungroup() %>%
	mutate(CVH_cat = case_when(LE8score < 50 & LE8score >= 0 ~ 0, LE8score < 80 & LE8score >= 50 ~ 1, LE8score >= 80 ~ 2 ))
LE8 <- dat %>% select(eid, LE8score, CVH_cat, dashpts, PA_pts, smoke_pts, sleep_pts, bmi_pts, nonhdl_pts, hba1c_pts, BP_pts)


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
	names(prot)[-1] <- gsub("^", "prot_", names(prot)[-1])
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
	names(hla)[-1]=paste0("hla_", names(hla)[-1])
	saveRDS(hla, paste0(indir,"/Rdata/ukb.hla.rds"))
pca <- read.table(paste0(indir,"/rap/pca.tab"), header=T)
	colnames(pca) <- c("eid", paste0("PC", 1:40))
	saveRDS(pca, paste0(indir,"/Rdata/ukb.pca.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GEN 和 PRS
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
prs0 <- read.table("D:/data/ukb/prs/all.prs.txt.gz", header=T, as.is=T)
	prs <- subset(prs0, select=grepl("^eid$|score_sum$|allele_cnt$", names(prs0)))
	saveRDS(prs, paste0(indir,"/Rdata/ukb.prs.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge all together
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (t in c("phe", "icd", "srd", "pca", "gen", "prs", "hla", "prot", "met")) {
	assign(t, readRDS(paste0(indir, '/Rdata/ukb.', t, '.rds')))
}
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
	sum(!is.na(dat0$icdDate_t2dm.2))	
	saveRDS(dat0, paste0(indir,"/Rdata/all.rds"))
dat0.plus <- Reduce(function(x,y) merge(x,y,by="eid",all=T), list(dat0, hla, prot, met))
saveRDS(dat0.plus, paste0(indir,"/Rdata/all.plus.rds"))


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
		bald1ab = ifelse(bald !=1, NA, ifelse(seq_len(nrow(dat)) %% 2 ==1, "a", "b")),
		abo.o.qt = ifelse(blood_group != "OO", NA, prot_abo),
		abo.a.qt = ifelse(!(blood_group %in% c("AO", "AA")), NA, prot_abo),
		abo.b.qt = ifelse(!(blood_group %in% c("BO", "BB")), NA, prot_abo),
		abo.ao.qt = ifelse(!(blood_group %in% c("OO", "AO", "AA")), NA, prot_abo)
		# lifespan里面的"绝大多数"是NA，第一个ifelse()如果不用 !is.na(lifespan)，那"绝大多数"就直接变成NA了，后面的">=80"就没用了
	) %>%
	dplyr::select(FID, IID, age, sex, PC1, PC2, PC3, PC4, bmi, height, leg, chunk, leg_ratio, leg_ratio_adj, chunk_ratio, chunk_ratio_adj, abo.o.qt,abo.a.qt,abo.b.qt,abo.ao.qt, bald12,bald13,bald14,bald134,bald1234, sexp25, sexp35, death60)
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