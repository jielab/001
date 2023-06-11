setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, dplyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main GEN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
imp <- read.table(gzfile('D:/data/ukb/gen/imp/vip.raw.gz','r'), header=T, as.is=T) %>% rename(eid=IID) 
abo <- read.table('D:/data/ukb/phe/hes/covid19_misc.txt', header=T)
apoe <- read.table("D:/data/ukb/gen/hap/apoe.hap", header=T, as.is=T) %>% rename(eid=IID) 
sqc <- read.table("D:/data/ukb/phe/common/ukb_sqc_v2.txt", header=T, as.is=T) %>%
	rename(garray=genotyping.array, batch=Batch, in.british=in.white.British.ancestry.subset, in.pca=used.in.pca.calculation, aneuploidy=putative.sex.chromosome.aneuploidy)
sqc <- subset(sqc, select=grepl("eid|garray|batch|in\\.|aneuploidy|kinship|excess|PC", names(sqc)))
# 先在LINUX里面跑：cat ukb1941_rel_s488366.dat | awk 'NR>1 {print $2}' | sort | uniq | awk 'BEGIN{print "eid related"}{print $1,"Y"}' > ukb.related
rel <- read.table("D:/data/ukb/phe/common/ukb.related", header=T, as.is=T)
dat0 = Reduce(function(x,y) merge(x,y,by="eid",all=T), list(imp, abo, apoe, sqc, rel))
dat <- dat0 %>%
	mutate(
	abo=recode_factor(blood_group, "AA"="A", "BB"="B", "OO"="O", "AO"="A", "BO"="B"),
		o_type=ifelse(blood_group=="OO", "O", "non-O"),	
	apoe = factor(apoe, levels=c("e3e3", "e2e2", "e2e3", "e2e4", "e3e4", "e4e4")),
		e4 = ifelse(apoe=="e4e4", 2, ifelse(grepl("e4", apoe), 1, 0)),
		e4_f= factor(e4, levels=0:2, labels=c("zero", "one", "two")),
		e4_yes = ifelse(grepl("e4",apoe), 1, 0),
	sp1 = ifelse(sp1.Z.rs28929474_C==0, "ZZ", 
		ifelse (sp1.S.rs17580_T ==0, "SS", 
		ifelse( (sp1.Z.rs28929474_C==1 & sp1.S.rs17580_T==1), "SZ", 
		ifelse( (sp1.Z.rs28929474_C==1 & sp1.S.rs17580_T==2), "MZ", 
		ifelse( (sp1.Z.rs28929474_C==2 & sp1.S.rs17580_T==1), "MS", 
		ifelse( (sp1.Z.rs28929474_C==2 & sp1.S.rs17580_T==2), "MM", 
		NA))))))
	)
prop.table(table(dat$sp1)) 
saveRDS(dat, file="D:/data/ukb/Rdata/ukb.gen.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HLA 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source(paste0(indir,"raw-50136/hla.r")) 
	alleles <- read.table(paste0(indir,"raw-50136/ukb_hla_v2.txt"), header=T, as.is=T)
	hla <- bd %>% separate(f.22182.0.0, into=names(alleles), remove=T, convert=T, sep=","); str(hla)
names(hla)=gsub("hla_f.eid", "eid", paste0("hla_", names(hla)))
saveRDS(hla, file="D:/data/ukb/Rdata/ukb.hla.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ABO 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# abo.O2.rs41302905_T (G/T, T for O2), abo.A|B.rs8176746_T, abo.O1.rs8176719_TC (TC/T, for ins/del)
On = ifelse( (abo.O1.rs8176719_T==2 | abo.O2.rs41302905_C==0 | (abo.O1.rs8176719_T==1 & abo.O2.rs41302905_C<2)), 2, 
	ifelse( (abo.O1.rs8176719_T==1 | abo.O2.rs41302905_C==1), 1, 0)),
abo = ifelse(On==2, "OO", 
	ifelse(On==1 & abo.AB.rs8176746_G==2, "OA", 
	ifelse(On==1, "OB", 
	ifelse(abo.AB.rs8176746_G==0, "BB", 
	ifelse(abo.AB.rs8176746_G==1, "AB", "AA"))))),