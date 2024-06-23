setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, tidyverse)
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# main GEN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
imp <- read.table('D:/data/ukb/gen/imp/vip.raw.gz', header=T, as.is=T) %>% rename(eid=IID) 
abo <- read.table('D:/data/ukb/phe/hes/covid19_misc.txt', header=T)
apoe <- read.table("D:/data/ukb/gen/hap/apoe.hap", header=T, as.is=T) %>% rename(eid=IID) 
dat0 = Reduce(function(x,y) merge(x,y,by="eid",all=T), list(imp, abo, apoe))
dat <- dat0 %>%
	mutate(
	abo=recode_factor(blood_group, "AA"="A", "BB"="B", "OO"="O", "AO"="A", "BO"="B"),
	abo.a=ifelse(abo=="A", 1, 0), abo.o=ifelse(abo=="O", 1, 0), abo.a_b=ifelse(abo=="O", NA, ifelse(abo=="A", 1, 0)),
	abo.se=ifelse(fut2.rs601338_A==2, 0, 1), # [PMID: 30345375]
	vte.f2 = hardcall(vte.F2.rs1799963_G), vte.f5= hardcall(vte.F5.rs6025_C),
#	vte.ff = ifelse( (vte.f2==0 | vte.f5==0), "homozygote", ifelse((vte.f2==1 & vte.f5==1), "compound", ifelse((vte.f2==2 & vte.f5==2), "wild-type", "heterozygote"))), 
	vte.ff = ifelse( (vte.f2!=2 | vte.f5!=2), 1,0),
	apoe = factor(apoe, levels=c("e3e3", "e2e2", "e2e3", "e2e4", "e3e4", "e4e4")),
		apoe.e4 = ifelse(grepl("e4", apoe), 1, 0),
	sp1 = ifelse(sp1.Z.rs28929474_C==0, "ZZ", 
		ifelse (sp1.S.rs17580_T ==0, "SS", 
		ifelse( (sp1.Z.rs28929474_C==1 & sp1.S.rs17580_T==1), "SZ", 
		ifelse( (sp1.Z.rs28929474_C==1 & sp1.S.rs17580_T==2), "MZ", 
		ifelse( (sp1.Z.rs28929474_C==2 & sp1.S.rs17580_T==1), "MS", 
		ifelse( (sp1.Z.rs28929474_C==2 & sp1.S.rs17580_T==2), "MM", 
		NA)))))),
	sp1.s.012 = ifelse(grepl("SS", sp1), 2, ifelse(grepl("S", sp1), 1, 0)),
	sp1.z.012 = ifelse(grepl("ZZ", sp1), 2, ifelse(grepl("Z", sp1), 1, 0)),
	sp1.s = ifelse(sp1.s.012==0, 0, 1), 
	sp1.z = ifelse(sp1.z.012==0, 0, 1)
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