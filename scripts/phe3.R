setwd("D:/")
pacman::p_load(data.table, dplyr, tidyverse, lubridate, naniar)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# merge all together 并生成GWAS需要的.pheno文件
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phe <- readRDS("D:/data/ukb/Rdata/ukb.phe.rds")
pc  <- readRDS("D:/data/ukb/Rdata/ukb.pc.rds")
icd <- readRDS("D:/data/ukb/Rdata/ukb.icd.rds")
srd <- readRDS("D:/data/ukb/Rdata/ukb.srd.rds")
fod <- readRDS("D:/data/ukb/Rdata/ukb.fod.rds")
gen <- readRDS("D:/data/ukb/Rdata/ukb.gen.rds") 
hla <- readRDS("D:/data/ukb/Rdata/ukb.hla.rds") #hla <- hla %>% select(which(colMeans(.,na.rm=T) >=0.02))
prs0 <- read.table("D:/data/ukb/prs/all.prs.txt.gz", header=T, as.is=T); prs <- subset(prs0, select=grepl("^eid$|score_sum$|allele_cnt$", names(prs0)))
dat0 = Reduce(function(x,y) merge(x,y,by="eid",all=T), list(phe, pc, icd, srd, fod, gen, prs)) %>% filter(eid>0)
saveRDS(dat0, file="D:/data/ukb/Rdata/all.Rdata")


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
dat <- dat0 %>% filter(ethnic_cat=="White") %>% rename(IID=eid) %>% 
	mutate(
		FID = IID,
		death60a = ifelse(!is.na(lifespan) & lifespan <=60, 1, ifelse(age_2021>=80, 0, NA)), 
		death60b = ifelse(!is.na(lifespan) & lifespan <=60, 1, ifelse(age_2021>=75, 0, NA)), 
		death60c = ifelse(!is.na(lifespan) & lifespan <=65, 1, ifelse(age_2021>=80, 0, NA))
		# lifespan里面的"绝大多数"是NA，第一个ifelse()如果不用 !is.na(lifespan)，那"绝大多数"就直接变成NA了，后面的">=80"就没用了
	) %>%
	dplyr::select(FID, IID, age, sex, PC1, PC2, PC3, PC4, bmi, height, leg, chunk, leg_ratio, leg_ratio_adj, chunk_ratio, chunk_ratio_adj, walk_pace, walk_brisk, death60a, death60b, death60c)
	write.table(dat, "ukb.pheno", append=FALSE, quote=FALSE, row.names=FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# find trio
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ped <- subset(dat0, ethnicity_gen==1, select=c("eid", "age", "sex", "birth_year"))
kin <- read.table("D:/data/ukb/phe/common/ukb.kin0", header=T, as.is=T) %>% 
	subset(ID1>0 & InfType=="PO", select=c("ID1", "ID2"))
dat <- merge(ped, kin, by.x="eid", by.y="ID1")
dat <- merge(ped, dat, by.x="eid", by.y="ID2") %>% 
	subset(abs(age.x - age.y)> 18) %>% rename(eid.x=eid) %>%
	mutate(
		child = ifelse(age.x > age.y, eid.y, eid.x),
		father = ifelse(eid.x==child & sex.y==1, eid.y, ifelse(eid.y==child & sex.x==1, eid.x, NA)),
		mother = ifelse(eid.x==child & sex.y==0, eid.y, ifelse(eid.y==child & sex.x==0, eid.x, NA))
	)
father <- subset(dat, !is.na(father), select=c("child", "father"))
mother <- subset(dat, !is.na(mother), select=c("child", "mother"))
trio <- merge(father, mother)
