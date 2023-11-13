setwd("D:/101")
pacman::p_load(data.table, lubridate, tidyverse, dplyr, ggplot2, CMplot, TwoSampleMR, MendelianRandomization, survival, survminer, mediation)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read and check data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata")
naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=T)), facet=sex)
grep("abo", names(dat0), value=T) # 找变量名字
table(dat0$abo, dat0$abo.O1.rs8176719); table(dat0$abo, dat0$abo.AB.rs8176746) # ABO血型
	cor(dat0$abo.microbe.rs545971, dat0$abo.microbe.rs8176645, use="complete.obs") # 0.85
	cor(dat0$abo.microbe.rs545971, dat0$abo.microbe.rs550057, use="complete.obs") # 0.83
	cor(dat0$abo.microbe.rs8176645, dat0$abo.microbe.rs550057, use="complete.obs") # 0.72
prop.table(table(dat0$abo, dat0$fut2.rs601338_A),1) # wild-type G allele encodes the "secretor" (Se) allele
hist(dat0$lct.MCM6.rs4988235_A, freq=F) # In Europe, T (或者说A) allele is common and lactase persistence, C/C（或者说G/G）is lactose intolerant.
	cor(dat0$lct.MCM6.rs4988235_A, dat0$lct.microbe.rs3940549_A, use="complete.obs") # 0.94 选这个就可以了
	cor(dat0$lct.MCM6.rs4988235_A, dat0$lct.microbe.rs182549_T, use="complete.obs") # 0.99
	cor(dat0$lct.microbe.rs182549_T, dat0$lct.microbe.rs3940549_A, use="complete.obs") # 0.94
table(dat0$sp1); table(dat0$sp1.M); table(dat0$sp1.S); table(dat0$sp1.Z)
hist(dat0$age_mn) #初潮年龄
table(dat0$sport.ACTN3.rs1815739); hist(dat0$sport.RBFOX1.rs7191721) 
table(dat0$icd_lungcancer); hist(dat0$icdDate_lungcancer, breaks="months") # 癌症
table(dat0$icd_chd); hist(dat0$icdDate_chd, breaks="months") # 心血管
table(dat0$icd_covid); hist(dat0$icdDate_covid, breaks="months") # COVID-19


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 万有引力之关联分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnicity_gen==1) # 欧洲白人
Xs <- grep(".rs", names(dat), value=T)
Ys <- grep("^icdDate|^fod_", names(dat), value=T)
for (Y in Ys) {
	print(paste("Y变量:", Y))
	dat1 <- dat %>% 
	mutate(
		outcome_date = dat[[Y]],
		outcome_yes = ifelse( is.na(dat[[Y]]), 0,1),
		follow_end_day = data.table::fifelse(!is.na(outcome_date), outcome_date, fifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25
	) %>% filter( follow_years >0 )
	print(table(dat1$outcome_yes))
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$outcome_yes)
	for (X in Xs) {
		print(paste("X变量:", X))
		dat1$X <- dat1[[X]]
		fit.suv <- survival::survfit(surv.obj ~ X + sex, data=dat1)
			#print(ggsurvplot(fit.suv, ylim=c(0.95,1), data=dat1)); dev.off() 
		fit.cox <- coxph(surv.obj ~ X + age+sex, data=dat1)
			#fit.cox <- coxph(surv.obj ~ X + X * fut2.rs601338_A + age+sex, data=dat1)
			print(summary(fit.cox))
			png(file=paste(X,Y,"frt.png",sep="."), w=1200, h=1600)
			print(ggforest(fit.cox, main=paste("X:", X, "| Y:", Y), fontsize=2.2, data=dat1)); dev.off()
	}
}
for (Y in grep("^bb_|^bc_", names(dat), value=T)) {
	print(noquote(Y))
	dat[[Y]] = inormal(dat[[Y]]) # normal transformation
	dat$Y = dat[[Y]] # assign temporary variable name "Y"
	for (X in Xs) {
		dat$X = dat[[X]]
		fit.lm <- lm(Y ~ X +age+sex, data=dat)
		print(noquote(paste("  --", X, signif(coef(summary(fit.lm))[2,4],3))))
	}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 北医合作研究案例：8 + happy 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# happy_general 20458; happy_health 20459; happy_life 20460; satisfy_bowel 21040
library(readxl)
wang <- read_excel("D:/data/ukb/phe/le8_data.xlsx") %>% rename(eid=n_eid, sex.y=sex)
dat <- dat0 %>% filter(ethnicity_gen==1) %>% 
	merge(wang, by="eid")
Y="fod_chd"
dat1 <- dat %>% 
	mutate( 
		outcome_date = dat[[Y]],
		outcome_yes = ifelse( is.na(outcome_date), 0,1),
		follow_end_day = ifelse(!is.na(outcome_date), outcome_date, ifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25
	) %>% filter( follow_years >0 )
	print(table(dat1$outcome_yes))
dat1$happy = ifelse(dat1$happy_health>=4,0, ifelse(dat1$happy_general>=1,1,NA))
#dat1$happy = ifelse(dat1$satisfy_bowel<=4,1, ifelse(dat1$satisfy_bowel>=6,0,NA))
surv.obj <- Surv(time=dat1$follow_years, event=dat1$outcome_yes)
fit.cox <- coxph(surv.obj ~ happy + cad.score_sum + dashpts+PA_pts+smoke_pts+sleep_pts+bmi_pts+nonhdl_pts+hba1c_pts+BP_pts +age+sex+PC1+PC2, data=dat1); summary(fit.cox)
}

