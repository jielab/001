setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, lubridate, tidyverse, dplyr, ggplot2, CMplot, TwoSampleMR, MendelianRandomization, survival, survminer, mediation)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read and check data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata")
	naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=T)), facet=sex)
	prop.table(table(dat0$abo, dat0$fut2.rs601338_A),1)
	table(dat0$sp1)
	hist(dat0$lac.rs4988235_A)
	hist(dat0$age_mn)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 万有引力之生存分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnicity_gen==1) 
Xs <- c("o_type", "fut2.rs601338_A", "lac.rs4988235_A", "sp1.Z")
Ys <- c("icdDate_copd", "icdDate_covid", "fod_t1dm", "fod_t2dm", "fod_dm")
for (Y in Ys) {
	print(paste("Y变量:", varY))
	dat1 <- dat %>% 
	mutate(
		outcome_date = dat[[varY]],
		outcome_yes = ifelse( is.na(dat[[varY]]), 0,1),
		follow_end_day = data.table::fifelse(!is.na(outcome_date), outcome_date, fifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25
	) %>% filter( follow_years >0 )
	print(table(dat1$outcome_yes))
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$outcome_yes)
	for (varX in varXs) {
		print(paste("X变量:", varX))
		dat1$varX <- dat1[[varX]]
		fit.suv <- survival::survfit(surv.obj ~ varX + sex, data=dat1)
			#print(ggsurvplot(fit.suv, ylim=c(0.95,1), data=dat1)); dev.off() 
		fit.cox <- coxph(surv.obj ~ varX + age+sex, data=dat1)
			#fit.cox <- coxph(surv.obj ~ varX + varX * fut2.rs601338_A + age+sex, data=dat1)
			print(summary(fit.cox))
			png(file=paste(varX,varY,"frt.png",sep="."), w=600, h=800)
			print(ggforest(fit.cox, data=dat1)); dev.off()
	}
}
for (varY in grep("^bb_|^bc_", names(dat), value=T)) {
	print(noquote(varY))
	dat[[varY]] = inormal(dat[[varY]]) # normal transformation
	dat$varY = dat[[varY]] # assign temporary variable name "varY"
	for (varX in varXs) {
		dat$varX = dat[[varX]]
		fit.lm <- lm(varY ~ varX +age+sex, data=dat)
		print(noquote(paste("  --", varX, signif(coef(summary(fit.lm))[2,4],3))))
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
varY="fod_chd"
dat1 <- dat %>% 
	mutate( 
		outcome_date = dat[[varY]],
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

