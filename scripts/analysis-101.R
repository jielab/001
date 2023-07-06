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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 万有引力之生存分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnicity_gen==1) %>% 
	mutate (
		sp1.M = ifelse(sp1=="MM", 2, ifelse(grepl("M", sp1), 1, 0)),
		sp1.S = ifelse(grepl("SS", sp1), 2,  ifelse(grepl("S", sp1), 1, 0)),
		sp1.Z = ifelse(grepl("ZZ", sp1), 2,  ifelse(grepl("Z", sp1), 1, 0))
	)
varXs <- c("o_type", "fut2.rs601338_A", "lac.rs4988235_A", "sp1.Z")
for (varY in c("icdDate_copd", "icdDate_covid", "fod_t1dm", "fod_t2dm", "fod_dm")) {
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
# Mediation 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnicity_gen==1) %>% dplyr::select(grep("age|sex|bmi|covid|^bb_|^bc|age_mn",names(dat0),value=T)) %>% 
	mutate (
	age_mn_3p = ifelse(age_mn %in% 12:14, "normal", ifelse(age_mn>14, "late", "early")),
	age_mn_3p = factor(age_mn_3p, levels=c("normal", "early", "late")),
	outcome = ifelse(!is.na(icdDate_covid), 1, ifelse(covid_inf==1,0, NA))
	) %>% rename(CYS=bb_CYS)
varX="age_mn_3p"; varY="outcome"
form <- formula(paste(varY, "~", varX, "+age+sex+bmi+CYS"))
form <- formula(paste(varY, "~", varX, "+age+sex+bmi+CYS+", paste(grep("^bb_", names(dat), value=T), sep="", collapse="+")))
summary(glm(form, data=dat))
for (varM in grep("CYS|^bb_|^bc_", names(dat), value=T)) {
	print(paste("M变量:", varM))
	dat1 <- subset(dat, select=c(varX, varY, varM, "age", "sex", "bmi")) %>% na.omit()
	dat1[[varM]] <- inormal(dat1[[varM]]) # normal transformation
	names(dat1) <- c("varX", "varY", "varM", "age","sex", "bmi")
	dat1$varX <- as.factor(dat1$varX)
	fit.med = lm(varM ~ varX + age+sex, data=dat1)
	fit.out = glm(varY ~ varX + varM +age+sex, data=dat1, family=binomial("probit"))
	med.out = mediation::mediate(fit.med, fit.out, treat="varX", mediator="varM", sims=100, boot=T)
	print(summary(med.out))
}
fit.med = lm(varM ~ varX, data=dat1)
fit.out = glm(varY ~ varX + varM, data=dat1, family=binomial("probit"))
med.out = mediation::mediate(fit.med, fit.out, treat="varX", mediator="varM", sims=100, boot=T)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# meta-analysis
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(meta)
data("Fleiss1993cont"); head(Fleiss1993cont)
res <- metacont(n.psyc, mean.psyc, sd.psyc, n.cont, mean.cont, sd.cont,
	comb.fixed=T, comb.random=T, studlab=study, sm="SMD", data=Fleiss1993cont) 
res
forest(res, leftcols = c('studlab'))
funnel(res)
metabias(res, method.bias = 'linreg', k.min = 5, plotit = T)# Egger


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Circular Manhattan Plot 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- read.table('D:/data/gwas/bb/bb_ALP.p03', header=T) %>% dplyr::select(SNP, CHR, POS, P) %>% filter(!is.na(SNP)) %>% mutate(P=ifelse(P<1e-20,1e-20,P)) %>% rename(ALP=P)
	dat2 <- read.table('D:/data/gwas/bb/covid_severe.p03', header=T) %>% 
	dplyr::select(rsid, CHR, POS, all_inv_var_meta_p) %>% filter(!is.na(rsid)) %>% mutate(all_inv_var_meta_p=ifelse(all_inv_var_meta_p<1e-20,1e-20,all_inv_var_meta_p)) %>% rename(SNP=rsid, CHR.2=CHR, POS.2=POS, COVID_severe=all_inv_var_meta_p)
	dat3 <- read.table('D:/data/gwas/bb/covid_inpatient.p03', header=T) %>% 
	dplyr::select(rsid, CHR, POS, all_inv_var_meta_p) %>% filter(!is.na(rsid)) %>% mutate(all_inv_var_meta_p=ifelse(all_inv_var_meta_p<1e-20,1e-20,all_inv_var_meta_p)) %>% rename(SNP=rsid, CHR.3=CHR, POS.3=POS, COVID_inpatient=all_inv_var_meta_p)
	dat4 <- read.table('D:/data/gwas/bb/covid_infected.p03', header=T) %>% 
	dplyr::select(rsid, CHR, POS, all_inv_var_meta_p) %>% filter(!is.na(rsid)) %>% mutate(all_inv_var_meta_p=ifelse(all_inv_var_meta_p<1e-20,1e-20,all_inv_var_meta_p)) %>% rename(SNP=rsid, CHR.4=CHR, POS.4=POS, COVID_infected=all_inv_var_meta_p)
dat0 <- Reduce(function(x,y) merge(x,y,by='SNP',all=T, no.dups=T), list(dat1, dat2, dat3, dat4)) # !! Must Remove NA SNPs before merging
dat <- dat0 %>%
	mutate (chrom=apply( dat0[,grepl('CHR',names(dat0))], 1, FUN=min, na.rm=T ), position=apply( dat0[,grepl('POS',names(dat0))], 1, FUN=min, na.rm=T )) %>% 
	dplyr::select(SNP, chrom, position, COVID_severe, COVID_inpatient, COVID_infected, ALP)
	write.table(dat, "plot.dat.txt", quote=F, row.names=F, append=F)
SNPs <- subset(dat, (chrom==1 & position >=21835857 & position <=21904905) | (chrom==9 & position >=136130562 & position<=136150630))$SNP
CMplot(dat, type="p", plot.type="c", r=0.4, col=matrix(c("gold","orange", "gray","lightblue", "olivedrab3", "orange", "grey30","grey60"), nrow=4, byrow=T),
	highlight=SNPs, highlight.col="red",highlight.cex=1, highlight.pch=10,
	threshold=5e-8, cir.chr.h=1, amplify=F, threshold.col="red", signal.line=1, signal.col="black",
	bin.size=1e6, outward=T, file='jpg', memo='', dpi=600, file.output=T, verbose=T, width=12, height=12
)