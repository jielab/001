setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, lubridate, dplyr, tidyverse, ggiraphExtra, ggplot2, scatterplot3d, survival, survminer)

dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata")
dates <- names(dat0)[sapply(dat0, is.Date)]; summary(dat0[,dates]) # check invalid dates


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sanity check by BMI ~ FTO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% rename(bmi.prs1=bmi.EUR77.score_sum, bmi.prs2=bmi.EUR941.score_sum, bmi.prs3=bmi.EUR2446.score_sum, fto=bmi.FTO.rs9939609_T,
			height.prs1=height.EUR697.score_sum, height.prs2=height.EUR3290.score_sum, height.snp1=height.ZBTB38.rs724016_A, height.snp2=height.ZBTB38.rs724016_A) %>%
	filter( ethnicity_gen==1 & is.na(related) & aneuploidy==0 & excess.relatives==0 ) %>%
	mutate(sex.fa=as.factor(sex))
cor(dat$fto, dat$bmi, use="complete.obs") # 一些简单的分组统计
	aggregate(gbmi ~ fto, data=dat, FUN=function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2)} )
	round(prop.table(table(dat$bmi_cat, dat$fto, dnn=c("bmi","FTO")), 2),4)
	round(prop.table(ftable(table(dat$sex, dat$bmi_cat, dat$fto, dnn=c("sex","bmi","FTO"))), 1),4)
fit.lm <- lm(bmi ~ bmi.prs3 + sex.fa, data=dat, na.action=na.exclude) # 一些简单的线性回归和预测
	summary(fit.lm)$coef; signif(coef(fit.lm),3); confint(fit.lm, level=0.95)
	dat$bmi_pred = predict.lm(fit.lm); table(dat$bmi_pred, dat$fto)
	pred.lm <- ggiraphExtra::ggPredict(fit.lm, terms="bmi.prs3")
plot(pred.lm, se=T, interactive=T) # 一些简单的图
	plot(dat$fto + runif(nrow(dat), -0.1, 0.1), dat$bmi)
	coplot(dat$bmi.prs1 ~ dat$bmi | dat$sex)
	with(dat, scatterplot3d::scatterplot3d(x=age, y=fto, z=bmi, xlab="Age", ylab="FTO", zlab="BMI", colvar=bmi))
dat1 <- subset(dat, select=c("bmi", "bmi.prs3")) %>% rename(prs=bmi.prs3) %>% na.omit() # 下面叠加两个图
	myhist = hist(dat1$prs, breaks=20) 
	avg = by(dat1$bmi, cut(dat1$prs, breaks=myhist$breaks), function(x) mean(x,na.rm=T))
	par(new=T); plot(myhist$mids, avg, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2)
	axis(side=4); mtext(side=4, line=3, 'measured')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check Dementia (PMID 31302669, JAMA 2019) results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% rename(outcome_date=dementia_date, exposure=smoke_status, gen=dementia.score_sum) %>%
	filter( ethnicity_gen==1 & age >=60 ) %>%
	mutate( 
		outcome_yes = ifelse( is.na(outcome_date), 0,1),
		follow_end_day = data.table::fifelse(!is.na(outcome_date), outcome_date, fifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25
	)
	table(dat$outcome_yes); nrow(subset(dat, outcome_date <= attend_date))
	fit.glm <- glm(outcome_yes ~ gen, data=dat, family="binomial")
	ggPredict(fit.glm, se=T, interactive=T)
dat <- dat %>% 
	filter( follow_years >0 ) %>% 
	mutate (
		gen_5p = cut(gen, breaks=quantile(gen, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		gen_3p = ifelse(gen_5p=="q1", "low", ifelse(gen_5p=="q5", "high", "middle")),
		gen_3p = factor(gen_3p, levels=c("low","middle","high"))
	)
	
	table(dat$outcome_yes, useNA="always"); hist(dat$follow_years); sum(dat$follow_years) #1,545,433 in paper
	prop.table(table(dat$gen_3p, dat$outcome_yes), 1)  # 1.23% vs. 0.63% in paper
surv.obj <- Surv(time=dat$follow_years, event=dat$outcome_yes)
fit.suv <- survival::survfit(surv.obj ~ exposure, data=dat) # compute a survival curve
	summary(fit.suv, data=dat, times=1) 
	survminer::ggsurvplot(fit.suv, ylim=c(0.95,1), data=dat) 
	survdiff(surv.obj ~ gen_3p, data=dat) # log-rank test
fit.cox <- coxph(surv.obj ~ exposure + gen + age+sex, data=dat) # quantify effect size
	summary(fit.cox) 
	ggforest(fit.cox, data=dat)


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