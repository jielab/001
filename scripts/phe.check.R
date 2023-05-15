setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, lubridate, dplyr, ggplot2, tidyverse, naniar, survival, ggsurvfit, gtsummary, survminer, multcomp, scatterplot3d, mediation)

dat0 <- readRDS(dat0, file="D:/data/ukb/Rdata/all.Rdata")
dates <- names(dat0)[sapply(dat0, is.Date)]; summary(dat0[,dates]) # check invalid dates


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sanity check by BMI ~ FTO
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% rename(bmi.prs1=bmi.EUR77.score_sum, bmi.prs2=bmi.EUR941.score_sum, bmi.prs3=bmi.EUR2446.score_sum, fto=bmi.FTO.rs9939609_T,
			height.prs1=height.EUR697.score_sum, height.prs2=height.EUR3290.score_sum, height.snp1=height.ZBTB38.rs724016_A, height.snp2=height.ZBTB38.rs724016_A) %>%
	filter( ethnic_cat=="White" & is.na(related) & aneuploidy==0 & excess.relatives==0 ) # ethnicity_gen==1 
cor(dat$fto, dat$bmi, use="complete.obs")
aggregate(bmi ~ fto, data=dat, mean); # bmi.prs1/2/3  height.prs1/2  height.snp1/2
coplot(dat$bmi.prs1 ~ dat$bmi | dat$sex)
lmfit <- lm(bmi ~ fto, data=dat, na.action=na.exclude)
summary(lmfit); summary(lmfit)$coef; coef(lmfit); confint(lmfit, level=0.95)
dat$bmi_pred = predict.lm(lmfit); table(dat$bmi_pred, dat$fto)
round(prop.table(table(dat$bmi_cat, dat$fto, dnn=c("bmi","FTO")), 2),4)
round(prop.table(ftable(table(dat$sex, dat$bmi_cat, dat$fto, dnn=c("sex","bmi","FTO"))), 1),4)
lmsum <- summary(lmfit <- lm(bmi ~ age +sex + as.factor(fto), data=dat)); signif(lmsum$coef, 3) 
plot(dat$bmi ~ dat$fto)
boxplot(dat$bmi ~ dat$fto)
plot(dat$fto + runif(nrow(dat), -0.1, 0.1), dat$bmi)
ggline(dat, x="fto", y="bmi", xlab="FTO", ylab="BMI", combine=T, add="mean_se", color="sex", palette="jco", scales="free_y") +rotate_x_text(45)
with(dat, scatterplot3d::scatterplot3d(x=age, y=fto, z=bmi, xlab="Age", ylab="FTO", zlab="BMI", colvar=bmi))
dat$prs <- dat$bmi.prs3 # make a nice bmi-PRS correlation plot
	dat1 <- subset(dat, !is.na(bmi) & !is.na(prs))
	dat1 <- subset(dat1, prs <=mean(prs)+3*sd(prs) & prs >=mean(prs)-3*sd(prs))
	myhist = hist(dat1$prs, breaks=20) # labels=T
	avg = by(dat1$bmi, cut(dat1$prs, breaks=myhist$breaks), function(x) mean(x,na.rm=T))
	par(new=T); plot(myhist$mids, avg, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2)
	axis(side=4); mtext(side=4, line=3, 'measured')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check Dementia (PMID 31302669, JAMA 2019) results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% rename(outcome_date=dementia_date, exposure=e4_f, gen=dementia.score_sum) %>%
	filter( ethnicity_gen==1 & age >=60 )
dat <- dat %>% 
	mutate( 
		outcome_yes = ifelse( is.na(outcome_date), 0,1),
		follow_end_day = data.table::fifelse(!is.na(outcome_date), outcome_date, fifelse(!is.na(death_date), death_date, as.Date("2018-02-28"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25
	) %>% filter( follow_years >0 ) %>% 
	mutate (
		gen_5p = cut(gen, breaks=quantile(gen, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		gen_3p = ifelse(gen_5p=="q1", "low", ifelse(gen_5p=="q5", "high", "middle")),
		gen_3p = factor(gen_3p, levels=c("low","middle","high"))
	)
summary(dat$outcome_date); table(dat$outcome_yes, useNA="always"); hist(dat$follow_years); sum(dat$follow_years) #1,545,433 in paper
prop.table(table(dat$gen_3p, dat$outcome_yes), 1)  # 1.23% vs. 0.63% in paper
	summary(glm(outcome_yes ~ age + sex + gen_3p, data=dat))
	ggplot(dat, aes(x=exposure, y=gen, fill=apoe)) + geom_bar(position="dodge", stat="identity")
	bp <- boxplot(gen ~ exposure*apoe, data=dat, xlab="", ylab="", main="", las=2, col=rainbow(6), font=2); bp$stats
	aggregate(gen ~ exposure*apoe, data=dat, FUN=function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2)} )
surv.obj <- Surv(time=dat$follow_years, event=dat$outcome_yes); surv.obj[1:100]
fit.1 <- survival::survfit(surv.obj ~ smoke_status, data=dat); summary(fit.1, data=dat, times=1) # ~ 1; ~ apoe
	survminer::ggsurvplot(fit.1, ylim=c(0.95,1), data=dat) 
	survdiff(surv.obj ~ gen_3p, data=dat) # log-rank test conducts between-group significance tests
	fit.1 %>% gtsummary::tbl_survfit(times=1, label_header = "**1-year survival (95% CI)**")
	dat %>% filter(outcome_yes==1) %>% summarize(median_surv=median(follow_years)) # when censoring is not considered
fit.2 <- ggsurvfit::survfit2(surv.obj ~ gen_3p, data=dat) 
	fit.2 %>% ggsurvfit() + labs(x="Days", y = "Survival probability") + add_confidence_interval() + add_risktable()
fit.cox <- coxph(surv.obj ~ age+sex + exposure + gen + exposure * gen, data=dat); summary(fit.cox) # quantify effect size
	fit.cox %>% gtsummary::tbl_regression(exp=T) 
	ggforest(fit.cox, data=dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check Arterial Stiffness (PMID 35076696, Diabetes care 2022) results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% rename(outcome_date=t2d_date, exposure=stiffness, gen=t2d.score_sum) %>% # stiffness
	filter( ethnic_cat=="White" & !is.na(exposure) ) # N=152,611 with ASI, excluded diabetes (n=26,736) or CVD (n=21,652) at baseline 
hist(dat$outcome_date, "years")
hist(dat$attend_date, "months")
dat <- dat %>% 
	mutate(
		outcome_yes = ifelse( is.na(outcome_date), 0,1),
		follow_end_day = fifelse(!is.na(outcome_date), outcome_date, fifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25,
	) %>% filter( follow_years >0 )
table(dat$outcome_yes) # 9.5 years, 3,000 developed T2D
	summary(glm(outcome_yes ~ age + sex + exposure + gen + exposure * gen, data=dat))
surv.obj <- Surv(time=dat$follow_years, event=dat$outcome_yes)
	#ggsurvplot(survfit(surv.obj ~ exposure_3p, data=dat), ylim=c(0.9,1)) 
	summary(cox.fit <- coxph(surv.obj ~ age+sex + exposure + gen + exposure * gen, data=dat))
	ggforest(cox.fit, data=dat) # ASI / 1-SD -> 3% T2D risk; HR=1.58 highest vs. lowest ASI; interaction betweeen ASI and PRS


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# find trio
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ped <- subset(dat0, race=="White", select=c("IID", "age", "sex", "birth_year", "race"))
kin <- read.table("D:/files/ukb.kin0", header=T, as.is=T)
kin <- subset(kin, ID1>0 & InfType=="PO", select=c("ID1", "ID2"))
dat <- merge(ped, kin, by.x="IID", by.y="ID1")
dat <- merge(ped, dat, by.x="IID", by.y="ID2")
summary(abs(dat$age.x - dat$age.y));  dat <- subset(dat, abs(dat$age.x - dat$age.y)> 18)
colnames(dat)[colnames(dat)=="IID"] <- "IID.x"
dat$child <- ifelse(dat$age.x > dat$age.y, dat$IID.y, dat$IID.x)
dat$father <- ifelse(dat$IID.x==dat$child & dat$sex.y==1, dat$IID.y, ifelse(dat$IID.y==dat$child & dat$sex.x==1, dat$IID.x, NA))
dat$mother <- ifelse(dat$IID.x==dat$child & dat$sex.y==0, dat$IID.y, ifelse(dat$IID.y==dat$child & dat$sex.x==0, dat$IID.x, NA))
father <- subset(dat, !is.na(father), select=c("child", "father"))
mother <- subset(dat, !is.na(mother), select=c("child", "mother"))
trio <- merge(father, mother)
### merge trio phenotypes 
phe <- subset(dat0, select=c("IID", "age", "sex", "bmi","height", paste0("PRS.",0:8), paste0("PC",1:40)))
dat <- merge(trio, phe, by.x="child", by.y="IID")
dat <- merge(dat, phe, by.x="father", by.y="IID") # ".x" for child, ".y" for father
dat <- merge(dat, phe, by.x="mother", by.y="IID") %>% rename(IID=child, IID.f=father, IID.m=mother)
dat <- na.omit(dat) # %>% filter_at(grep("height|prs", names(dat),value=T), all_vars(!is.na(.)))
dat <- subset(dat, select=c("IID","IID.f","IID.m", outer(c("age","sex","bmi","height", paste0("PRS.",0:8), paste0("PC",1:40) ), c("",".x",".y"), FUN=paste0)))
names(dat) <- c("IID","IID.f","IID.m", outer(c("age","sex","bmi","height", paste0("PRS.",0:8), paste0("PC",1:40)), c(".m","",".f"), FUN=paste0))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 大批量 pair-wise associations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter( aneuploidy==0 & excess.relatives==0 & ethnic_cat=="White" & is.na(related)) 
dat <- dat %>%
	mutate(
		age_mn_3p = ifelse(age_mn %in% 12:14, "normal", ifelse(age_mn>14, "late", "early")),
		age_mn_3p = factor(age_mn_3p, levels=c("normal", "early", "late")),
		sp1_s=as.factor(ifelse(grepl("SS", sp1), 2,  ifelse(grepl("S", sp1), 1, 0))),
		sp1_z=as.factor(ifelse(grepl("ZZ", sp1), 2,  ifelse(grepl("Z", sp1), 1, 0))), # Use "grepl" NOT "grep"
	)
Yvars = grep("inf$|^date_|^fod_",names(dat), value=T); Yvars=Yvars[Yvars !="date_attend"]
Xvars = grep("^age_mn|^abo1$|^muscle_|^apoe$|^e4|^sp1_|^hla_|^bb_|^bc_", names(dat), value=T)
Yvars = grep("inf$|^date_|^fod_",names(dat), value=T); Yvars=Yvars[Yvars !="date_attend"]
Xvars = grep("abo1|^sp1_", names(dat), value=T)
df = NULL
for (Yvar in Yvars) {
	print(Yvar)
	for (Xvar in Xvars) {
		if (is.Date(dat[[Yvar]])) {
			dat1 <- dat %>% mutate(
			Xvar = dat[[Xvar]],
			outcome_date = dat[[Yvar]],
			outcome_yes = ifelse( is.na(outcome_date), 0,1),
			follow_end_day = fifelse(!is.na(outcome_date), outcome_date, fifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
			follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25		
			) %>% filter(follow_years >0)
			surv.obj <- Surv(time=dat1$follow_years, event=dat1$outcome_yes)
			fit = coxph(surv.obj ~ Xvar +age+sex+PC1+PC2+PC3+PC4, data=dat1)
			p_value = summary(fit)$coef[1,5]
		} else if (n_distinct(dat[[Yvar]],na.rm=T) > 2) {
			fit = lm(dat[[Yvar]] ~ dat[[Xvar]] +age+sex+PC1+PC2+PC3+PC4, data=dat)
			p_value = summary(fit)$coef[2,4]
		} else if (n_distinct(dat[[Yvar]],na.rm=T) ==2){
			fit = glm(dat[[Yvar]] ~ dat[[Xvar]] +age+sex+PC1+PC2+PC3+PC4, data=dat, family=binomial)
			p_value = summary(fit)$coef[2,4]
		} else {
			next
		}
		p_value=signif(p_value, 2)
		print(paste(" ->>", Xvar, p_value, sep=", "))
		df = rbind(df, paste(Yvar, Xvar, p_value, sep=", "))
	}
}
df2 = reshape(as.data.frame(df), idvar="V1", timevar="V2", direction="wide")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 五个 CVD 亚型
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cvd_names = c("rheumatic", "hypertensive", "ischaemic", "pulmonary", "other", "cerebro", "artery", "vein", "unspecified")	
dat$cvd_cnt = rowSums(dat[, paste0("icd_",cvd_names)], na.rm=T)
for (trait in cvd_names) {
	dat$tmp = dat[[paste0("icd_",trait)]]
	dat[[paste0(trait,"_raw")]] = ifelse(dat$icd_cvd==0, 0, ifelse(dat$tmp==1, 1, NA))
	dat[[paste0(trait,"_1")]] = ifelse(dat$icd_cvd==0, 0, ifelse(dat$tmp==1 & dat$cvd_cnt==1, 1, NA))
	dat[[paste0(trait,"_2")]] = ifelse(dat$icd_cvd==0, 0, ifelse(dat$tmp==1 & (dat$cvd_cnt==1 | (dat$cvd_cnt==2 & dat$icd_hypertensive==1)), 1, NA))
}
#venn.diagram(list(MI=dat$IID[!is.na(dat$date_MI)], STEMI=dat$IID[!is.na(dat$date_STEMI==1)], NSTEMI=dat$IID[!is.na(dat$date_NSTEMI==1)]), scaled=T, imagetype="tiff", filename="before.png", fill=rainbow(3), main.cex=2)
#venn.diagram(list(MI=na.omit(dat$IID[dat$MI==1]), STEMI=na.omit(dat$IID[dat$STEMI==1]), NSTEMI=na.omit(dat$IID[dat$NSTEMI==1])), scaled=T, imagetype="png", filename="after.png", fill=rainbow(3), main.cex=2)
dat1 <- subset(dat, select=paste0("icd_",cvd_names)); dat1[is.na(dat1)] =0; names(dat1) = cvd_names
corrplot(cor(dat1), col=terrain.colors(100), cl.pos="b", tl.pos="d", tl.srt=60)
corrplot.mixed(cor(dat1), number.cex=1, number.digits=2, tl.cex=1, cl.cex=1, cl.lim=c(0,1))
dat1.mat <- matrix(as.numeric(as.matrix(dat1)), nrow = nrow(dat1))
dat1.sum <- crossprod(dat1.mat)
rownames(dat1.sum) <- cvd_names; colnames(dat1.sum) <- cvd_names
dat1.sum # sum displayed in the diag line
colsum <- apply(dat1.mat, 1, sum); dat1.mat[colsum !=1, ] =0
rowsum <- apply(dat1.mat, 2, sum); diag(dat1.sum) <- rowsum
rownames(dat1.sum) <- cvd_names; colnames(dat1.sum) <- cvd_names
dat1.sum # uniq displayed in the diag line
dat1 = subset(dat, select=c("IID","IID", "copd01", "asthma01", "asthma_dx", "allergy_dx", "MI","STEMI","NSTEMI","otherMI", "icd_cvd", "icd_sca", "icd_arrhythmias", paste0(cvd_names,"_raw"),  paste0(cvd_names,"_1"),  paste0(cvd_names,"_2")))
dat1$hypertensive_2 = NULL
write.table(dat1, paste("CVD.pheno",sep="."), na="NA", append=F, quote=F, row.names=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stratify by sex vs. Adjust for sex
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# height>=140 & height <=210 & height_sitting>=70 & height_sitting <=110
dat$trait=dat$height
aggregate(trait ~ sex, data=dat, mean)
mu = ddply(dat, "sex", summarise, grp.mean=mean(trait))
h = ggplot(dat, aes(x=trait)) + geom_histogram(aes(color=sex, fill=sex), position="identity", bins=30, alpha=0.4) 
h + geom_vline(data=mu, aes(xintercept=grp.mean, color=sex), linetype="dashed") + theme(legend.position="top")
ggplot(dat1, aes(x=trait)) + geom_histogram(aes(color=sex), fill = "white", position = "identity")	
m = subset(dat, sex=="m"); f = subset(dat, sex=="f"); boxplot(m$trait, f$trait)
plot(density(m$trait), xlab="", col="blue", main=paste0("N=",nrow(dat)))
lines(density(f$trait), xlab="", col="red", main="")
legend("topright", c("male", "female"), cex=1.5, col=c("blue", "red"), lty=1, pch=1)
dat$res_AgeSex <- residuals(lm(trait ~ age + sex, data=dat, na.action=na.exclude))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate files for GWAS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- subset(dat, select=grepl("sex|bb_|bc_", names(dat))); gg_miss_var(dat1, facet=sex)
dat1 <- subset(dat, select=grep("rheu|oest", grep("IID|^age$|sex|^PC|bmi|bc_", names(dat), value=T), invert=T, value=T)); 
dat_all <- subset(dat, select=c("IID", "sex"))
for (t in grep("^bmi$|height$|height_sitting|leg|height_r|bb_|bc_", names(dat), value=T)) {
	dat$trait = dat[[t]]
	dat1 = subset(dat, !is.na(trait))
	pdf(paste0(t,".pdf"), h=12, w=8); par(mfrow=c(3,2))
	m = subset(dat1, sex=="m"); f = subset(dat1, sex=="f")
	ggplot(dat1, aes(x=trait)) + geom_histogram(aes(color=sex), fill = "white", position = "identity")	
	boxplot(m$trait, f$trait)
	plot(density(m$trait), xlab="", col="blue", main=paste0("N=",nrow(dat1)))
	lines(density(f$trait), xlab="", col="red", main="")
	legend("topright", c("male", "female"), cex=1.5, col=c("blue", "red"), lty=1, pch=1)
	df <- NULL
	for (sex in c("m","f")) {
		dat2 = get(eval(sex)) %>%
	    mutate(
		trait_res = round(residuals(lm(trait ~ age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, na.action=na.exclude)),5),
		trait_inv = round(qnorm((rank(trait_res,na.last="keep")-0.5)/length(na.omit(trait_res))), 5)
		) %>% select(IID,trait_inv) %>% set_colnames(c("IID", t)) # setnames(dat2, old=c("IID","IID.1"), new=c("FID","IID"))
		df <- rbind(df, dat2)
		write.table(dat2, paste(t,sex,"pheno.txt",sep="."), na="NA", append=F, quote=F, row.names=F)
		hist(dat2$trait, xlab="", main=paste(t,sex,"raw",nrow(dat2),sep=", "))
		hist(dat2$trait_inv, xlab="", main=paste(t,sex,"inv-norm",nrow(dat2),sep=", "))
		dat_all = merge(dat_all, df, by="IID", all=T)
	}
	dev.off()
	print(t)
}
write.table(dat_all, "bbc.pheno", na="NA", append=F, quote=F, col.names=T, row.names=F)

