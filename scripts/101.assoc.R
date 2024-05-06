setwd("D:/")
pacman::p_load(data.table, readxl, lubridate, tidyverse, dplyr, survival, survminer, reshape2, psych)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))

dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata") %>% mutate(birth_month=factor(birth_month))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X-Y 或 X-Y-Z交互作用 批量分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% drop_na(age, sex) %>% filter(ethnic_cat=="White")
Ys <- grep("^icdDate_", names(dat), value=TRUE)
Xs <- grep("^age_sex|age_m|^edu_score|^birth_weight|birth_month|^height$|^chunk|^leg|^hippo_|^fev1fvc|^stiffness|score_sum$", names(dat), value=TRUE) 
Zs <- grep("^o$|^se$", names(dat), value=TRUE) # |^rh|shbg|^apoe$|\\.rs
outfile="101.assoc.tsv"; file.create(outfile); #sink("101.assoc.log")
for (Y in Ys) {
	writeLines(paste('\n\n--> Run:', Y))
	dat1 <- dat %>%
	mutate(
		Y_date = dat[[Y]],
		Y_yes = ifelse( is.na(Y_date), 0,1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0 )
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	for (X in Xs) {
		if (X==Y) next
		print(paste("RUN", X, Y))
		dat1$X = inormal(dat1[[X]])
		fit.cox <- coxph(surv.obj ~ X +age+sex +PC1+PC2+weight+smoke_status+alcohol_status, data=dat1)
		print(res.cox <- coef(summary(fit.cox)))
		b=round(res.cox[1,1],4); se=round(res.cox[1,3],4); p=signif(res.cox[1,5],2)
		write.table(paste(Y, X, b, se, p), file=outfile, sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
		next
		dat1$X_qt <- cut(dat1$X, breaks=quantile(dat1$X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5))
		dat1$X_qt <- factor(ifelse(dat1$X_qt=="q1", "low", ifelse(dat1$X_qt=="q5", "high", "middle")), levels=c("low","middle","high"))
		for (Z in Zs) {
			print(paste("RUN", X, Y, Z))
			dat1$Z <- dat1[[Z]]
			if (Z %like% "\\.rs" & length(unique(dat1$Z)) >3) dat1$Z <- hardcall(dat1$Z) # 要不然没法画出Z的3个类型
			fit.cox <- coxph(surv.obj ~ X + Z + X*Z +age+sex +PC1+PC2, data=dat1)
			res.cox <- coef(summary(fit.cox)); print(res.cox)
			png(paste(X, Y, Z, "png", sep="."))
			print(ggsurvplot(survfit(surv.obj ~X_qt + Z, data=dat1), ylim=c(0.5,1), risk.table=FALSE))	
			dev.off()
		}
	}
}
#sink()
# 可用下面代码 transpose 汇总数据，画 heatmap 图
dat <- read.table('D:/analysis/ldsc/all.rg.res', header=T) 
dat$bp <- paste(round(dat$rg,3), "(", signif(dat$p,2), ")", sep="")
reshape_bp <- dat %>% reshape2::acast(p2 ~ p1, value.var='bp')
write.table(reshape_bp, file="101.new.tsv", sep='\t', row.names=TRUE, col.names=TRUE, append=FALSE, quote=FALSE)
library(ggcorrplot); plot_bp <- ggcorrplot::ggcorrplot(dat_b, lab=TRUE, p.mat=dat_p, sig.level=1e-4, insig ='blank') + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 针对某一*显著*结果的精细分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- dat %>% 
	mutate(
		walk_pace=4 - walk_pace, 
		across(grep("walk", names(dat0), value=T), ~factor(.x)),
		walk_pace=factor(walk_pace, labels=c("brisk","steady", "slow"))
	) %>%
	rename(X=walk_pace, Y_date=icdDate_vte, Z=vte.nf.score_sum) %>% drop_na(X, Z) 
dat1 <- dat1 %>% 
	mutate(
	#	X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		Z_qt = cut(Z, breaks=quantile(Z, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
	#	X_qt = factor(ifelse(X_qt=="q1", "loss", ifelse(X_qt=="q5", "gain", "same")), levels=c("loss", "same", "gain")),
		Z_qt = factor(ifelse(hardcall(vte.F5.rs6025_C)!=2, "F5", ifelse(vte.F2.rs1799963_G !=2, "F2", ifelse(Z_qt=="q1", "low", ifelse(Z_qt=="q5", "high", "middle")))), levels=c("low", "middle", "high", "F2", "F5")),
		X_Z = factor(paste(Z_qt, X, sep="|"), levels=c("low|brisk","low|steady","low|slow", "middle|brisk","middle|steady","middle|slow", "high|brisk","high|steady","high|slow",  "F2|brisk","F2|steady","F2|slow", "F5|brisk","F5|steady", "F5|slow")),
		Y_yes = ifelse(is.na(Y_date), 0, 1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0 )			
coef(summary(glm(Y_yes ~ X + age+sex+smoke_status+alcohol_status +PC1+PC2, data=dat1)))
	aggregate(Y_yes ~ X_qt, dat1, FUN=function(x) { paste( length(x), sum(x), round(sum(x)/length(x),3)) } )
surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
km.obj <- survfit(surv.obj ~ Z_qt, data=dat1)
	plot(km.obj, ylim=c(0.5,1)); plot(km.obj, fun=function(x) 1-x)
	ggsurvplot(km.obj, ylim=c(0,0.1), fun="event", break.time.by=2, risk.table=FALSE, censor=FALSE, pval=TRUE) #palette=c("") 
fit.cox <- coxph(surv.obj ~ X_Z +walk_time+walk_freq +age+bmi+PC1+PC2+ smoke_status+alcohol_status, data=dat1); print(coef(summary(fit.cox)))
	survminer::ggforest(fit.cox, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1) # 不能显示interaction值
	fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
		fit.glm <- glm(outcome_yes ~ gen, data=dat, family="binomial")
	survdiff(surv.obj ~ gen_3p, data=dat) # log-rank test


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 高阶分析和画图
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
pacman::p_load(vivid, randomForest, MASS)
set.seed(12345)
dat2 <- dat1 %>% drop_na(X, Z, walk_time,walk_freq, age, sex, bmi,PC1,PC2,smoke_status,alcohol_status)
rf_model <- randomForest(Y_yes ~ X+Z+ walk_time+walk_freq +age+sex+bmi+PC1+PC2+ smoke_status+alcohol_status, data=dat2)
VIVI_rf <- vivi(fit = rf_model, response = "Y_yes", data=dat2)
viviHeatmap(mat = VIVI_rf)
viviNetwork(mat = VIVI_rf)