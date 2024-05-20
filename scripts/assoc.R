setwd("D:/")
pacman::p_load(data.table, lubridate, dplyr, tidyr, survival, survminer, gtsummary, psych)
inormal <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))

dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata")
dat <- dat0 %>% filter(ethnic_cat=="White") 
dat1 <- dat %>% rename(X=walk_pace, Y_date=icdDate_vte2, F5=vte.F5.rs6025_C, F2=vte.F2.rs1799963_G) %>% 
	drop_na(X, F2, F5) 
dat1 <- dat1 %>% mutate(
	#	X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
	#	X_qt = factor(ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")), levels=c("low", "middle", "high")),
	#	Z = factor(ifelse(hardcall(F5)!=2, "F5", ifelse(F2 !=2, "F2", ifelse(Z_qt=="q1", "low", ifelse(Z_qt=="q5", "high", "middle")))), levels=c("low", "middle", "high", "F2", "F5")),
		Z = factor(ifelse(hardcall(F5)!=2, "F5", ifelse(hardcall(F2) !=2, "F2", "none")), levels=c("none","F2","F5")),
		X_Z = factor(paste(Z, X, sep="|"), levels=paste(rep(levels(Z),each=3), rep(levels(X),3), sep="|")),
		Y_yes = ifelse(is.na(Y_date), 0, 1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0)
	# dat2 <- dat1 %>% dplyr::select(eid, follow_years, Y_yes, X, Z, X_Z, walk_time,walk_freq,age,sex,bmi,PC1,PC2,smoke_status,alcohol_status) 
prop.table(table(dat1$X, dat1$Y_yes),1); prop.table(table(dat1$X_Z, dat1$Y_yes),1)	
	aggregate(Y_yes ~ X, dat1, FUN=function(x) { paste( length(x), sum(x), round(sum(x)/length(x),3)) } )
	coef(summary(glm(Y_yes ~ X+Z+ age+sex+smoke_status+alcohol_status +PC1+PC2, family="binomial", data=dat1)))
surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	km.obj <- survfit(surv.obj ~ Z, data=dat1) # KM 是不能带协变量的，Z会被作为分层变量
	survdiff(surv.obj ~ Z, data=dat1) # log-rank test
	plot(km.obj, fun=function(x) 1-x)
	survminer::ggsurvplot(km.obj, ylim=c(0,0.08), fun="event", pval=TRUE, risk.table=TRUE, risk.table.col="Z", ncensor.plot = TRUE, ggtheme = theme_bw(), palette=c("green","gray","orange"))
cox.fit <- coxph(surv.obj ~ X+Z+X*Z +walk_time+walk_freq +age+sex+bmi+PC1+PC2+ smoke_status+alcohol_status, data=dat1); print(coef(summary(cox.fit)))
	survminer::ggforest(cox.fit, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1) # 不能显示interaction值
	cox.fit %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
# 计算10年风险并画图 🏮
cox.fit <- coxph(surv.obj ~ X+Z, data=dat1) 
	dat2 <- expand.grid(X=levels(dat1$X), Z=levels(dat1$Z))
	surv.pred <- survfit(cox.fit, newdata=dat2)
	surv.10y <- summary(surv.pred, times=10)
dat2 <- dat2 %>% mutate(
	risk = 1-t(surv.10y$surv), ci_lower = t(1-surv.10y$upper), ci_upper = t(1-surv.10y$lower)
)
ggplot(dat2, aes(x=Z, y =risk, fill=X)) +
	geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7) +
	geom_errorbar(aes(ymin =ci_lower, ymax=ci_upper), width=0.2, position=position_dodge(width=0.7)) +
	geom_text(aes(label=sprintf("%.2f", risk), y=ci_upper), vjust=-1, position=position_dodge(width=0.7), size=2) +
	labs(x="Mendelian mutations", y="10-year risk of VTE (%)", title="") +
	scale_fill_manual(values=c("green", "gray", "orange"), name="Walk pace") + theme_minimal() #+ coord_cartesian(ylim=c(0.915, 0.99))