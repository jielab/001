setwd("D:/")
pacman::p_load(data.table, readxl, lubridate, tidyverse, dplyr, survival, survminer, gtsummary, reshape2, psych)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))

dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata") %>% mutate(birth_month=factor(birth_month))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X-Y 或 X-Y-Z交互作用 批量分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% mutate(
		icdDate_vte2=as.Date(ifelse(!is.na(icdDate_vte), as.character(icdDate_vte), ifelse(!is.na(srdYear_vte), (paste0(srdYear_vte,"-07-01")), ifelse(!is.na(srdAge_vte) & srdAge_vte >0, paste0(birth_year+srdAge_vte,"-07-01"), NA))))
	) %>% filter(ethnic_cat=="White") 
hist(dat$icdDate_vte2, breaks="year")
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
		cox.fit <- coxph(surv.obj ~ X +age+sex +PC1+PC2+weight+smoke_status+alcohol_status, data=dat1)
		print(res.cox <- coef(summary(cox.fit)))
		b=round(res.cox[1,1],4); se=round(res.cox[1,3],4); p=signif(res.cox[1,5],2)
		write.table(paste(Y, X, b, se, p), file=outfile, sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
		next
		dat1$X_qt <- cut(dat1$X, breaks=quantile(dat1$X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5))
		dat1$X_qt <- factor(ifelse(dat1$X_qt=="q1", "low", ifelse(dat1$X_qt=="q5", "high", "middle")), levels=c("low","middle","high"))
		for (Z in Zs) {
			print(paste("RUN", X, Y, Z))
			dat1$Z <- dat1[[Z]]
			if (Z %like% "\\.rs" & length(unique(dat1$Z)) >3) dat1$Z <- hardcall(dat1$Z) # 要不然没法画出Z的3个类型
			cox.fit <- coxph(surv.obj ~ X + Z + X*Z +age+sex +PC1+PC2, data=dat1)
			res.cox <- coef(summary(cox.fit)); print(res.cox)
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
dat1 <- dat %>% mutate(
		walk_pace=4 - walk_pace, 
		across(grep("walk", names(dat0), value=T), ~factor(.x)),
		walk_pace=factor(walk_pace, labels=c("brisk","steady", "slow"))
	) %>% rename(X=walk_pace, Y_date=icdDate_vte2, F5=vte.F5.rs6025_C, F2=vte.F2.rs1799963_G) %>% 
	drop_na(X, F2, F5) 
dat1 <- dat1 %>% mutate(
	#	X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
	#	X_qt = factor(ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")), levels=c("low", "middle", "high")),
	#	Z = factor(ifelse(hardcall(F5)!=2, "F5", ifelse(F2 !=2, "F2", ifelse(Z_qt=="q1", "low", ifelse(Z_qt=="q5", "high", "middle")))), levels=c("low", "middle", "high", "F2", "F5")),
		Z = factor(ifelse(hardcall(F5)!=2, "F5", ifelse(hardcall(F2) !=2, "F2", "none")), levels=c("none","F2","F5")),
		Z_X = factor(paste(Z, X, sep="|"), levels=paste(rep(levels(Z),each=3), rep(levels(X),3), sep="|")),
		Y_yes = ifelse(is.na(Y_date), 0, 1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0)
	#dat2 <- dat1 %>% dplyr::select(eid, follow_years, Y_yes, X, Z, Z_X, walk_time,walk_freq,age,sex,bmi,PC1,PC2,smoke_status,alcohol_status) 
prop.table(table(dat1$X, dat1$Y_yes),1); prop.table(table(dat1$X_Z, dat1$Y_yes),1)	
	aggregate(Y_yes ~ X, dat1, FUN=function(x) { paste( length(x), sum(x), round(sum(x)/length(x),3)) } )
	coef(summary(glm(Y_yes ~ X+Z+ age+sex+smoke_status+alcohol_status +PC1+PC2, family="binomial", data=dat1)))
surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
km.obj <- survfit(surv.obj ~ Z, data=dat1) # KM 是不能带协变量的，Z会被作为分层变量
	survdiff(surv.obj ~ Z, data=dat1) # log-rank test
	plot(km.obj, fun=function(x) 1-x)
	ggsurvplot(km.obj, ylim=c(0,0.08), fun="event") # linetype=rep(1:3,3), palette=rep(c("green","orange","red"),3))
cox.fit <- coxph(surv.obj ~ X+Z+X*Z +walk_time+walk_freq +age+sex+bmi+PC1+PC2+ smoke_status+alcohol_status, data=dat1); print(coef(summary(cox.fit)))
	survminer::ggforest(cox.fit, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1) # 不能显示interaction值
	cox.fit %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
cox.fit <- coxph(surv.ojb ~ X+Z, data = dat1)
#	surv_pred <- survfit(cox.fit, newdata = new_data)
#	surv_summary <- summary(surv_pred, times = 10)
	base_haz <- basehaz(cox.fit, centered = FALSE)
	cum_base_haz <- approx(base_haz$time, base_haz$hazard, xout = 10, method = "linear")$y
	new_data <-  expand.grid(X = levels(dat$X), Z = levels(dat$Z))
	new_data_num <- model.matrix(~ X + Z, data = new_data)[,-1]
	Ht <- exp(new_data_num %*% cox.fit$coefficients) * cum_base_haz
	St <- exp(-Ht)
	var_beta_x <- new_data_num%*%diag(cox.fit$var) 
	se_Ht <- sqrt(var_beta_x) * Ht
	se_St <- St * se_Ht
	z_value <- qnorm(0.975)  # 正态分布的97.5%分位数
	lower_CI <- St * exp(-z_value * se_St)
	upper_CI <- St * exp(z_value * se_St)
res <- data.frame(X = new_data$X, Z = new_data$Z, risk_10y=1-St, LowerCI =1-upper_CI, UpperCI=1-lower_CI)
	ggplot(res, aes(x = Z, y =risk_10y, fill = X)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
	geom_errorbar(aes(ymin =UpperCI, ymax = LowerCI), width = 0.2, position = position_dodge(width = 0.7)) +
	geom_text(aes(label = sprintf("%.2f", UpperCI), y = UpperCI), vjust = -0.5, position = position_dodge(width = 0.7), size = 2) +
	geom_text(aes(label = sprintf("%.2f", LowerCI), y = LowerCI), vjust = 1.5, position = position_dodge(width = 0.7), size = 2) +
	labs(x = "Mendelian mutations", y = "10-year risk", title = "") +
	scale_fill_manual(values=c("green", "yellow", "orange"), name="Walk pace") + theme_minimal() #+ coord_cartesian(ylim = c(0.915, 0.99))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# K-M & COX 科普
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- NULL
	dat$years <- c(2, 3, 4, 6, 9, 12, 15, 18, 20, 22)
	dat$status <- c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1) 
	dat$ages <- seq(85, 40, -5)
	dat$smokes <- rep(1, 1, 1, 1, 0, 0, 0, 1, 1, 1)
surv.obj <- Surv(time=dat$years, event=dat$status, type="right")
	km.obj <- survfit(surv.obj ~ 1, data=dat); round(km.obj$surv,2)
	year <- 10; # 查找距离时间点“10年” 最近的记录点，找到5
	index <- which(km.obj$time <= year)[length(which(km.obj$time <= year))]; index
	1-km.obj$surv[index]; 1-km.obj$upper[index]; 1-km.obj$lower[index]
cox.fit <- coxph(surv.obj ~ ages, data=dat)
	base_hazard <- basehaz(cox.fit, centered = FALSE)
	BETA <- cox.fit$coef
	RR <- exp(sum(BETA * dat$smokes))
	cumulative_risk <- 1 - exp(-base_hazard$hazard * RR) # 计算每个时间点的累积风险