pacman::p_load(data.table, tidyverse, lubridate, survival, randomForest, randomForestSRC, vivid)
inormal <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))

std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
remove_outlier <- function(x) (ifelse((x > (mean(x,na.rm=TRUE) + 3*sd(x,na.rm=TRUE)) | x < (mean(x,na.rm=TRUE) - 3*sd(x,na.rm=TRUE))), NA, x))
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))
expo <- function(x) 1.1^x

run_vivid="NO"; cal_10y_risk="NO"; 
dat0 <- readRDS(file="/work/sph-huangj/data/ukb/Rdata/all.Rdata")
dat <- dat0 %>% filter(ethnic_cat=="White") 
	#rename(F5=vte.F5.rs6025_C, F2=vte.F2.rs1799963_G) % drop_na(???)

Xs <- "?" #grep("^age_sex|age_m|^edu_score|^birth_weight|^height$|^chunk|^leg|^hippo_|^fev1fvc|^stiffness|score_sum$", names(dat), value=TRUE)
Ys <- "?" # grep("^icdDate_", names(dat), value=TRUE)
Ms <- "?" # grep("^bmi$|bb_|bc_", names(dat), value=TRUE) # grep("^o$|^se$", names(dat), value=TRUE) # |^rh|shbg|^apoe$|\\.rs
sink("?.log")

for (Y in Ys) { # 🙍
#	writeLines(paste('\n\n--> Run:', Y))
	dat1 <- dat %>% mutate(
		Y_date=dat[[paste0('icdDate_',Y)]],
		Y_yes=ifelse( is.na(Y_date), 0,1),
		follow_end_day=fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years=(as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25
	) %>% filter( follow_years >0 )

	for (X in Xs) { # 🍷
		if (X==Y) next
		dat1$X=dat1[[X]]
		if (X=="walk_pace") {
			dat1$X[dat1$X=="steady"] <- NA; dat1$X=droplevels(dat1$X); table(dat1$X) # 对walk_pace 🔔
		}
		dat1 <- dat1 %>% mutate(
		#	X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		#	X_qt = factor(ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")), levels=c("low", "middle", "high"))
		) # filter(Y_yes==1 | as.numeric(rownames(dat1)) %%5==0) %>% 
		# names(dat1) <- gsub("bb_", "", names(dat1))
		surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	#	fit.cox <- coxph(surv.obj ~ . - Y_yes - follow_years, data=dat1); print(coef(summary(fit.cox)))
	
	##	一般的 survival 分析 🔦
		#km.obj <- survfit(surv.obj ~ X_qt, data=dat1) # KM 是不能带协变量的，M会被作为分层变量
		#	survdiff(surv.obj ~ X_qt, data=dat1) # log-rank test
		#	plot(km.obj, fun=function(x) 1-x)
		#	survminer::ggsurvplot(km.obj, ylim=c(0,0.08), fun="event", pval=TRUE, risk.table=TRUE, risk.table.col="M", ncensor.plot = TRUE, ggtheme = theme_bw(), palette=c("green","gray","orange"))
		#	survminer::ggforest(fit.cox, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1) # 不能显示interaction值
		#	fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
	
	##	众多因素 vivid 交互作用分析 💃
		if (run_vivid == "YES") {
			dat1$Y_yes <- as.factor(dat1$Y_yes); dat1$follow_years <- NULL
			fit.rforest <- randomForest(Y_yes ~ ., na.action=na.omit, data=dat1)
			#fit.rforest <- rfsrc(surv.obj ~ ., data=dat1)
			fit.vivi <- vivi(fit=fit.rforest, response="Y_yes", data=dat1); print(fit.vivi, digits=1)
			pdf(paste(Y,X,'heatmap.pdf',sep='.')); print(viviHeatmap(mat=fit.vivi)); dev.off()
			pdf(paste(Y,X,'network.pdf',sep='.')); print(viviNetwork(mat=fit.vivi)); dev.off()
		}
		
		#next # 🛑
		for (M in Ms) { # 🐎 
			dat1$M=dat1[[M]]
			if (cal_10y_risk == "YES") {
				dat1 <- dat1 %>% mutate(
				#	M = factor(ifelse(hardcall(F5)!=2, "F5", ifelse(hardcall(F2) !=2, "F2", "none")), levels=c("none","F2","F5")),
				#	XnM = factor(paste(M, X, sep="|"), levels=paste(rep(levels(M),each=3), rep(levels(X),3), sep="|")),
				)
				## 下面计算10年风险并画图 🏮
				fit.cox <- coxph(surv.obj ~ X+M, data=dat1) 
				dat1$risk <- nricens::get.risk.coxph(mdl, 10) # 每人的risk 🔔
				dat2 <- expand.grid(X=levels(dat1$X), M=levels(dat1$M))
				surv.pred <- survfit(fit.cox, newdata=dat2)
				surv.10y <- summary(surv.pred, times=10)
				dat2 <- dat2 %>% mutate(risk = 1-t(surv.10y$surv), ci_lower = t(1-surv.10y$upper), ci_upper = t(1-surv.10y$lower))
				ggplot(dat2, aes(x=M, y =risk, fill=X)) + geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7) +
					geom_errorbar(aes(ymin =ci_lower, ymax=ci_upper), width=0.2, position=position_dodge(width=0.7)) +
					geom_text(aes(label=sprintf("%.2f", risk), y=ci_upper), vjust=-1, position=position_dodge(width=0.7), size=2) +
					labs(x="X labels", y="10-year risk (%)", title="") + scale_fill_manual(values=c("green", "gray", "orange"), name="Legend name") + theme_minimal() 
			}
			## 下面进行 mediation 分析
			fit.X2Y <- survreg(Surv(follow_years, Y_yes) ~ X +age+sex+PC1+PC2, data=dat1)
			fit.X2M <- lm(M ~ X +age+sex+PC1+PC2, data=dat1); res.X2M=coef(summary(fit.X2M))
			fit.M2Y <- survreg(Surv(follow_years, Y_yes) ~ M + X +age+sex+PC1+PC2, data=dat1); res.M2Y=summary(fit.M2Y)$table # 🐕 这儿用surveg
			fit.medi <- mediation::mediate(fit.X2M, fit.M2Y, treat="X", mediator="M"); res <- summary(fit.medi) 
			print(paste("RES: XY|", X,M,Y, nrow(dat1), round(res$tau.coef,3), paste(round(res$tau.ci,3),collapse='-'), signif(res$tau.p,2), "|", 
				nrow(dat1), round(res.X2M[2,1],3), round(res.X2M[2,2],3), signif(res.X2M[2,4],2), "|", 
				nrow(dat1),  round(res.M2Y[2,1],3), round(res.M2Y[2,2],3), signif(res.M2Y[2,4],2), 
					round(res.M2Y[3,1],3), round(res.M2Y[3,2],3), signif(res.M2Y[3,4],2), "|", 
				round(res$n.avg,3), round(res$d.avg,3), paste(round(res$d.avg.ci,3),collapse='-'), signif(res$d.avg.p,2))
			)
		}
	}
}
sink()
