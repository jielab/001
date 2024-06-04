pacman::p_load(data.table, lubridate, dplyr, tidyr, survival, survminer, randomForest, randomForestSRC, vivid)
inormal <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
remove_outlier <- function(x) (ifelse((x > (mean(x,na.rm=TRUE) + 3*sd(x,na.rm=TRUE)) | x < (mean(x,na.rm=TRUE) - 3*sd(x,na.rm=TRUE))), NA, x))
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))
expo <- function(x) 1.1^x

dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata")
dat <- dat0 %>% mutate(bb_TEST=runif(nrow(dat0))) %>% dplyr::select(-c("bb_oestradiol","bb_shbg.tp53.rs1042522_C")) %>% filter(ethnic_cat=="White") 
	#rename(F5=vte.F5.rs6025_C, F2=vte.F2.rs1799963_G) % drop_na(???)

Xs <- "?" #grep("^age_sex|age_m|^edu_score|^birth_weight|^height$|^chunk|^leg|^hippo_|^fev1fvc|^stiffness|score_sum$", names(dat), value=TRUE)
Ys <- "?" # grep("^icdDate_", names(dat), value=TRUE)
Zs <- "?" # grep("^bmi$|bb_|bc_", names(dat), value=TRUE) # grep("^o$|^se$", names(dat), value=TRUE) # |^rh|shbg|^apoe$|\\.rs
sink("?.log")

for (Y in Ys) { # 🙍🙍
	writeLines(paste('\n\n--> Run:', Y))
	dat1 <- dat %>% mutate(
		Y.date=dat[[Y]],
		Y.yes=ifelse( is.na(Y.date), 0,1),
		follow_end_day=fifelse(!is.na(Y.date), Y.date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years=(as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25
	) %>% filter( follow_years >0 )

	for (X in Xs) { # 🍷🍷
		if (X==Y) next
		print(paste("RUN", X, Y))
		dat1$X=dat1[[X]]
		if (X=="walk_pace") {
			dat1$X[dat1$X=="steady"] <- NA; dat1$X=droplevels(dat1$X); table(dat1$X) # 对walk_pace 🔔
		}
		dat1 <- dat1 %>% mutate(
			X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
			X_qt = factor(ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")), levels=c("low", "middle", "high"))
		) %>% filter(Y.yes==1 | as.numeric(rownames(dat1)) %%5==0) %>% 
		dplyr::select(Y.yes, follow_years, X, age, sex, bmi, height, smoke_status, alcohol_status, PC1, PC2, grep("^bb_",names(dat1),value=T)) %>%
		na.omit() #%>% slice_sample(n=1000)
		names(dat1) <- gsub("bb_", "", names(dat1))
		surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y.yes)
	
	##	一般的 survival 分析 🔦🔦
		#km.obj <- survfit(surv.obj ~ X_qt, data=dat1) # KM 是不能带协变量的，Z会被作为分层变量
		#	survdiff(surv.obj ~ X_qt, data=dat1) # log-rank test
		#	plot(km.obj, fun=function(x) 1-x)
		#	survminer::ggsurvplot(km.obj, ylim=c(0,0.08), fun="event", pval=TRUE, risk.table=TRUE, risk.table.col="Z", ncensor.plot = TRUE, ggtheme = theme_bw(), palette=c("green","gray","orange"))
		#	survminer::ggforest(fit.cox, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1) # 不能显示interaction值
		#	fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
	
	##	机器学习 survival 分析 🤖🤖
	
	##	众多因素 vivid 交互作用分析 💃💃
		fit.cox <- coxph(surv.obj ~ . - Y.yes - follow_years, data=dat1); print(coef(summary(fit.cox)))
		dat1$Y.yes <- as.factor(dat1$Y.yes); dat1$follow_years <- NULL
		fit.rforest <- randomForest(Y.yes ~ ., na.action=na.omit, data=dat1)
		#fit.rforest <- rfsrc(surv.obj ~ ., data=dat1)
		fit.vivi <- vivi(fit=fit.rforest, response="Y.yes", data=dat1); print(fit.vivi, digits=1)
		pdf(paste(Y,X,'heatmap.pdf',sep='.')); print(viviHeatmap(mat=fit.vivi)); dev.off()
		pdf(paste(Y,X,'network.pdf',sep='.')); print(viviNetwork(mat=fit.vivi)); dev.off()

		#next # 🛑
		for (Z in Zs) { # ⛵ 
			dat1$Z=dat1[[Z]]
		#	dat1 <- dat1 %>% mutate(
		#		Z = factor(ifelse(hardcall(F5)!=2, "F5", ifelse(hardcall(F2) !=2, "F2", "none")), levels=c("none","F2","F5")),
		#		X_Z = factor(paste(Z, X, sep="|"), levels=paste(rep(levels(Z),each=3), rep(levels(X),3), sep="|")),
		#	)
			## 下面计算10年风险并画图 🏮🏮
		#	fit.cox <- coxph(surv.obj ~ X+Z, data=dat1) 
		#	dat2 <- expand.grid(X=levels(dat1$X), Z=levels(dat1$Z))
		#	surv.pred <- survfit(fit.cox, newdata=dat2)
		#	surv.10y <- summary(surv.pred, times=10)
		#	dat2 <- dat2 %>% mutate(risk = 1-t(surv.10y$surv), ci_lower = t(1-surv.10y$upper), ci_upper = t(1-surv.10y$lower))
		#	ggplot(dat2, aes(x=Z, y =risk, fill=X)) + geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7) +
		#		geom_errorbar(aes(ymin =ci_lower, ymax=ci_upper), width=0.2, position=position_dodge(width=0.7)) +
		#		geom_text(aes(label=sprintf("%.2f", risk), y=ci_upper), vjust=-1, position=position_dodge(width=0.7), size=2) +
		#		labs(x="X labels", y="10-year risk (%)", title="") + scale_fill_manual(values=c("green", "gray", "orange"), name="Legend name") + theme_minimal() 
			## 下面进行 mediation 分析 🐎🐎
			fit.X2M <- lm(Z ~ X +age+sex+PC1+PC2, data=dat1) 
			fit.M2Y <- survreg(Surv(follow_years, Y.yes) ~ Z + X +age+sex+PC1+PC2, data=dat1) # 🐕 这儿用surveg
			res <- mediation::mediate(fit.X2M, fit.M2Y, treat="X", mediator="Z")
			print(SUM <- summary(res)) 
			print(paste(c("RES", round(c(SUM$tau.coef, SUM$z.avg, SUM$n.avg, SUM$d.avg, as.numeric(SUM$d.avg.ci)),2), signif(SUM$d.avg.p,2)), collapse=", ")) # 分别表示 Total, ADE, Prop, ACME(B,CI,P) 
		}
	}
}
sink()
