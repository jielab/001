pacman::p_load(data.table, tidyverse, lubridate, survival, randomForest, randomForestSRC, vivid, rcssci)
inormal <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
remove_outlier <- function(x) (ifelse((x > (mean(x,na.rm=TRUE) + 3*sd(x,na.rm=TRUE)) | x < (mean(x,na.rm=TRUE) - 3*sd(x,na.rm=TRUE))), NA, x))
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))
expo <- function(x) 1.1^x
rb <- function(x) (round(x,3)); rp <- function(x) (signif(x,2))

dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata")
dat <- dat0 %>% filter(ethnic_cat=="White") %>% rename(F5=vte.F5.rs6025_C, F2=vte.F2.rs1799963_G) # %>% drop_na(???)
dat$abo.a = ifelse(dat$blood_group =="OO", 0, ifelse(dat$blood_group =="AO", 1, ifelse(dat$blood_group =="AA", 2, NA)))
dat$abo.ao = ifelse(dat$abo =="O", 0, ifelse(dat$abo =="A", 1, NA))
dat$blood_group <- factor(dat$blood_group, levels=c("OO","AO","AA","BO","BB","AB"))
dat$fut2 <- factor(dat$fut2.rs601338_A)
dat$rh <- ifelse(dat$rhd.rs590787_A >0.5, 1, 0)

X <- "age_fsex" #grep("^age_sex|age_m|^edu_score|^birth_weight|^height$|^chunk|^leg|^hippo_|^fev1fvc|^stiffness|score_sum$", names(dat), value=TRUE)
Y <- "bald" # grep("^icdDate_", names(dat), value=TRUE) 去掉 icdDate_ 前缀
M <- "bb_GLU" # grep("^bmi$|bb_|bc_", names(dat), value=TRUE) # grep("^o$|^se$", names(dat), value=TRUE) # |^rh|shbg|^apoe$|\\.rs

for (Y in Ys) { # 🙍
#	writeLines(paste('\n\n--> Run:', Y))
	dat1 <- dat %>% mutate(
		Y_date=dat[[paste0('icdDate_',Y)]],
		Y_yes=ifelse( is.na(Y_date), 0,1),
		follow_end_day=fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years=(as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25
	) %>% filter( follow_years >0 )

	dat1$X=dat1[[X]]
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	fit.X2Y <- coxph(surv.obj ~ X + fut2 + X*fut2 + rh + age+sex+bmi+smoke_status+alcohol_status +PC1+PC2+PC3+PC4, data=dat1); print(coef(summary(fit.X2Y)))	
	
	##	一般的 survival 分析 🔦
	km.obj <- survfit(surv.obj ~ X_qt, data=dat1) # KM 是不能带协变量的，M会被作为分层变量
		survdiff(surv.obj ~ X_qt, data=dat1) # log-rank test
		plot(km.obj, fun=function(x) 1-x)
		survminer::ggsurvplot(km.obj, ylim=c(0,0.08), fun="event", pval=TRUE, risk.table=TRUE, risk.table.col="M", ncensor.plot = TRUE, ggtheme = theme_bw(), palette=c("green","gray","orange"))
		survminer::ggforest(fit.cox, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1)
		fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
		rcssci_cox(data=dat, time="time", y="status", x="sbp", covs=c("age","gender"),  prob=0.1, filepath= 'D:') # 🏮
	
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
				## 10年风险 🏮
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
			## mediation 分析。 mediate函数只能实现survreg 拟合的参数化生存回归模型，无法实现coxph拟合的半参数的生存回归模型
			set.seed(12345)
			fit.X2Y <- survreg(Surv(follow_years, Y_yes) ~ X +age+sex+PC1+PC2, data=dat1); res.X2Y=summary(fit.X2Y)$table
			fit.X2M <- lm(M ~ X +age+sex+PC1+PC2, data=dat1); res.X2M=coef(summary(fit.X2M))
			fit.M2Y <- survreg(Surv(follow_years, Y_yes) ~ M + X +age+sex+PC1+PC2, data=dat1); res.M2Y=summary(fit.M2Y)$table # 🐕 这儿用surveg
			fit.medi <- mediation::mediate(fit.X2M, fit.M2Y, treat="X", mediator="M", boot=F, outcome="status", sims=1000); res <- summary(fit.medi) 
			print(paste("RES:", X, M, Y, nrow(dat1), rb(res.X2Y[2,1]), rp(res.X2Y[2,4]), 
				"|X2M", rb(res.X2M[2,1]), rp(res.X2M[2,4]), "|M2Y", rb(res.M2Y[2,1]), rp(res.M2Y[2,4]), 
				"|Me", rb(res$tau.coef), rp(res$tau.p), rb(res$d.avg), rp(res$d.avg.p), rb(res$z.avg), rp(res$z.avg.p), rb(res$n.avg), rp(res$n.avg.p) # TOT, ACME, ADE, Prop
			))
		}
	}
}
