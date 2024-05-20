pacman::p_load(data.table, lubridate, dplyr, tidyr, survival, survminer, randomForest)
inormal <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))

dat0 <- readRDS(file="/work/sph-huangj/data/ukb/Rdata/all.Rdata")
dat <- dat0 %>% mutate(bb_TEST <- runif(nrow(dat0))) %>% filter(ethnic_cat=="White")

Ys <- grep("^icdDate_", names(dat), value=TRUE)
Xs <- grep("^age_sex|age_m|^edu_score|^birth_weight|^height$|^chunk|^leg|^hippo_|^fev1fvc|^stiffness|score_sum$", names(dat), value=TRUE)
Ms <- grep("^bmi$|bb_|bc_", names(dat), value=TRUE)
Zs <- grep("^o$|^se$", names(dat), value=TRUE) # |^rh|shbg|^apoe$|\\.rs
for (Y in Ys) { # 🙍‍
	writeLines(paste('\n\n--> Run:', Y))
	dat1 <- dat %>%
	mutate(
		Y_date=dat[[Y]],
		Y_yes=ifelse( is.na(Y_date), 0,1),
		follow_end_day=fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years=(as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0 )
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	for (X in Xs) { # 🍷
		if (X==Y) next
		print(paste("RUN", X, Y))
		dat1$X=dat1[[X]]
		cox.fit <- coxph(surv.obj ~ X +age+sex +PC1+PC2+weight+smoke_status+alcohol_status, data=dat1)
		print(res.cox <- coef(summary(cox.fit)))
		dat1$M=inormal(dat1$M) # 🔔 这儿更改变量
		fit.X2M <- lm(M ~ X, data=dat1)
		fit.M2Y <- survreg(Surv(follow_years, Y_yes) ~ M + X, data=dat1) # 🐕 这儿用surveg
		res <- mediation::mediate(fit.X2M, fit.M2Y, treat="X", mediator="M")
		logfile=paste(X,M,Y,'.log',sep='.'); sink(logfile) # file.create(logfile)
		print(SUM <- summary(res)) 
		print(paste("RES:", SUM$tau.coef, SUM$z.avg, SUM$n.avg, SUM$d.avg, SUM$d.avg.ci, SUM$d.avg.p)) # 分别表示 Total, ADE, Prop, ACME(B,CI,P) 
		sink()
		fit.rforest <- randomForest(Y_yes ~ X+M+ age+sex+smoke_status+alcohol_status, na.action=na.rm, data=dat1)
		fit.vivi <- vivi(fit=fit.rforest, response="Y_yes", data=dat1)
		print(fit.vivi, digits=1)
		pdf(paste(X,M,Y,'.heatmap.pdf',sep='.')); viviHeatmap(mat=fit.vivi); dev.off()
		pdf(paste(X,M,Y,'.network.log',sep='.')); viviNetwork(mat=fit.vivi); dev.off()
	}
}
