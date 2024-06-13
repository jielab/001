setwd("C:/Users/jiehu/Desktop")
pacman::p_load(dplyr, tidyr, TwoSampleMR, MendelianRandomization, RadialMR, psych, bda)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TwoSampleMR最简单的方法
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat.X.raw <- read.table('D:/data/gwas/ppp/ABO.gz', header=T)
	dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
	dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
	dat.X <- dat.X.raw %>% merge(dat.X.iv, by="SNP") %>% 
	format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', beta_col='BETA', se_col='SE', pval_col='P') 
dat.Y.raw <- read.table('D:/data/gwas/main/y.cad.gz', header=T) 
	dat.Y <- dat.Y.raw %>% merge(dat.X.iv, by="SNP") %>% 
	format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', eaf_col='A1FREQ', beta_col='BETA', se_col='SE', pval_col='P') 
dat <- harmonise_data(dat.X, dat.Y, action=1) 
dat.radial <- format_radial(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome, dat$SNP)
	ivw.radial <- ivw_radial(dat.radial, 0.05/nrow(dat.radial), 3, 0.0001)
	egg.radial <- egger_radial(dat.radial, 0.05/nrow(dat.radial), 3); egg.radial$outliers
	plot_radial(c(ivw.radial, egg.radial), T, F, F)
res <- mr(dat); res; # generate_odds_ratios(res) # 置信区间
	mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_egger_regression_bootstrap"))
	mr_pleiotropy_test(dat); mr_heterogeneity(dat)
	mr_scatter_plot(res, dat)[[1]]
	res.single <- mr_singlesnp(dat)
	mr_forest_plot(res.single)[[1]]
	mr_funnel_plot(res.single)
dat.bsmr <- dat_to_MRInput(dat); dat.bsmr <- dat.bsmr[[1]]
	mr_plot(dat.bsmr, interactive=F)
	mr_forest(dat.bsmr, snp_estimates=F, methods=c("ivw", "median", "wmedian", "egger", "maxlik", "conmix")) + scale_colour_manual(values = c("IVW estimate"="red")) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=15, face="bold"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR：从个体数据到summary数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR的关键词是 genetically determined X。
# 对于一个很强的SNP或PRS，它基本就可代表“X本身”；但对于一个很弱的，它所决定的X，跟X本身差别很大。
library(OneSampleMR)
dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata")
dat <- dat0 %>% filter(ethnic_cat=="White") %>% rename(G=walk_pace.score_sum)
dat1$X.pred = predict.lm( lm( X ~ G, data=dat1)) 
summary(lm(Y ~ X.pred, dat1))
fit1 <- ivreg::ivreg(Y ~ X | G, data = dat1); summary(fit1) # 跟上面的结果一样。 Z可以是 G1+G2+G3
#下面展示基于summary数据的方法，所得结果也是一样
	beta.G2X = summary(lm( X ~ G, data=dat1))$coef[2,1]
	beta.G2Y = summary(lm( Y ~ G, data=dat1))$coef[2,1]
	se.G2X = summary(lm( X ~ G, data=dat1))$coef[2,2]
	se.G2Y = summary(lm( Y ~ G, data=dat1))$coef[2,2]
	beta.wald = beta.G2Y / beta.G2X; beta.wald
	se.wald = se.G2Y / beta.G2X; se.wald
	signif(2*pnorm(-abs(beta.wald/se.wald)), 2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mediation：从个体数据到summary数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#示例原文出处：https://rpubs.com/MichaelGilbertUCR/732775
set.seed(12334) 
dat1 <- iris %>% rename(X=Sepal.Length) %>% 
	mutate(
		X.qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		X.qt = factor(ifelse(X.qt=="q1", "low", ifelse(X.qt=="q5", "high", "middle")), levels=c("low", "middle", "high")),
		random1 = runif(nrow(iris), min=min(X), max=max(X)),
		M = X*0.35+random1*0.65,
		random2 = runif(nrow(iris), min=min(M), max=max(M)),
		Y = M*0.35 + random2*.65,
		Y_yes = c(rep(c(0,0,1),15), rep(0:1,30), rep(c(1,1,0),15)),
		follow_years=X*10 -40,
		G1 = ifelse(X > quantile(X, probs=0.80), 2, ifelse(X < quantile(X, probs=0.20), 0, 1)), # a strong SNP
		G2 = ifelse(X > quantile(X, probs=0.98), 2, ifelse(X < quantile(X, probs=0.02), 0, 1))  # a weak SNP		
) %>% dplyr::select(X, M, Y, Y_yes, follow_years, G1, G2) 
	dat1 %>% psych::pairs.panels()
	psych::mediate(y="Y", x="X", m="M", data=dat1, n.iter=1000) %>% 
		print(short=FALSE) # 🎉 其实，这一行代码就可以了。
	bda::mediation.test(dat1$M, dat1$X, dat1$Y)
fit.X2Y <- lm(Y ~ X, data=dat1); summary(fit.X2Y)
	coef.X2Y <- coef(summary(fit.X2Y)); coef.X2Y 
	beta.X2Y <- coef.X2Y[2,1]; se.X2Y <- coef.X2Y[2,2]; p.X2Y <- signif(coef.X2Y[2,4],2)
fit.X2M <- lm(M ~ X, data=dat1); summary(fit.X2M)
	coef.X2M <- coef(summary(fit.X2M)); coef.X2M
	beta.X2M <- coef.X2M[2,1]; se.X2M <- coef.X2M[2,2]; p.X2M <- signif(coef.X2M[2,4],2)
fit.M2Y <- lm(Y ~ M+X, data=dat1); summary(fit.M2Y) # 🐖 这儿必须捎带上X。 X的BETA=0.01667，是ADE，不再显著，表明 complete mediation。
#fit.M2Y <- survreg(Surv(follow_years, Y_yes) ~ M + X, data=dat1)
	coef.M2Y <- coef(summary(fit.M2Y)); coef.M2Y 
	beta.M2Y <- coef.M2Y[2,1]; se.M2Y <- coef.M2Y[2,2]; signif(coef.M2Y[2,4],2)
	beta.X2Y.adjM <- coef.M2Y[3,1]; se.X2Y.adjM <- coef.M2Y[3,2]; signif(coef.M2Y[3,4],2)
res <- mediation::mediate(fit.X2M, fit.M2Y, treat="X", mediator="M", boot=T); print(SUM <- summary(res))
	print(c(round(SUM$tau.coef,3), round(SUM$z.avg,3), round(SUM$n.avg,3))) # Total, ADE, prop
	print(c(round(SUM$d.avg,3), round(SUM$d.avg.ci,3), signif(SUM$d.avg.p,2))) # ACME
	beta.me <- beta.X2M * beta.M2Y; beta.me; beta.me / beta.X2Y # 🐕 基于summary数据做乘法Product
	# beta.me <- beta.X2Y - beta.X2Y.adjM; beta.me # 🐕 基于summary数据做减法Difference
	se.me <- sqrt(beta.X2M^2 * se.X2M^2 + beta.M2Y^2 * se.M2Y^2); se.me # 🐕 POE法，也称Delta法 
	signif(2*pnorm(-abs(beta.me/se.me)),2)	
