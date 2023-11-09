#install_github("WSpiller/MVMR", build_opts=c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)
pacman::p_load(data.table, dplyr, tidyverse, TwoSampleMR, MendelianRandomization, mediation, MVMR, RMediation)

set.seed(12334)
dat <- iris %>% rename(X=Sepal.Length) # 鸢尾属植物 花萼
dat <- dat %>% 
	mutate(
		random1 = runif(nrow(dat), min=min(X),max=max(X)),
		C = X*0.00 + random1*1.00, # 一个 confounder
		M = X*0.35 + random1*0.65,
		random2 = runif(nrow(dat), min=min(M), max=max(M)),
		Y = M*0.35 + random2*.65,
		G = ifelse(X > quantile(X, probs=0.90), 2, ifelse(X < quantile(X, probs=0.10), 0, 1)),
		G1 = ifelse(X > quantile(X, probs=0.99), 2, ifelse(X < quantile(X, probs=0.01), 0, 1))
	)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 模拟基于individual或summary数据的 MR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(lm(X ~ G, dat)); summary(lm(X ~ G1, dat)) # G强，G1弱
summary(lm(M ~ G, dat)); summary(lm(Y ~ G, dat)) #由于样本量不够大，G和Y没有显著性
dat$X.byG = predict.lm( lm( X ~ G, data=dat)) 
dat$X.byG1 = predict.lm( lm( X ~ G1, data=dat)) 
summary(lm(M ~ X, dat));  summary(lm(M ~ X.byG, dat)) # 对于一个很强的G，G基本就可代表X
summary(lm(M ~ X, dat));  summary(lm(M ~ X.byG1, dat)) # 但对于一个弱的G1，“G1所决定的XX” 跟 “XX”本身，差别很大
summary(lm(Y ~ X.byG, data=dat))
	beta.G2X = summary(lm( X ~ G, data=dat))$coef[2,1]
	beta.G2Y = summary(lm( Y ~ G, data=dat))$coef[2,1]
	se.G2X = summary(lm( X ~ G, data=dat))$coef[2,2]
	se.G2Y = summary(lm( Y ~ G, data=dat))$coef[2,2]
	beta.wald = beta.G2Y / beta.G2X; beta.wald
	se.wald = se.G2Y / beta.G2X; se.wald 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 模拟基于个体数据或2-step的 Mediation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_X2Y <- lm(Y ~ X, dat); summary(fit_X2Y)
	fit_X2M <- lm(M ~ X, dat); summary(fit_X2M)
	fit_M2Y <- lm(Y ~ M, dat); summary(fit_M2Y)
	fit_X2Y.adjM <- lm(Y ~ X+M, dat); summary(fit_X2Y.adjM)
	res_med <- mediation::mediate(fit_X2M, fit_X2Y.adjM, treat='X', mediator='M', boot=T); summary(res_med) # “金标准”
	beta_X2M <- coef(summary(fit_X2M))[2,1]; se_X2M <- coef(summary(fit_X2M))[2,2] 
	beta_M2Y <- coef(summary(fit_M2Y))[2,1]; se_M2Y <- coef(summary(fit_M2Y))[2,2] 
	beta_2step <- beta_X2M * beta_M2Y
	CIs = RMediation::medci(beta_X2M, beta_M2Y, se_X2M, se_M2Y); se_2step = CIs$SE # 不能用 (se_X2M^2 + se_M2Y^2)^0.5
	p_2step <- 2*pnorm(abs(beta_2step/se_2step), lower.tail=F); p_2step	


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 模拟基于基于MVMR的 Mediation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~