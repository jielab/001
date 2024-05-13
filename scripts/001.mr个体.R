pacman::p_load(tidyverse, TwoSampleMR, MendelianRandomization, mediation, RMediation, psych, bda)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mediation: 从individual数据到summary数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#示例原文出处：https://rpubs.com/MichaelGilbertUCR/732775
set.seed(12334) 
dat <- iris %>% rename(X=Sepal.Length) %>% 
	mutate(
		random1 = runif(nrow(iris), min=min(X), max=max(X)),
		M = X*0.35+random1*0.65,
		random2 = runif(nrow(iris), min=min(M), max=max(M)),
		Y = M*0.35 + random2*.65,
		G1 = ifelse(X > quantile(X, probs=0.80), 2, ifelse(X < quantile(X, probs=0.20), 0, 1)), # a strong SNP
		G2 = ifelse(X > quantile(X, probs=0.98), 2, ifelse(X < quantile(X, probs=0.02), 0, 1))  # a weak SNP		
) %>% dplyr::select(X, M, Y, G1, G2) 
	dat %>% psych::pairs.panels()
	psych::mediate(y="Y", x="X", m="M", data=dat, n.iter=1000) %>% print(short=FALSE)
	bda::mediation.test(dat$M, dat$X, dat$Y)
fit.X2Y <- lm(Y ~ X, data=dat); summary(fit.X2Y)
	coef.X2Y <- coef(summary(fit.X2Y)); coef.X2Y 
	beta.X2Y <- coef.X2Y[2,1]; se.X2Y <- coef.X2Y[2,2]; p.X2Y <- signif(coef.X2Y[2,4],2)
fit.X2M <- lm(M ~ X, data=dat); summary(fit.X2M)
	coef.X2M <- coef(summary(fit.X2M)); coef.X2M
	beta.X2M <- coef.X2M[2,1]; se.X2M <- coef.X2M[2,2]; p.X2M <- signif(coef.X2M[2,4],2)
fit.M2Y <- lm(Y ~ M+X, data=dat); summary(fit.M2Y) # 🐖 这儿必须捎带上X。 X的BETA=0.01667，是ADE，不再显著，表明 complete mediation。
	coef.M2Y <- coef(summary(fit.M2Y)); coef.M2Y 
	beta.M2Y <- coef.M2Y[2,1]; se.M2Y <- coef.M2Y[2,2]; signif(coef.M2Y[2,4],2)
	beta.X2Y.adjM <- coef.M2Y[3,1]; se.X2Y.adjM <- coef.M2Y[3,2]; signif(coef.M2Y[3,4],2)
res <- mediation::mediate(fit.X2M, fit.M2Y, treat="X", mediator="M", boot=T)
	print(SUM <- summary(res))
	print(c(round(SUM$tau.coef,3), round(SUM$z.avg,3), round(SUM$n.avg,3))) # Total, ADE, prop
	print(c(round(SUM$d.avg,3), round(SUM$d.avg.ci,3), signif(SUM$d.avg.p,2))) # ACME
beta.me <- beta.X2M * beta.M2Y; beta.me # 🐕 基于summary数据做乘法Product
# beta.me <- beta.X2Y - beta.X2Y.adjM; beta.me # 🐕 基于summary数据做减法Difference
	se.me <- sqrt(beta.X2M^2 * se.X2M^2 + beta.M2Y^2 * se.M2Y^2); se.me # 🐕 Propagation of errors (POE) 法 
	CIs <- RMediation::medci(beta.X2M, beta.M2Y, se.X2M, se.M2Y, type="dop"); se.me=CIs$SE; se.me # Delta法
	signif(2*pnorm(-abs(beta.me/se.me)),2)	


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MVMR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MRMVInputObject <- mr_mvinput(bx=cbind(ldlc, hdlc, trig), bxse=cbind(ldlcse, hdlcse, trigse), by=chdlodds, byse=chdloddsse); MRMVInputObject
MRMVObject <- mr_mvivw(MRMVInputObject, model="default", correl=FALSE, distribution="normal", alpha=0.05)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR: 从individual数据到summary数据
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(lm(X ~ G1, dat)); summary(lm(X ~ G2, dat)) # G1强，G2弱
summary(lm(M ~ G1, dat)); summary(lm(M ~ G2, dat)) # G2和M没有显著性了
summary(lm(Y ~ G1, dat)); summary(lm(Y ~ G2, dat)) # G1和Y也快没有显著性了
dat$X_byG1 = predict.lm( lm( X ~ G1, data=dat)) 
dat$X_byG2 = predict.lm( lm( X ~ G2, data=dat)) 
summary(lm(Y ~ X, dat));  summary(lm(Y ~ X_byG1, dat)) # 对于一个很强的G1，G1基本就可代表“X本身”
summary(lm(Y ~ X, dat));  summary(lm(Y ~ X_byG2, dat)) # 对于一个很弱的G2，“G2所决定的XX” 跟 ”X本身“，差别很大
	beta_G2X = summary(lm( X ~ G1, data=dat))$coef[2,1]
	beta_G2Y = summary(lm( Y ~ G1, data=dat))$coef[2,1]
	se_G2X = summary(lm( X ~ G1, data=dat))$coef[2,2]
	se_G2Y = summary(lm( Y ~ G1, data=dat))$coef[2,2]
	beta_wald = beta_G2Y / beta_G2X; beta_wald
	se_wald = se_G2Y / beta_G2X; se_wald
	signif(2*pnorm(-abs(beta_wald/se_wald)), 2)