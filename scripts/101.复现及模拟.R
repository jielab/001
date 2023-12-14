pacman::p_load(readxl, tidyverse, TwoSampleMR, MendelianRandomization, RMediation, mediation)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mediation个体数据分析示例
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#示例原文出处：https://rpubs.com/MichaelGilbertUCR/732775
set.seed(12334) 
dat <- iris %>% rename(X=Sepal.Length) %>% 
	mutate(
		random1=runif(nrow(iris),min=min(X),max=max(X)),
		M=X*0.35+random1*0.65,
		random2=runif(nrow(iris),min=min(M),max=max(M)),
	
	Y=M*0.35+random2*.65,
		G1 = ifelse(X > quantile(X, probs=0.90), 2, ifelse(X < quantile(X, probs=0.10), 0, 1))
)
dat %>% dplyr::select(random1, random2, X, M, Y) %>% psych::pairs.panels()
psych::mediate(y="Y", x="X", m="M", data=dat, n.iter=10000) %>% print(short=FALSE)
bda::mediation.test(dat$M, dat$X, dat$Y)
# ACME (Average Causal Mediation Effects) = 0.1132; ACME + ADE = Total effect
fit_X2Y <- lm(Y ~ X, dat)
	res_X2Y <- coef(summary(fit_X2Y)); res_X2Y # beta_X2Y=0.12984，大约0.35 * 0.35。
	beta_X2Y <- res_X2Y[2,1]; se_X2Y <- res_X2Y[2,2]; p_X2Y <- signif(res_X2Y[2,4],2)
fit_X2M <- lm(M ~ X, dat)
	res_X2M <- coef(summary(fit_X2M)); res_X2M # beta_X2M=0.30429, 大约0.35
	beta_X2M <- res_X2M[2,1]; se_X2M <- res_X2M[2,2]; p_X2M <- signif(res_X2M[2,4],2)
fit_M2Y.adjX <- lm(Y ~ M+X, dat) 
	res_M2Y.adjX <- coef(summary(fit_M2Y.adjX)); res_M2Y.adjX # beta_M=0.37194，大约0.35。beta_X=0.01667，是 ADE（Average Direct Effect），不再显著，表明 complete mediation。
	beta_M2Y.adjX <- res_M2Y.adjX[2,1]; se_M2Y.adjX <- res_M2Y.adjX[2,2]; p_M2Y.adjX <- signif(res_M2Y.adjX[2,4],2)
	beta_X2Y.adjM <- res_M2Y.adjX[3,1]; se_X2Y.adjM <- res_M2Y.adjX[3,2]; p_X2Y.adjM <- signif(res_M2Y.adjX[3,4],2)
res <- mediation::mediate(fit_X2M, fit_M2Y.adjX, treat='X', mediator = 'M', boot=T); summary(res)
beta = beta_X2Y - beta_X2Y.adjM; beta # 减法 Difference
beta = beta_X2M * beta_M2Y.adjX; beta # 乘法 Product
CIs = RMediation::medci(beta_X2M, beta_M2Y.adjX, se_X2M, se_M2Y.adjX, type="dop"); CIs
se = CIs$SE; se # Delta 法
se = sqrt(se_X2M^2 + se_M2Y.adjX^2); se # Propagation of errors (POE) 法
se = sqrt(beta_X2M^2 * se_X2M^2 + beta_M2Y.adjX^2 * se_M2Y.adjX^2); se # ?? 法
pval = signif(2*pnorm(abs(beta/se), lower.tail=F),2); pval
print(paste(round(beta_X2Y,3),round(se_X2Y,3),p_X2Y, round(beta_X2M,3),round(se_X2M,3),p_X2M, round(beta_M2Y.adjX,3),round(se_M2Y.adjX,3),p_M2Y.adjX, round(beta,3),round(se,3), pval))


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