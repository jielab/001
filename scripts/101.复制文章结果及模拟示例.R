pacman::p_load(readxl, tidyverse, TwoSampleMR, MendelianRandomization, mediation)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 重现阜外医院NAFLD文章结果
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_X <- as.data.frame(read_excel("D:/Downloads/2023阜外医院NAFLD附件.xlsx", sheet="Table S3", skip=1)) %>% rename(EA="Effect Allele", EAF="Effect Allele Frequency", P="P value")
dat_Y <- extract_outcome_data(dat_X$SNP, "ebi-a-GCST90091033")
#dat <- harmonise_data(dat_X, dat_Y); dat %>% mr()
dat0 <- merge(dat_X, dat_Y, by="SNP") 
	dat0$beta.outcome.NEW <- ifelse(dat0$EA==dat0$effect_allele.outcome, dat0$beta.outcome, 0-dat0$beta.outcome)
	mr_ivw(mr_input(dat0$Beta, dat0$SE, dat0$beta.outcome, dat0$se.outcome))


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
		Y=M*0.35+random2*.65
)
fit_X2Y <- lm(Y ~ X, dat)
	res_X2Y <- coef(summary(fit_X2Y)) # beta_X2Y=0.12984，大约0.35 * 0.35。
	beta_X2Y <- res_X2Y[2,1]; se_X2Y <- res_X2Y[2,2]; p_X2Y <- signif(res_X2Y[2,4],2)
fit_X2M <- lm(M ~ X, dat)
	res_X2M <- coef(summary(fit_X2M)) # beta_X2M=0.30429, 大约0.35
	beta_X2M <- res_X2M[2,1]; se_X2M <- res_X2M[2,2]; p_X2M <- signif(res_X2M[2,4],2)
fit_M2Y.adjX <- lm(Y ~ M+X, dat)
	res_M2Y.adjX <- coef(summary(fit_M2Y.adjX)) # beta_M=0.37194，大约0.35。beta_X=0.01667，是 ADE（Average Direct Effect），不再显著，表明 complete mediation。
	beta_M2Y.adjX <- res_M2Y.adjX[2,1]; se_M2Y.adjX <- res_M2Y.adjX[2,2]; p_M2Y.adjX <- signif(res_M2Y.adjX[2,4],2)
beta_2step <- round(beta_X2M * beta_M2Y.adjX, 4) # beta_2step=0.1132，是 ACME（Average Causal Mediation Effects）。ACME + ADE = Total effect
	se=(se_X2M^2 + se_M2Y.adjX^2)^0.5 #有的地方建议用 CIs = RMediation::medci(beta_X2M, beta_M2Y.adjX, se_X2M, se_M2Y.adjX, type='MC'); se_2step = CIs$SE
	p_2step <- signif(2*pnorm(abs(beta_2step/se_2step), lower.tail=F),2) # 目前有点对不上？？
print(paste(round(beta_X2Y,3),round(se_X2Y,3),p_X2Y, round(beta_X2M,3),round(se_X2M,3),p_X2M, round(beta_M2Y.adjX,3),round(se_M2Y.adjX,3),p_M2Y.adjX, round(beta_2step,3),round(se_2step,3),p_2step))
res = mediation::mediate(fit_X2M, fit_M2Y.adjX, treat='X', mediator = 'M', boot=T); summary(res)
	se=(upper limit – lower limit) / 3.92
	confint(res)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 更多模拟数据分析展示
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(12334)
round(cor(iris[1:4]), 2) # 鸢尾属植物，花萼Sepal的Length和Width之间基本没有关联
dat <- iris %>% rename(X=Sepal.Length) %>%
	mutate(
		ran1 = runif(nrow(iris), min=min(X),max=max(X)),
		ran2 = runif(nrow(iris), min=min(X),max=max(X)),
		X1 = X*0.5 + ran1*0.5,
		X2 = X*0.5 + ran2*0.5,
		M = X1*0.3 + ran1*0.7,
		Y = M*0.3 + ran2*0.7,
		G1 = ifelse(X1 > quantile(X1, probs=0.90), 2, ifelse(X1 < quantile(X1, probs=0.10), 0, 1)),
		G2 = ifelse(X2 > quantile(X2, probs=0.90), 2, ifelse(X2 < quantile(X2, probs=0.10), 0, 1))
	)
round(cor(subset(dat, select=c(ran1, ran2, X, X1, X2, M, Y, G1, G2))), 2)
fit.X2Y <- lm(Y ~ X1+X2, data=dat); summary(fit.X2Y)	
ggplot(mtcars, aes(x=wt, y=mpg, color=cyl, shape=cyl)) + geom_point() + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)
 

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