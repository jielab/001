pacman::p_load(tidyverse, TwoSampleMR, MendelianRandomization, psych, RMediation)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 模拟基于individual或summary数据的 MR
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
		G1 = ifelse(X > quantile(X, probs=0.80), 2, ifelse(X < quantile(X, probs=0.20), 0, 1)), # a strong SNP
		G2 = ifelse(X > quantile(X, probs=0.98), 2, ifelse(X < quantile(X, probs=0.02), 0, 1))  # a weak SNP		
)
dat %>% dplyr::select(random1, random2, X, M, Y) %>% psych::pairs.panels()
psych::mediate(y="Y", x="X", m="M", data=dat, n.iter=10000) %>% print(short=FALSE)
bda::mediation.test(dat$M, dat$X, dat$Y)
fit.dv<-glm(dementia_status~CRP+BC_ra+age+sex+edu_level+eth_white+TDI+smoking+alcohol+TPA+TV+fruit+vegetable+oilyfish+processed_meat+red_meat,family=binomial(link=logit),merge_B_dementia)
	fit.mediator<-lm(CRP~BC_ra+age+sex+edu_level+eth_white+TDI+smoking+alcohol+TPA+TV+fruit+vegetable+oilyfish+processed_meat+red_meat,merge_B_dementia)
	results<-mediate(fit.mediator,fit.dv,treat="BC_ra",mediator="CRP",sims=1000)
	SUM<-summary(results)
	ACME<-c(SUM$d.avg,SUM$d.avg.ci,SUM$d.avg.p)
	ADE<-c(SUM$z.avg,SUM$z.avg.ci,SUM$z.avg.p)
	Prop<-c(SUM$n.avg,SUM$n.avg.ci,SUM$n.avg.p)
	Total<-c(SUM$tau.coef,SUM$tau.ci,SUM$tau.p)
	res<-rbind(ACME,ADE,Prop,Total); res
	results<-rbind(coef(summary(fit.mediator))[2,], coef(summary(fit.dv))[2,], res)
	colnames(results)<-c("Estimate","SE/95L","Z/95U","P")

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
# 1. 乘法 Product
beta_medi = beta_X2M * beta_M2Y.adjX; beta_medi # 也有人简单的 beta_X2M * beta_M2Y，不考虑 .adjX
se_medi = sqrt(beta_X2M^2 * se_X2M^2 + beta_M2Y.adjX^2 * se_M2Y.adjX^2); se_medi # Delta 法，可通过R包实现 CIs = RMediation::medci(beta_X2M, beta_M2Y.adjX, se_X2M, se_M2Y.adjX, type="dop"); CIs; se_medi = CIs$SE; se_medi 
# 2. 减法 Difference
beta_medi = beta_X2Y - beta_X2Y.adjM; beta_medi 
se_medi = sqrt(se_X2Y^2 + se_X2Y.adjM^2); se_medi # Propagation of errors (POE) 法 
pval_medi = signif(2*pnorm(-abs(beta_medi/se_medi)),2); pval_medi
print(paste(round(beta_X2Y,3),round(se_X2Y,3),p_X2Y, round(beta_X2M,3),round(se_X2M,3),p_X2M, round(beta_M2Y.adjX,3),round(se_M2Y.adjX,3),p_M2Y.adjX, round(beta_medi,3),round(se_medi,3), pval_medi))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MVMR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MRMVInputObject <- mr_mvinput(bx=cbind(ldlc, hdlc, trig), bxse=cbind(ldlcse, hdlcse, trigse), by=chdlodds, byse=chdloddsse); MRMVInputObject
MRMVObject <- mr_mvivw(MRMVInputObject, model="default", correl=FALSE, distribution="normal", alpha=0.05)
