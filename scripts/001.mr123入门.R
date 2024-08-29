pacman::p_load(data.table, stringi, tidyverse, OneSampleMR, TwoSampleMR, MendelianRandomization, mrMed, RadialMR, psych, bda)
pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')
rb <- function(x) (round(x,3)); rp <- function(x) (signif(x,2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR：从"一"开始 (OneSampleMR)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR的关键词是 genetically determined X。
# 对于一个很强的SNP或PRS，它基本就可代表“X本身”；但对于一个很弱的SNP，它所决定的X，跟X本身差别很大。
dat <- readRDS(file="D:/data/ukb/Rdata/all.Rdata") %>% filter(ethnic_cat=="White")
	dat <- dat %>% rename(X=bmi, Y=icdCt_cad, G=bmi.EUR941.score_sum)
	dat$X.pred = predict.lm( lm( X ~ G, data=dat, na.action=na.exclude)) # 🏮
	summary(lm(Y ~ X.pred, data=dat, na.action=na.exclude))
fit.mr <- ivreg::ivreg(Y ~ X | G, data = dat) # 可以是 G1+G2+G3
	coef(summary(fit.mr))  # 跟上面的结果一样
# 从个体数据到summary数据的飞跃 🚀
	beta.G2X = summary(lm( X ~ G, data=dat))$coef[2,1]
	beta.G2Y = summary(lm( Y ~ G, data=dat))$coef[2,1]
	se.G2X = summary(lm( X ~ G, data=dat))$coef[2,2]
	se.G2Y = summary(lm( Y ~ G, data=dat))$coef[2,2]
	beta.wald = beta.G2Y / beta.G2X; beta.wald
	se.wald = se.G2Y / beta.G2X; se.wald
	rp(2*pnorm(-abs(beta.wald/se.wald)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR: 因"二"流行 (TwoSampleMR)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat.X.file <- 'D:/data/gwas/main/ppp_ABO.gz'
dat.Y.file <- 'D:/data/gwas/main/covid_B2.EUR.gz'
dat.X.raw <- read.table(dat.X.file, header=T)
	names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
	dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
	dat.X <- dat.X.raw %>% merge(dat.X.iv, by="SNP") %>% 
	format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') 
dat.Y.raw <- read.table(dat.Y.file, header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
 	dat.Y <- dat.Y.raw %>% merge(dat.X.iv, by="SNP") %>% 
	format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') 
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
	mr_plot(mr_allmethods( dat.bsmr, method="main", iterations = 50)) # method="all"
	mr_forest(dat.bsmr, snp_estimates=F, methods=c("ivw", "median", "wmedian", "egger", "maxlik", "conmix")) 
	+ scale_colour_manual(values = c("IVW estimate"="red")) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=15, face="bold"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR："三"角关系复杂 (MVMR)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 示例原文出处：https://rpubs.com/MichaelGilbertUCR/732775
set.seed(12334) # 🏮
dat <- iris %>% rename(X=Sepal.Length) %>% 
	mutate(
		random1 = runif(nrow(iris), min=min(X), max=max(X)),
		M = X*0.35+random1*0.65,
		random2 = runif(nrow(iris), min=min(M), max=max(M)),
		Y = M*0.35 + random2*.65,
		Y_yes = c(rep(c(0,0,1),15), rep(0:1,30), rep(c(1,1,0),15)),
		G = ifelse(X > quantile(X, probs=0.80), 2, ifelse(X < quantile(X, probs=0.20), 0, 1)) # a strong SNP
) %>% dplyr::select(X, M, Y, Y_yes, G) 
	dat %>% psych::pairs.panels()
	psych::mediate(y="Y", x="X", m="M", data=dat, n.iter=1000) %>% print(short=FALSE) # 🎉 其实，这一行代码就可以了。
# 基于个体数据的mediation
fit.X2Y <- lm(Y ~ X, data=dat); summary(fit.X2Y)
fit.X2M <- lm(M ~ X, data=dat); summary(fit.X2M)
fit.M2Y <- lm(Y ~ M+X, data=dat); summary(fit.M2Y) # 🐖 这儿必须捎带上X
res <- mediation::mediate(fit.X2M, fit.M2Y, treat="X", mediator="M", boot=T); print(SUM <- summary(res))
	c(rb(SUM$d.avg), rb(SUM$z.avg), rb(SUM$tau.coef), rb(SUM$n.avg)) # ACME, ADE, Total, Prop
# 从个体数据到summary数据的飞跃 🚀，但是这里没有用到SNP作为IV
coef.X2Y <- coef(summary(fit.X2Y)); coef.X2Y
	beta.X2Y <- coef.X2Y[2,1]; se.X2Y <- coef.X2Y[2,2]; p.X2Y <- rp(coef.X2Y[2,4])
coef.X2M <- coef(summary(fit.X2M)); coef.X2M
	beta.X2M <- coef.X2M[2,1]; se.X2M <- coef.X2M[2,2]; p.X2M <- rp(coef.X2M[2,4])
coef.M2Y <- coef(summary(fit.M2Y)); coef.M2Y 
	beta.M2Y <- coef.M2Y[2,1]; se.M2Y <- coef.M2Y[2,2]; rp(coef.M2Y[2,4])
	beta.X2Y.dr <- coef.M2Y[3,1]; se.X2Y.dr <- coef.M2Y[3,2]; rp(coef.M2Y[3,4])
print(beta.me <- beta.X2M * beta.M2Y) # 乘法
	print(se.me <- sqrt(beta.X2M^2 * se.X2M^2 + beta.M2Y^2 * se.M2Y^2)) # SE算法不太对
	print(p.me <- rp(2*pnorm(-abs(beta.me/se.me))))
print(beta.dr <- beta.X2Y - beta.me) # 减法
	print(se.dr <- sqrt(se.X2Y^2 + se.me^2))
	print(p.dr <- rp(2*pnorm(-abs(beta.dr/se.dr))))
print(beta.prop <- beta.me / beta.X2)
	print(se.prop <- sqrt(se.me^2/beta.X2Y^2 + beta.me^2*se.X2Y^2/beta.X2Y^4))
	print(p.prop <- rp(2*pnorm(-abs(beta.prop/se.prop))))