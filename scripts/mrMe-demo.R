setwd("C:/Users/jiehu/Desktop")
# devtools::install_github("scllin/mrMed")
pacman::p_load(data.table, stringi, dplyr, tidyverse, OneSampleMR, TwoSampleMR, MendelianRandomization, mrMed, RadialMR, psych, bda)

pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')
rb <- function(x) (round(x,3)); rp <- function(x) (signif(x,2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR：从"一"开始 (OneSampleMR)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR的关键词是 genetically determined X。
# 对于一个很强的SNP或PRS，它基本就可代表“X本身”；但对于一个很弱的SNP，它所决定的X，跟X本身差别很大。
dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata")
	dat <- dat0 %>% filter(ethnic_cat=="White") %>% rename(G=walk_pace.score_sum)
	dat1$X.pred = predict.lm( lm( X ~ G, data=dat1)) 
	summary(lm(Y ~ X.pred, dat1))
fit1 <- ivreg::ivreg(Y ~ X | G, data = dat1); summary(fit1) # 跟上面的结果一样，Z可以是 G1+G2+G3
# 从个体数据到summary数据的飞跃 🚀
	beta.G2X = summary(lm( X ~ G, data=dat1))$coef[2,1]
	beta.G2Y = summary(lm( Y ~ G, data=dat1))$coef[2,1]
	se.G2X = summary(lm( X ~ G, data=dat1))$coef[2,2]
	se.G2Y = summary(lm( Y ~ G, data=dat1))$coef[2,2]
	beta.wald = beta.G2Y / beta.G2X; beta.wald
	se.wald = se.G2Y / beta.G2X; se.wald
	rp(2*pnorm(-abs(beta.wald/se.wald)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR: 因"二"流行 (TwoSampleMR)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iv.file <- 'D:/data/gwas/main/ppp_ABO.top.snp'
dat.X.file <- 'D:/data/gwas/main/ppp_ABO.gz'
dat.Y.file <- 'D:/data/gwas/main/covid_B2.EUR.gz'
dat.X.raw <- read.table(dat.X.file, header=T)
	names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	if (file.exists(iv.file)) {
		dat.X.iv <- read.table(iv.file, header=T); names(dat.X.iv) <- "SNP" 
	} else {
		dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
		dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
	}
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
set.seed(12334) 
dat1 <- iris %>% rename(X=Sepal.Length) %>% 
	mutate(
		random1 = runif(nrow(iris), min=min(X), max=max(X)),
		M = X*0.35+random1*0.65,
		random2 = runif(nrow(iris), min=min(M), max=max(M)),
		Y = M*0.35 + random2*.65,
		Y_yes = c(rep(c(0,0,1),15), rep(0:1,30), rep(c(1,1,0),15)),
		G1 = ifelse(X > quantile(X, probs=0.80), 2, ifelse(X < quantile(X, probs=0.20), 0, 1)) # a strong SNP
) %>% dplyr::select(X, M, Y, Y_yes, G1, G2) 
	dat1 %>% psych::pairs.panels()
	psych::mediate(y="Y", x="X", m="M", data=dat1, n.iter=1000) %>% print(short=FALSE) # 🎉 其实，这一行代码就可以了。
	bda::mediation.test(dat1$M, dat1$X, dat1$Y)
fit.X2Y <- lm(Y ~ X, data=dat1); summary(fit.X2Y)
fit.X2M <- lm(M ~ X, data=dat1); summary(fit.X2M)
fit.M2Y <- lm(Y ~ M+X, data=dat1); summary(fit.M2Y) # 🐖 这儿必须捎带上X
res <- mediation::mediate(fit.X2M, fit.M2Y, treat="X", mediator="M", boot=T); print(SUM <- summary(res))
	c(rb(SUM$tau.coef), rb(SUM$d.avg), rb(SUM$z.avg), rb(SUM$n.avg)) # Total, ACME, ADE, Prop
# 从个体数据到summary数据的飞跃 🚀
	coef.X2Y <- coef(summary(fit.X2Y)); coef.X2Y 
		beta.X2Y <- coef.X2Y[2,1]; se.X2Y <- coef.X2Y[2,2]; p.X2Y <- rp(coef.X2Y[2,4])
	coef.X2M <- coef(summary(fit.X2M)); coef.X2M
		beta.X2M <- coef.X2M[2,1]; se.X2M <- coef.X2M[2,2]; p.X2M <- rp(coef.X2M[2,4])
	coef.M2Y <- coef(summary(fit.M2Y)); coef.M2Y 
		beta.M2Y <- coef.M2Y[2,1]; se.M2Y <- coef.M2Y[2,2]; rp(coef.M2Y[2,4])
		beta.X2Y.dr <- coef.M2Y[3,1]; se.X2Y.dr <- coef.M2Y[3,2]; rp(coef.M2Y[3,4])
	print(beta.me <- beta.X2M * beta.M2Y) # 乘法
		print(se.me <- sqrt(beta.X2M^2 * se.X2M^2 + beta.M2Y^2 * se.M2Y^2)
		print(p.me <- rp(2*pnorm(-abs(beta.me/se.me))))
	print(beta.dr <- beta.X2Y - beta.me) # 减法
		print(se.dr <- sqrt(se.X2Y^2 + se.me^2))
		print(p.dr <- rp(2*pnorm(-abs(beta.dr/se.dr))))
	print(beta.prop <- beta.me / beta.X2)
		print(se.prop <- sqrt(se.me^2/beta.X2Y^2 + beta.me^2*se.X2Y^2/beta.X2Y^4))
		print(p.prop <- rp(2*pnorm(-abs(beta.prop/se.prop))))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 实战 ABO-ALP-COVID
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 参考https://rdrr.io/github/scllin/toypack/src/R/mrMed.R
dir.Y = 'D:/data/gwas/main/clean'; Y = 'covid_C2' 
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.EUR.gz'), header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
dir.X = dir.Y; X = 'ABO'
	dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
	names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
	dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
	dat.X <- dat.X.raw %>% merge(dat.X.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
	dat.Y <- dat.Y.raw %>% merge(dat.X.iv, by="SNP") %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
# Total effect
	dat.XY <- harmonise_data(dat.X, dat.Y, action=1) 
	fit.X2Y <- mr(dat.XY, method_list=c("mr_wald_ratio", "mr_ivw")); fit.X2Y 
	beta.X2Y <- fit.X2Y$b; se.X2Y <- fit.X2Y$se; p.X2Y <- rp(fit.X2Y$pval)
dir.M = 'D:/data/gwas/bb/clean'; M = 'bb_ALP'
	dat.M.raw <- read.table(paste0(dir.M, '/', M, '.gz'), header=T)
	names(dat.M.raw) <- stri_replace_all_regex(toupper(names(dat.M.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	if (file.exists(paste0(dir.M, '/', M, '.top.snp'))) {
		dat.M.iv <- read.table(paste0(dir.M, '/', M, '.top.snp'), header=T); names(dat.M.iv) <- "SNP" 
	} else {
		dat.M.sig <- dat.M.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
		dat.M.iv <- dat.M.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
	}
# Step1: X --> M
	dat.M.4X <- dat.M.raw %>% merge(dat.X.iv) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
	dat.XM <- harmonise_data(dat.X, dat.M.4X, action=1)
	fit.X2M <- mr(dat.XM, method_list=c("mr_wald_ratio", "mr_ivw")); fit.X2M 
	beta.X2M <- fit.X2M$b; se.X2M <- fit.X2M$se; p.X2M <- rp(fit.X2M$pval)
# 合并 dat.X 和 dat.M 和 dat.Y
	dat.XnM.iv <- rbind(dat.X.iv, dat.M.iv) %>% unique() # 最理想的是把 dat.X和dat.M中所有的显著性snp合并然后 %>% clump_data() %>% select(SNP)
	dat.X.mv <- dat.X.raw %>% merge(dat.XnM.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
	dat.M.mv <- dat.M.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =M)
	dat.Y.mv <- dat.Y.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
	dat.XY.mv <- harmonise_data(dat.X.mv, dat.Y.mv, action=1)
	dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action=1); names(dat.XM.mv) <- gsub("outcome", "mediator", names(dat.XM.mv))
	dat <- merge(dat.XM.mv, dat.XY.mv, by="SNP")
	bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
	names(dat) <- gsub("\\.x$", "", names(dat)); dat <- subset(dat, select=!grepl("\\.y", names(dat)))
# Step 2: M (+X) --> Y 
	dat.M2Y <- MendelianRandomization::mr_mvinput(bx=cbind(dat$beta.mediator, dat$beta.exposure), bxse=cbind(dat$se.mediator, dat$se.exposure), by=dat$beta.outcome, byse=dat$se.outcome)
	fit.M2Y <- MendelianRandomization::mr_mvivw(dat.M2Y, model="default", correl=FALSE, distribution="normal", alpha=0.05); fit.M2Y
	beta.M2Y <- fit.M2Y$Estimate[["Bx1"]]; se.M2Y <- fit.M2Y$StdError[["Bx1"]] ; p.M2Y <- fit.M2Y$Pvalue[["Bx1"]]
	beta.X2Y.dr <- fit.M2Y$Estimate[["Bx2"]]; se.X2Y.dr <- fit.M2Y$StdError[["Bx2"]] ; p.X2Y.dr <- fit.M2Y$Pvalue[["Bx2"]]
# 自创 Mediation 方法
	beta.me <- beta.X2M * beta.M2Y # 🐕 乘法
	se.me <- sqrt(beta.X2M^2 * se.X2M^2 + beta.M2Y^2 * se.M2Y^2); p.me <- 2*pnorm(-abs(beta.me/se.me))
	beta.dr <- beta.X2Y - beta.me # 减法
	se.dr <- sqrt(se.X2Y^2 + se.me^2); p.dr <- 2*pnorm(-abs(beta.dr/se.dr))
	beta.prop <- beta.me / beta.X2Y
	se.prop <- sqrt(se.me^2/beta.X2Y^2 + beta.me^2*se.X2Y^2/beta.X2Y^4); p.prop <- 2*pnorm(-abs(beta.prop/se.prop))
# mrMed 方法
	dat1 <- dat
	names(dat1) <- stri_replace_all_regex(names(dat1), pattern=c("exposure", "mediator", "outcome"), replacement=c("X", "M", "Y"), vectorize_all=FALSE)
	dat1 <- dat1 %>% mutate (
		Gx = ifelse(pval.X >5e-08 | pval.M >5e-08, 0, 1), # 仅仅用于 *_0 的 function 里面
		Gx_plum = Gx, # (Gx==1 & [XY]keep=TRUE), 用于 X -> M或Y
		Gm = ifelse(pval.M <5e-08, 1, 0), # (... & [XM]keep=TRUE), 没用 
		Gm_plum = Gm, # (... & [MY]keep=TRUE), 用于 M -> Y
		G_mvmr = ifelse(Gx_plum==0 & Gm_plum==0, 0, 1)
	)
	res <- mrMed(dat_mrMed=dat1, method_list="Prod_IVW") # 还有 Diff_IVW, Prod_Median
# 汇总分析结果
	X_str=paste0(X, "(", nrow(dat.X.iv), ")"); M_str=paste0(M, "(", nrow(dat.M.iv), ")");
	print(paste("RES:", X_str, M_str, Y, "|[X2Y]", nrow(dat.XY), rb(beta.X2Y), rp(p.X2Y), 
		"[X2M]", nrow(dat.XM), rb(beta.X2M), rp(p.X2M), rb(beta.dr), rp(p.dr),
		"[M2Y]", nrow(dat), rb(beta.M2Y), rp(p.M2Y), 
		"[me]", rb(beta.me), rp(p.me), rb(beta.X2Y.dr), rp(p.X2Y.dr), rb(beta.prop), rp(p.prop), # ACME, ADE, Prop(rho) 
		"[mrMed]", rb(res$TE$b), rp(res$TE$p), rb(res$IE$b), rp(res$IE$p), rb(res$DE$b), rp(res$DE$p), rb(res$rho$b), rp(res$rho$p) # TE, ACME, ADE, Prop(rho) 
	))

