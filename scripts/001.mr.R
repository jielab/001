setwd("C:/Users/jiehu/Desktop")
pacman::p_load(dplyr, tidyr, TwoSampleMR, MendelianRandomization, MVMR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 基于individual数据的MR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MR的关键词是 genetically determined X。
# 对于一个很强的SNP或PRS，它基本就可代表“X本身”；但对于一个很弱的，它所决定的X，跟X本身差别很大。
library(OneSampleMR)
dat0 <- readRDS(file="D:/data/ukb/Rdata/all.Rdata")
dat <- dat0 %>% filter(ethnic_cat=="White") %>% rename(G=vte.EUR.score_sum)
dat1$X.pred = predict.lm( lm( X ~ G, data=dat1)) 
summary(lm(Y ~ X.pred, dat1))
fit1 <- ivreg::ivreg(Y ~ X | G, data = dat1); summary(fit1) # 跟上面的结果一样。 Z可以是 G1+G2+G3
#下面展示基于summary数据的方法，所得结果也是一样
	beta_G2X = summary(lm( X ~ G, data=dat1))$coef[2,1]
	beta_G2Y = summary(lm( Y ~ G, data=dat1))$coef[2,1]
	se_G2X = summary(lm( X ~ G, data=dat1))$coef[2,2]
	se_G2Y = summary(lm( Y ~ G, data=dat1))$coef[2,2]
	beta_wald = beta_G2Y / beta_G2X; beta_wald
	se_wald = se_G2Y / beta_G2X; se_wald
	signif(2*pnorm(-abs(beta_wald/se_wald)), 2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# TwoSampleMR最简单的方法
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iv.snp <- read.table('D:/data/gwas/main/walk_pace.top.snp', header=T)
dat_X <- read.table('D:/data/gwas/main/walk_pace.gz', header=T) %>% merge(iv.snp, by="SNP") %>% 
	format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') 
dat_Y <- read.table('D:/data/gwas/main/y.vte.gz', header=T) %>%
	merge(iv.snp, by="SNP") %>% 
	format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') 
dat <- harmonise_data(dat_X, dat_Y, action=1) 
res <- mr(dat); res; generate_odds_ratios(res) #置信区间
	res_single <- mr_singlesnp(dat)
	mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_egger_regression_bootstrap"))
mr_pleiotropy_test(dat); mr_heterogeneity(dat)
mr_scatter_plot(res, dat)[[1]]
	mr_forest_plot(res_single)[[1]]
	mr_funnel_plot(res_single)
bsmr_dat <- dat_to_MRInput(dat); bsmr_dat <- bsmr_dat[[1]]
	mr_plot(bsmr_dat, interactive=F)
	mr_forest(bsmr_dat, snp_estimates=F, methods=c("ivw", "median", "wmedian", "egger", "maxlik", "conmix")) + scale_colour_manual(values = c("IVW estimate"="red")) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=15, face="bold"))